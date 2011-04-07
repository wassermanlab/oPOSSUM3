package oPossumGeneSSAWeb;

use base 'CGI::Application';

use oPossumWebInclude;
use oPossumGeneWebInclude;

use Data::Dumper;    # for debugging only
use File::Temp qw/ tempdir /;

use OPOSSUM::Web::State;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::CombinedResultSet;

use CGI::Carp qw(carpout);    # fatalsToBrowser;

use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_gene_ssa';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_gene_ssa";
$log_file .= "_devel" if DEVEL_VERSION;
$log_file .= "_$USER" if $USER;
$log_file .= ".log";

open(LOG, ">>$log_file") || die "Error opening log file $log_file - $!\n";

carpout(\*LOG);

sub setup
{
    my $self = shift;

    #print STDERR "setup\n";
    
    $self->start_mode('input');
    $self->mode_param('rm');
    $self->run_modes(
        'input'             => 'input',
        'results'           => 'results',
        'tfbs_details'      => 'tfbs_details',
        'text_results'      => 'text_results',
        'text_tfbs_details' => 'text_tfbs_details'
    );

    my $q = $self->query();
    
    my $sid = $q->param('sid');

    my $state;
    if ($sid) {
        #
        # Existing session. Load state from file.
        #
        my $filename = _session_tmp_file($sid);
        $state = OPOSSUM::Web::State->new(__Fn => $filename);
    } else {
        #
        # New session. Create new session ID and state object.
        #
        $sid = $$ . time;
        my $filename = _session_tmp_file($sid);

        $state = OPOSSUM::Web::State->new(
            -sid => $sid,
            __Fn => $filename
        );

        $self->initialize_state($state);
    }

    $self->state($state);

    #printf STDERR sprintf("\n\noPOSSUM State:\n%s\n\n",
    #    Data::Dumper::Dumper($self->state())
    #);

    unless ($self->opossum_db_connect()) {
        return $self->error("Could not connect to oPOSSUM DB");
    }
    
    #printf STDERR "\n\nrun mode = %s\n\n", $q->param('rm');
}


sub teardown
{
    my $self = shift;

    #print STDERR "teardown\n";

    if ($self->opdba()) {
        if ($self->opdba()->dbc()) {
            $self->opdba()->dbc()->disconnect();
        }
        $self->{-opdba} = undef;
    }

    my $state = $self->state();
    if ($state) {
        $state->dumper->Purity(1);
        $state->dumper->Deepcopy(1);
        $state->commit();
    }

    $self->_clean_tempfiles;
    $self->_clean_resultfiles;
}

sub input
{
    my $self = shift;

    #print STDERR "input\n";
    
    my $q = $self->query;
    #print STDERR "input query:\n"
    #    . Data::Dumper::Dumper($q);

    my $state = $self->state();

    my $species = $state->species() || $self->param('species');

    #
    # Retrieve any previously entered values from state.
    #
    my $in_t_gene_id_type = $state->t_gene_id_type();
    my @in_t_gene_ids = @{$state->t_gene_ids()} if $state->t_gene_ids();

    #
    # If gene IDs passed in via the query string from an external
    # application, these take precedent. These are always assumed
    # to be Ensembl IDs (gene ID type = 0).
    #
    foreach my $id ($q->param("id")) {
        push @in_t_gene_ids, $id;
        $in_t_gene_id_type = 0;
    }

    # this should be later modified to give the user a choice
    my $biotype = DFLT_BIOTYPE;
    #my $biotype = 'protein_coding' if $species eq 'worm';
    $biotype = _parse_biotype($biotype);
    $state->biotype($biotype);

    my $total_genes;
    if ($biotype) {
        $total_genes = $self->fetch_gene_count("biotype = '$biotype'");
    } else {
        $total_genes = $self->fetch_gene_count();
    }
    my $db_info = $self->fetch_db_info();
    my $xgid_types = $self->fetch_external_gene_id_types();
    my $cl_hash = $self->fetch_conservation_levels();
    my $thl_hash = $self->fetch_threshold_levels();
    my $srl_hash = $self->fetch_search_region_levels();

    # This doesn't work properly, db_info methods seem to be lost when 
    # retrieving from state later on DJA 19/10/2010
    #$state->db_info($db_info);

    my @cl_levels  = sort keys %$cl_hash;
    my @thl_levels = sort keys %$thl_hash;
    my @srl_levels = sort keys %$srl_hash;

    my $tax_group       = $db_info->tax_group();
    my $min_ic          = $db_info->min_ic();
    my $max_upstream_bp = $db_info->max_upstream_bp();

    # for worms, tax_group include nematodes,vertebrates,insects
    # you could give this option for insects as well
    # in the input menu, you need to be able to choose from 3 tax groups
    # easier to deal with if you separate them
    my @tax_groups = split /\s*,\s*/, $tax_group;
    my $num_tax_groups = scalar @tax_groups;

    $state->has_operon($db_info->has_operon());
    $state->tax_groups(\@tax_groups);
    $state->num_tax_groups($num_tax_groups);

    #
    # Connect to JASPAR DB and retrieve TF info
    #
    $self->jaspar_db_connect();

    my %core_tf_sets;
    my %pending_tf_sets;
    my %pbm_tf_sets;
    foreach my $tax_group (@tax_groups)
    {
        my $core_tf_set = $self->fetch_tf_set(
            -collection => 'CORE',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $core_tf_sets{$tax_group} = $core_tf_set;
        
        my $pending_tf_set = $self->fetch_tf_set(
            -collection => 'PENDING',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $pending_tf_sets{$tax_group} = $pending_tf_set;
        
        unless ($species eq 'yeast') {
            my $pbm_tf_set = $self->fetch_tf_set(
                -collection => 'PBM',
                -tax_group  => $tax_group,
                -min_ic     => $min_ic
            );
            $pbm_tf_sets{$tax_group} = $pbm_tf_set;
        }
    }

    #my $fam_tf_set;
    #unless ($species eq 'yeast' or $species eq 'worm') {
    #    $fam_tf_set = $self->fetch_tf_set(-collection => 'FAM');
    #}
    
    # so that the second db access won't be necessary 
    # why do I get tainted data error?
    #$state->core_tf_sets(\%core_tf_sets);
    #$state->pbm_tf_sets(\%pbm_tf_sets);
    #$state->pending_tf_sets(\%pending_tf_sets);
    #$state->fam_tf_set($fam_tf_set);

    #printf STDERR "\ncore_tf_set:\n"
    #    . Data::Dumper::Dumper($core_tf_set) . "\n\n";

    #printf STDERR "input_gene_ids:\n" . Data::Dumper::Dumper(\@input_gene_ids);

    my $vars = {
        abs_htdocs_path         => ABS_HTDOCS_PATH,
        rel_htdocs_path         => REL_HTDOCS_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        bg_color_class          => $state->bg_color_class(),
        title                   => $state->title(),
        heading                 => $state->heading(),
        section                 => 'Select Analysis Parameters',
        version                 => VERSION,
        devel_version           => DEVEL_VERSION,
        nresults                => NUM_RESULTS,
        dflt_nresults           => DFLT_NUM_RESULTS,
        zcutoffs                => ZSCORE_CUTOFFS,
        fcutoffs                => FISHER_CUTOFFS,
        dflt_zcutoff            => DFLT_ZSCORE_CUTOFF,
        dflt_fcutoff            => DFLT_FISHER_CUTOFF,
        cl_dflt_level           => DFLT_CONSERVATION_LEVEL,
        thl_dflt_level          => DFLT_THRESHOLD_LEVEL,
        srl_dflt_level          => DFLT_SEARCH_REGION_LEVEL,
        dflt_bg_num_rand_genes  => DFLT_BG_NUM_RAND_GENES,
        sid                     => $state->sid(),
        species                 => $species,
        db_info                 => $db_info,
        total_genes             => $total_genes,
        xgid_types              => $xgid_types,
        cl_hash                 => $cl_hash,
        thl_hash                => $thl_hash,
        srl_hash                => $srl_hash,
        cl_levels               => \@cl_levels,
        thl_levels              => \@thl_levels,
        srl_levels              => \@srl_levels,
        tax_groups              => \@tax_groups,
        num_tax_groups          => $num_tax_groups,
        core_tf_sets            => \%core_tf_sets,
        pbm_tf_sets             => \%pbm_tf_sets,
        pending_tf_sets         => \%pending_tf_sets,
        #fam_tf_set              => $fam_tf_set,
        in_t_gene_id_type       => $in_t_gene_id_type,
        in_t_gene_ids           => \@in_t_gene_ids,
        var_template            => "input_gene_ssa.html"
    };

    my $output = $self->process_template('master.html', $vars);
    #print STDERR "input results:\n"
    #    . Data::Dumper::Dumper($output);

    return $output;
}


sub results
{
    my $self = shift;

    #print STDERR "results\n";

    my $q = $self->query;
    #print STDERR "results query:\n"
    #    . Data::Dumper::Dumper($q);

    my $state = $self->state();

    my $analysis_type = 'default';

    my $species = $state->species() || $self->param('species');

    if (!defined $species) {
        return $self->error("Species not specified");
    }

    my $opdba = $self->opdba();

    #
    # This doesn't work properly, db_info methods seem to get lost
    # DJA 19/10/2010
    #
    #my $db_info = $state->db_info();

    my $dbia = $opdba->get_DBInfoAdaptor();
    my $db_info = $dbia->fetch_db_info();

    #my $tax_group   = $db_info->tax_group();
    my $min_ic      = $db_info->min_ic();

    #
    # Fetch all the input form entries and store them in the state
    #

    my $t_id_input_method = $q->param('t_id_input_method');
    if (!$t_id_input_method) {
        return $self->error("Target dene ID input method not specified");
    }
    $state->t_id_input_method($t_id_input_method);

    my $t_gene_id_type = $q->param("t_gene_id_type");
    if (!defined $t_gene_id_type) {
        return $self->error("Target gene ID type not specified");
    }
    $state->t_gene_id_type($t_gene_id_type);

    my $bg_id_input_method = $q->param('bg_id_input_method');
    if (!$bg_id_input_method) {
        return $self->error("Target gene ID input method not specified");
    }
    $state->bg_id_input_method($bg_id_input_method);

    my $bg_gene_id_type = DFLT_GENE_ID_TYPE;
    unless ($bg_id_input_method eq 'all') {
        if (!defined $bg_gene_id_type) {
            return $self->error("Background gene ID type not specified");
        }
    }
    $state->bg_gene_id_type($bg_gene_id_type);

    my $t_gene_text;
    if ($t_id_input_method eq 'paste') {
        $t_gene_text = $q->param('t_gene_text');

        return $self->error("No target gene IDs pasted") unless $t_gene_text;
    } elsif ($t_id_input_method eq 'upload') {
        my $gene_file = $q->param('t_gene_file');
        my $gene_fh   = $q->upload('t_gene_file');

        while (my $line = <$gene_fh>) {
            $t_gene_text .= $line;
        }

        return $self->error("No target gene IDs uploaded") unless $t_gene_text;
    }

    my $t_gene_ids = $self->parse_gene_id_text($t_gene_text);

    return $self->error("Error parsing target gene IDs") unless $t_gene_ids;

    $state->t_gene_ids($t_gene_ids);
    
    # t_operon_gids would be empty for species without any operons
    my (
        $t_gids,
        $t_included_gene_ids,
        $t_missing_gene_ids,
        #$t_gid_ensids,
        $t_gid_gene_ids,
        #$t_gene_id_gids,
        $t_operon_gids,
        $t_operon_unique_gids
    ) = $self->fetch_opossum_gene_ids(
            $t_gene_id_type, $t_gene_ids, $state->has_operon()
    );

    if (!$t_gids || !$t_gids->[0]) {
        return $self->error(sprintf("%s",
              "There was a problem fetching gene data for"
            . " the given gene target IDs.<br><br>"
            . " Please make sure that the gene ID type"
            . " selected matches the actual type of gene ID"
            . " entered.<br><br>"
            . " It is also possible that none of the genes entered"
            . " are stored in the oPOSSUM database. This could"
            . " be due to one or more of the following reasons:<br>"
            . " 1) The gene IDs entered do not map to Ensembl IDs<br>"
            . " 2) No significant conservation was obtained for these"
            . " genes<br><br>"
            . "Please see the FAQ for more information."
        ));
    }

    printf STDERR "t_gids:\n%s\n\n", Data::Dumper::Dumper($t_gids);
    printf STDERR "t_operon_gids:\n%s\n\n",
        Data::Dumper::Dumper($t_operon_gids);
    printf STDERR "t_operon_unique_gids:\n%s\n\n",
        Data::Dumper::Dumper($t_operon_unique_gids);

    $state->t_gids($t_gids);
    $state->t_included_gene_ids($t_included_gene_ids);
    $state->t_missing_gene_ids($t_missing_gene_ids);
    #$state->t_gid_ensids($t_gid_ensids);
    $state->t_gid_gene_ids($t_gid_gene_ids);
    #$state->t_gene_id_gids($t_gene_id_gids);
    $state->t_operon_gids($t_operon_gids);
    $state->t_operon_unique_gids($t_operon_unique_gids);

    my $bg_gids;
    my $bg_gene_ids;
    my $bg_included_gene_ids;
    my $bg_missing_gene_ids;
    my $bg_gid_gene_ids;
    my $bg_operon_gids;
    my $bg_operon_unique_gids;
    if ($bg_id_input_method eq 'paste') {
        my $bg_gene_text = $q->param('bg_gene_text');

        return $self->error("No background gene IDs pasted")
            unless $bg_gene_text;

        $bg_gene_ids = $self->parse_gene_id_text($bg_gene_text);

        return $self->error("Error parsing pasted background gene IDs")
            unless $bg_gene_ids;
    } elsif ($bg_id_input_method eq 'upload') {
        my $gene_file = $q->param('bg_gene_file');

        my $bg_gene_text;
        while (my $line = <$gene_file>) {
            $bg_gene_text .= $line;
        }

        return $self->error("No background gene IDs uploaded")
            unless $bg_gene_text;

        $bg_gene_ids = $self->parse_gene_id_text($bg_gene_text);

        return $self->error("Error parsing uploaded background gene IDs")
            unless $bg_gene_ids;
    } elsif ($bg_id_input_method eq 'random') {
        my $bg_num_rand_genes = $q->param('bg_num_rand_genes');

        return $self->error(
            "Number of random background gene IDs not specified"
        ) unless $bg_num_rand_genes;

        ($bg_gids,
         $bg_operon_gids,
         $bg_operon_unique_gids
        ) = $self->fetch_random_opossum_gene_ids(
            $bg_num_rand_genes, $state->has_operon(), $state->biotype
        );

        return $self->error(
            "Error fetching $bg_num_rand_genes random background gene IDs"
        ) unless $bg_gids;
    }

    $state->bg_gene_ids($bg_gene_ids);

    if ($bg_gene_ids) {
#        printf STDERR "\nbg gene IDs:\n%s\n\n",
#            Data::Dumper::Dumper($bg_gene_ids) if $bg_gene_ids;

        ($bg_gids,
         $bg_included_gene_ids,
         $bg_missing_gene_ids,
         $bg_gid_gene_ids,
         $bg_operon_gids,
         $bg_operon_unique_gids
        ) = $self->fetch_opossum_gene_ids(
            $bg_gene_id_type, $bg_gene_ids, $state->has_operon()
        );

        if (!$bg_gids || !$bg_gids->[0]) {
            return $self->error(sprintf("%s",
                  "There was a problem fetching gene data for"
                . " the given background gene IDs.<br><br>"
                . " Please make sure that the gene ID type"
                . " selected matches the actual type of gene ID"
                . " entered.<br><br>"
                . " It is also possible that none of the genes entered"
                . " are stored in the oPOSSUM database. This could"
                . " be due to one or more of the following reasons:<br>"
                . " 1) The gene IDs entered do not map to Ensembl IDs<br>"
                . " 2) No significant conservation was obtained for these"
                . " genes<br><br>"
                . "Please see the FAQ for more information."
            ));
        }
    }
    
    #printf STDERR "\nbg oPOSSUM gene IDs:\n%s\n\n",
    #    Data::Dumper::Dumper($bg_gids) if $bg_gids;

    $state->bg_gids($bg_gids);
    $state->bg_included_gene_ids($bg_included_gene_ids);
    $state->bg_missing_gene_ids($bg_missing_gene_ids);
    $state->bg_operon_gids($bg_operon_gids);
    $state->bg_operon_unique_gids($bg_operon_unique_gids);
    
    my $tf_select_method = $q->param("tf_select_method");
    if (!$tf_select_method) {
        return $self->error("TFBS profiles not specified");
    }
    $state->tf_select_method($tf_select_method);
    
    #
    # Divine the collection name from the tf_select_method parameter. The tf_select_method
    # parameter should be something like, e.g.:
    #   core_min_ic, fam_specific, pbm etc.
    #
    my $collection;
    if ($tf_select_method =~ /(\w+?)_/) {
        $collection = uc $1;
    } else {
        $collection = uc $tf_select_method;
    }
    $state->collection($collection);

    my $tf_select_criteria;
    if ($tf_select_method =~ /min_ic/) {
        $tf_select_criteria = 'min_ic';
    } elsif ($tf_select_method =~ /specific/) {
        $tf_select_criteria = 'specific';
    } else {
        $tf_select_criteria = 'all';
    }
    $state->tf_select_criteria($tf_select_criteria);

    my $core_min_ic;
    my $pbm_min_ic;
    my $pending_min_ic;
    my @core_tax_groups;
    my @pbm_tax_groups;
    my @pending_tax_groups;
    my $core_tax_group_str;
    my @core_tf_ids;
    my @pbm_tf_ids;
    my @pending_tf_ids;
    if ($tf_select_criteria eq 'min_ic') {
        my $num_tax_groups = $state->num_tax_groups();
        if ($collection eq "CORE") {
            $core_min_ic = $q->param('core_min_ic');
            if (!$core_min_ic) {
                return $self->error("No JASPAR CORE minimum IC specified");
            }

            if ($num_tax_groups > 1) {
                @core_tax_groups = $q->param("core_tax_groups");
                if (scalar @core_tax_groups == 0) {
                    return $self->error("No JASPAR CORE tax group selected");
                }
                $core_tax_group_str = join(',', @core_tax_groups);
            } else {
                @core_tax_groups = @{$state->tax_groups()};
                $core_tax_group_str = $db_info->tax_group();
            }
        } elsif ($collection eq "PBM") {
            $pbm_min_ic = $q->param('pbm_min_ic');
            if (!$pbm_min_ic) {
                return $self->error("No JASPAR PBM minimum IC specified");
            }

            #@pbm_tax_groups = $q->param("pbm_tax_groups");
            #if (scalar @pbm_tax_groups == 0) {
            #    return $self->error("No JASPAR PBM tax group selected");
            #}
        } elsif ($collection eq "PENDING") {
            $pending_min_ic = $q->param('pending_min_ic');
            if (!$pending_min_ic) {
                return $self->error("No JASPAR PENDING minimum IC specified");
            }

            #@pending_tax_groups = $q->param("pending_tax_groups");
            #if (scalar @pending_tax_groups == 0) {
            #    return $self->error("No JASPAR PENDING tax group selected");
            #}
        }
    } elsif ($tf_select_criteria eq 'specific') {
        if ($collection eq 'CORE') {
            push @core_tf_ids, $q->param('core_tfs');
            if (scalar(@core_tf_ids) == 0) {
                return $self->error("No specific $collection TFs selected");
            }
        } elsif ($collection eq 'PBM') {
            push @pbm_tf_ids, $q->param('pbm_tfs');
            if (scalar(@pbm_tf_ids) == 0) {
                return $self->error("No specific $collection TFs selected");
            }
        } elsif ($collection eq 'PENDING') {
            push @pending_tf_ids, $q->param('pending_tfs');
            if (scalar(@pending_tf_ids) == 0) {
                return $self->error("No specific $collection TFs selected");
            }
        }
    }

    $state->core_tf_ids(\@core_tf_ids);
    $state->pbm_tf_ids(\@pbm_tf_ids);
    $state->pending_tf_ids(\@pending_tf_ids);

    #print STDERR "results: TF IDS: core = " . scalar(@core_tf_ids) . "\tpbm = "
    #. scalar(@pbm_tf_ids) . "\tpending = " . scalar(@pending_tf_ids) . "\n";
    
    $state->core_min_ic($core_min_ic);
    $state->pbm_min_ic($pbm_min_ic);
    $state->pending_min_ic($pending_min_ic);

    $state->core_tax_groups(\@core_tax_groups);

    my $conservation_level;
    if ($species eq 'yeast') {
        $conservation_level = 1;
    } else {
        $conservation_level = $q->param('conservation_level');
    }

    my $cl_hash = $self->fetch_conservation_levels();
    my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();

    $state->conservation_level($conservation_level);
    $state->min_conservation($min_conservation);

    my $threshold_level = $q->param('threshold_level');
    my $threshold       = $q->param('threshold') / 100;

    my $thl_hash = $self->fetch_threshold_levels();

    if ($threshold) {
        #
        # Custom threshold entered
        #
        $analysis_type = 'custom';

        $threshold_level = undef;
        $state->threshold_level($threshold_level);
        $state->threshold($threshold);

        #
        # Check if entered threshold corresponds to a default threshold level
        #
        foreach my $level (keys %$thl_hash) {
            if ($threshold == $thl_hash->{$level}->threshold()) {
                $state->threshold_level($level);
                $analysis_type = 'default';
                last;
            }
        }
    } else {
        $state->threshold_level($threshold_level);
        $state->threshold($thl_hash->{$threshold_level}->threshold());
    }

    my $search_region_level = $q->param('search_region_level');
    my $upstream_bp         = $q->param('upstream_bp');
    my $downstream_bp       = $q->param('downstream_bp');

    my $srl_hash = $self->fetch_search_region_levels();
    
    if (   (defined $upstream_bp && $upstream_bp ne '' && $upstream_bp > 0)
        || (defined $downstream_bp && $downstream_bp ne ''
            && $downstream_bp > 0)
       )
    {
        #
        # Custom search region entered
        #
        $analysis_type = 'custom';

        $upstream_bp = 0 if !defined $upstream_bp;
        $downstream_bp = 0 if !defined $downstream_bp && $species ne 'yeast';

        $search_region_level = undef;
        $state->search_region_level($search_region_level);
        $state->upstream_bp($upstream_bp);
        $state->downstream_bp($downstream_bp);

        #print STDERR "Upstream: $upstream_bp\tDownstream: $downstream_bp\n";
        
        #
        # Check if entered search region corresponds to a default search
        # region level
        #
        foreach my $level (keys %$srl_hash) {
            if ($species eq 'yeast') {
                if ($upstream_bp == $srl_hash->{$level}->upstream_bp()) {
                    $state->search_region_level($level);
                    $analysis_type = 'default';
                    last;
                }
            } else {
                if (   $upstream_bp == $srl_hash->{$level}->upstream_bp()
                    && $downstream_bp == $srl_hash->{$level}->downstream_bp()
                )
                {
                    $state->search_region_level($level);
                    $analysis_type = 'default';
                    last;
                }
            }
        }
    } else {
        #
        # default search regions
        #
        $state->search_region_level($search_region_level);

        $state->upstream_bp($srl_hash->{$search_region_level}->upstream_bp());
        $state->downstream_bp(
            $srl_hash->{$search_region_level}->downstream_bp()
        );
        
    }
    $state->analysis_type($analysis_type);

    
    my $result_type = $q->param('result_type');
    $state->result_type($result_type);
    if ($result_type eq 'top_x_results') {
        $state->num_display_results($q->param('num_display_results'));
        $state->zscore_cutoff(undef);
        $state->fisher_cutoff(undef);
    } elsif ($result_type eq 'significant_hits') {
        $state->num_display_results(undef);
        $state->zscore_cutoff($q->param('zscore_cutoff'));
        $state->fisher_cutoff($q->param('fisher_cutoff'));
    }

    $state->result_sort_by($q->param('result_sort_by'));


    #
    # Retrieve information from oPOSSUM and JASPAR based on user input values
    #

    #print STDERR "before tf check:\n" . Data::Dumper::Dumper($self);

    #
    # Connect to JASPAR DB and retrieve TF info
    #
    
    # Here, you are retrieving tfs for the second time.
    # this could be reworked to minimize db access

    $self->jaspar_db_connect();

    my %matrix_args = (
        -collection => $collection,
    );

    my $whole_tf_sets;
    if ($collection eq 'CORE') {
        $whole_tf_sets = $state->core_tf_sets;
        if (@core_tf_ids) {
            $matrix_args{-ID} = \@core_tf_ids;
        } else {
            $matrix_args{-tax_group} = \@core_tax_groups;
            $matrix_args{-min_ic} = $core_min_ic || $min_ic;
        }
    } elsif ($collection eq 'PBM') {
        $whole_tf_sets = $state->pbm_tf_sets;
        if (@pbm_tf_ids) {
            $matrix_args{-ID} = \@pbm_tf_ids;
        } else {
            $matrix_args{-min_ic} = $pbm_min_ic || $min_ic;
        }
    } elsif ($collection eq 'PENDING') {
        $whole_tf_sets = $state->pending_tf_sets;
        if (@pending_tf_ids) {
            $matrix_args{-ID} = \@pending_tf_ids;
        } else {
            $matrix_args{-min_ic} = $pending_min_ic || $min_ic;
        }
    }

    my $tf_set = $self->fetch_tf_set(%matrix_args);
    # let's try minimizing db access
    # but I keep getting tainted data error
    #my $tf_set = OPOSSUM::TFSet->new();
    #foreach my $tax_group (keys %$whole_tf_sets)
    #{
    #    my $whole_set = $$whole_tf_sets{$tax_group};
    #    my $sub_tf_set = $whole_set->subset(
    #        -ids => $matrix_args{-ID},
    #        -collections => $collection,
    #        -tax_groups => $tax_group,
    #        -min_ic => $matrix_args{-min_ic}
    #    );
    #    $tf_set->add_matrix_list($sub_tf_set->get_matrix_list);
    #}
    
    if (!$tf_set) {
        return $self->error("Error fetching JASPAR TFBS profiles");
    }

    #
    # Don't store. Cause of taint errors with persitent state object?
    #
    #$state->tf_set($tf_set);

    my $aca = $opdba->get_AnalysisCountsAdaptor();
    if (!$aca) {
        return $self->error("Could not get AnalysisCountsAdaptor");
    }

    #printf STDERR sprintf("\n\noPOSSUM State:\n%s\n\n",
    #    Data::Dumper::Dumper($self->state())
    #);
    
    #
    # fetch target and background conserved region GC content
    #
    my $t_cr_gc_content = $self->fetch_cr_gc_content(
        $state->t_gids, $state->conservation_level,
        $state->upstream_bp, $state->downstream_bp,
        $state->biotype);
    
    my $bg_cr_gc_content = $self->fetch_cr_gc_content(
        $state->bg_gids, $state->conservation_level,
        $state->upstream_bp, $state->downstream_bp,
        $state->biotype);
    
    # if operon genes present, the retrieved gene counts are all based on
    # the first gene search region, taken care by the CountsAdaptor.
    # no further action necessary.
    my $t_counts;
    my $bg_counts;
    my $t_cr_length;
    my $bg_cr_length;
    my $crla = $opdba->get_ConservedRegionLengthAdaptor();
    if ($analysis_type eq 'default') {
        $t_counts = $self->fetch_analysis_counts(
            $analysis_type,
            -gene_ids               => $state->t_gids(),
            -has_operon             => $state->has_operon(), #optional
            -operon_gene_ids        => $state->t_operon_gids(),
            -tf_ids                 => $tf_set->tf_ids(),
            -conservation_level     => $state->conservation_level(),
            -threshold_level        => $state->threshold_level(),
            -search_region_level    => $state->search_region_level()
        );

        return $self->error("Error fetching target gene TFBS counts")
            unless $t_counts;

        $bg_counts = $self->fetch_analysis_counts(
            $analysis_type,
            -gene_ids               => $state->bg_gids(),
            -has_operon             => $state->has_operon(), #not optional, if choosing all genes
            -operon_gene_ids        => $state->bg_operon_gids(),
            -tf_ids                 => $tf_set->tf_ids(),
            -conservation_level     => $state->conservation_level(),
            -threshold_level        => $state->threshold_level(),
            -search_region_level    => $state->search_region_level()
        );

        return $self->error("Error fetching background gene TFBS counts")
            unless $bg_counts;

        $t_cr_length = $crla->fetch_total_length(
            -conservation_level     => $state->conservation_level(),
            -search_region_level    => $state->search_region_level(),
            -gene_ids               => $state->t_gids(),
            -operon_gene_ids        => $state->t_operon_gids(),
            -has_operon             => $state->has_operon()
        );

        return $self->error(
            "Error fetching target gene total conserved region length"
        ) unless $t_cr_length;

        $bg_cr_length = $crla->fetch_total_length(
            -conservation_level     => $state->conservation_level(),
            -search_region_level    => $state->search_region_level(),
            -gene_ids               => $state->bg_gids(),
            -operon_gene_ids        => $state->bg_operon_gids(),
            -has_operon             => $state->has_operon()
        );
        
        return $self->error(
            "Error fetching background gene total conserved region length"
        ) unless $bg_cr_length;

    } elsif ($analysis_type eq 'custom') {

        $t_counts = $self->fetch_analysis_counts(
            $analysis_type,
            -gene_ids           => $state->t_gids(),
            -has_operon         => $state->has_operon(), # optional
            -operon_gene_ids    => $state->t_operon_gids(),
            -tf_ids             => $tf_set->tf_ids(),
            -conservation_level => $state->conservation_level(),
            -threshold          => $state->threshold(),
            -upstream_bp        => $state->upstream_bp(),
            -downstream_bp      => $state->downstream_bp()
        );

        return $self->error("Error fetching target gene TFBS counts")
            unless $t_counts;

        $bg_counts = $self->fetch_analysis_counts(
            $analysis_type,
            -gene_ids           => $state->bg_gids(),
            -has_operon         => $state->has_operon(), # not optional
            -operon_gene_ids    => $state->bg_operon_gids(),
            -tf_ids             => $tf_set->tf_ids(),
            -conservation_level => $state->conservation_level(),
            -threshold          => $state->threshold(),
            -upstream_bp        => $state->upstream_bp(),
            -downstream_bp      => $state->downstream_bp()
        );

        return $self->error("Error fetching background gene TFBS counts")
            unless $bg_counts;

        $t_cr_length = $crla->fetch_total_length(
            -conservation_level     => $state->conservation_level(),
            -upstream_bp            => $state->upstream_bp(),
            -downstream_bp          => $state->downstream_bp(),
            -gene_ids               => $state->t_gids(),
            -operon_gene_ids        => $state->t_operon_gids(),
            -has_operon             => $state->has_operon()
        );

        return $self->error(
            "Error fetching target gene total conserved region length"
        ) unless $t_cr_length;

        $bg_cr_length = $crla->fetch_total_length(
            -conservation_level     => $state->conservation_level(),
            -upstream_bp            => $state->upstream_bp(),
            -downstream_bp          => $state->downstream_bp(),
            -gene_ids               => $state->bg_gids(),
            -operon_gene_ids        => $state->bg_operon_gids(),
            -has_operon             => $state->has_operon()
        );

        return $self->error(
            "Error fetching background gene total conserved region length"
        ) unless $bg_cr_length;
    }

    #printf STDERR "\n\nt_counts:\n" . Data::Dumper::Dumper($t_counts) . "\n\n";

    #printf STDERR "\n\nbg_counts:\n" . Data::Dumper::Dumper($bg_counts)
    #    . "\n\n";

    #$state->t_counts($t_counts);
    #$state->bg_counts($bg_counts);

    my $fisher = OPOSSUM::Analysis::Fisher->new();
    return $self->error("Error initializing Fisher analysis") unless $fisher;

    my $fresult_set = $fisher->calculate_Fisher_probability(
        $bg_counts,
        $t_counts
    );
    return $self->error("Error performing Fisher analysis")
        unless $fresult_set;

    my $zscore = OPOSSUM::Analysis::Zscore->new();
    return $self->error("Error initializing z-score analysis") unless $zscore;

    my $zresult_set = $zscore->calculate_Zscore(
        $bg_counts,
        $t_counts,
        $bg_cr_length,
        $t_cr_length,
        $tf_set
    );
    return $self->error("Error computing z-score") unless $zresult_set;

    my $result_set = OPOSSUM::Analysis::CombinedResultSet->new(
        -fisher_result_set  => $fresult_set,
        -zscore_result_set  => $zresult_set
    );

    return $self->error("Error combining Fisher and z-score result_set")
        unless $result_set;

    #printf STDERR "\nresult set:\n%s\n\n", Data::Dumper::Dumper($result_set);

    $state->result_set($result_set);

    #
    # Fisher scores are now -ln() values so all results are reverse sorted
    # (largest value first).
    # DJA 14/03/2011
    #
    #my $result_sort_reverse = 0;
    #if ($state->result_sort_by() eq 'zscore') {
    #    $result_sort_reverse = 1;
    #}
    my $result_sort_reverse = 1;

    my $results = $result_set->get_list(
        -num_results    => $state->num_display_results(),
        -zscore_cutoff  => $state->zscore_cutoff(),
        -fisher_cutoff  => $state->fisher_cutoff(),
        -sort_by        => $state->result_sort_by(),
        -reverse        => $result_sort_reverse
    );

    #printf STDERR "\nresult list:\n%s\n\n", Data::Dumper::Dumper($results);

    #
    # Create a temporary working directory for all the temp. input files and
    # output results file as a sub-directory under the defined temp. dir.
    #
    my $results_dir = tempdir(DIR => ABS_HTDOCS_RESULTS_PATH);

    my $results_subdir = $results_dir;
    $results_subdir =~ s/.*\///;

    $state->results_subdir($results_subdir);

    my $abs_result_file = "$results_dir/" . RESULTS_TEXT_FILENAME;

    if ($results && $results->[0]) {
        $self->write_results($abs_result_file, $results, $tf_set);
    }
    
    #print STDERR "results self:\n" . Data::Dumper::Dumper($self);

    # passed through to tfbs_details routine
    #my $gene_id_url = '&gene_id=';
    #$gene_id_url .= join '&gene_id=', @{$self->target_gene_ids};

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        bg_color_class      => $state->bg_color_class(),
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        rel_htdocs_results_path => REL_HTDOCS_RESULTS_PATH,
        result_file         => RESULTS_TEXT_FILENAME,
        jaspar_url          => JASPAR_URL,
        low_matrix_ic       => LOW_MATRIX_IC,
        high_matrix_ic      => HIGH_MATRIX_IC,
        low_matrix_gc       => LOW_MATRIX_GC,
        high_matrix_gc      => HIGH_MATRIX_GC,
        low_seq_gc          => LOW_SEQ_GC,
        high_seq_gc         => HIGH_SEQ_GC,
        title               => $state->title(),
        heading             => $state->heading(),
        section             => 'Analysis Results',
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        result_retain_days  => REMOVE_TEMPFILES_OLDER_THAN,
        analysis_type       => $analysis_type,
        sid                 => $state->sid(),
        species             => $species,
        gene_id_type        => $state->t_gene_id_type(),
        tf_select_criteria  => $state->tf_select_criteria(),
        collection          => $state->collection(),
        core_min_ic         => $state->core_min_ic(),
        pbm_min_ic          => $state->pbm_min_ic(),
        pending_min_ic      => $state->pending_min_ic(),
        core_tax_group_str  => $core_tax_group_str,
        conservation_level  => $state->conservation_level(),
        threshold_level     => $state->threshold_level(),
        search_region_level => $state->search_region_level(),
        min_conservation    => $state->min_conservation(),
        threshold           => $state->threshold(),
        upstream_bp         => $state->upstream_bp(),
        downstream_bp       => $state->downstream_bp(),
        result_type         => $state->result_type(),
        num_display_results => $state->num_display_results(),
        result_sort_by      => $state->result_sort_by(),
        zscore_cutoff       => $state->zscore_cutoff(),
        fisher_cutoff       => $state->fisher_cutoff(),
        t_gids              => $state->t_gids(),
        t_gene_ids          => $state->t_gene_ids(),
        t_operon_gids       => $state->t_operon_gids(),
        t_included_gene_ids => $state->t_included_gene_ids(),
        t_missing_gene_ids  => $state->t_missing_gene_ids(),
        bg_gids             => $state->bg_gids(),
        bg_gene_ids         => $state->bg_gene_ids(),
        bg_included_gene_ids=> $state->bg_included_gene_ids(),
        bg_missing_gene_ids => $state->bg_missing_gene_ids(),
        num_t_gids          => scalar @{$state->t_gids()},
        num_t_gene_ids      => scalar @{$state->t_gene_ids()},
        num_t_included_gene_ids => scalar @{$state->t_included_gene_ids()},
        num_t_missing_gene_ids  => $state->t_missing_gene_ids()
                                    ? scalar @{$state->t_missing_gene_ids()}
                                    : 0,
        tf_set              => $tf_set,
        t_cr_gc_content     => $t_cr_gc_content,
        bg_cr_gc_content    => $bg_cr_gc_content,
        target_counts       => $t_counts,
        results_subdir      => $results_subdir,
        results             => $results,
        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'N/A'
                               },
        formatg             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*g", $dec, $f)
                                        : 'N/A'
                               },
        var_template        => "results_gene_ssa.html"
    };

    #print STDERR "results vars:\n" . Data::Dumper::Dumper($vars);

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

sub tfbs_details
{
    my $self = shift;

    #print STDERR "tfbs_details\n";
    #print STDERR "tfbs_details self:\n" . Data::Dumper::Dumper($self);

    my $q = $self->query;

    my $tf_id = $q->param('tf_id');

    #print STDERR "tfbs_details query\n" . Data::Dumper::Dumper($q);

    my $state = $self->state();

    my $species             = $state->species();
    my $collection          = $state->collection();
    my $conservation_level  = $state->conservation_level();
    my $min_conservation    = $state->min_conservation();
    my $threshold_level     = $state->threshold_level();
    my $threshold           = $state->threshold();
    my $search_region_level = $state->search_region_level();
    my $upstream_bp         = $state->upstream_bp();
    my $downstream_bp       = $state->downstream_bp();
    my $gene_id_type        = $state->t_gene_id_type();
    my $t_gids              = $state->t_gids();
    my $t_gene_ids          = $state->t_gene_ids();
    my $gid_gene_ids        = $state->t_gid_gene_ids();
    #my $gene_id_gids        = $state->t_gene_id_gids(); # is this even used?
    #my $gid_ensids          = $state->t_gid_ensids();
    my $has_operon          = $state->has_operon();
    my $t_operon_gids       = $state->t_operon_gids();
    my $t_operon_unique_gids= $state->t_operon_unique_gids();
    my $results_subdir      = $state->results_subdir();

    #
    # Cause of persitent state taint errors?
    #
    #my $tf_set              = $state->tf_set();
    #my $tf = $tf_set->get_tf($tf_id);

    my $jdbh = $self->jdbh();
    unless ($jdbh) {
        $jdbh = $self->jaspar_db_connect();
    }
    my $tf = $jdbh->get_Matrix_by_ID($tf_id, 'PFM');

    my $opdba = $self->opdba();
    unless ($opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $ctfsa = $opdba->get_ConservedTFBSAdaptor();
    my $ga = $opdba->get_GeneAdaptor();

    #
    # Fetch all TFBSs for this TF for the set of input genes
    #
    #my $tfbss = $ctfsa->fetch_all_by_tf_id(
    #    -tf_id              => $tf_id,
    #    -gene_ids           => $t_gids,
    #    -conservation_level => $conservation_level,
    #    -threshold          => $threshold,
    #    -upstream_bp        => $upstream_bp,
    #    -downstream_bp      => $downstream_bp
    #);

    #printf STDERR "\ntfbs genes:\n" . Data::Dumper::Dumper($tfbss) . "\n\n";
    
    #
    # Get unique list of genes which have TFBSs for this TF (not all input
    # genes may have a TFBS for this TF). Also map TFBSs to oPOSSUM gene IDs.
    #
    # for operon genes, key = gid, val = first gene gid
    # t_gids contain trailing operon genes as well as first genes
    # tf_genes will be passed to write_tfbs_details
    # for operon genes, tf_genes should contain only the first genes
    # operon_genes, keyed by the first gene id, contains all operon genes
    # that are in t_gids. If the first gene is not in t_gids, it will not be
    # included. (but key = still first gene id)
    # 
    #my @tf_genes;
    #my %gid_tfbss;
    #my %operon_genes; # key = first gene id, val = operon genes
    #foreach my $tfbs (@$tfbss) {
    #    my $gid = $tfbs->gene_id();

    #    unless ($gid_tfbss{$gid}) {
    #        my $gene = $ga->fetch_by_id($gid);
    #        if ($$t_operon_gids{$gid}) {
    #            my $fgid = $$t_operon_gids{$gid};
    #            if ($gid == $fgid) { # this gene is the first operon gene
    #                push @tf_genes, $gene;
    #                push @{$gid_tfbss{$gid}}, $tfbs;
    #            } elsif (!$operon_genes{$fgid}) {
    #                # new operon
    #                my $fgene = $ga->fetch_by_id($fgid);
    #                push @tf_genes, $fgene;
    #            }
    #            push @{$operon_genes{$fgid}}, $gene;
    #            
    #        } else {
    #            push @tf_genes, $gene;
    #        }
    #    }

    #    # if operon gene, already taken care of in the above unless {}
    #    unless ($$t_operon_gids{$gid}) {
    #        push @{$gid_tfbss{$gid}}, $tfbs;
    #    }
    #}

    my @tf_genes;
    my %gid_tfbss;
    my %operon_genes; # key = first gene id, val = operon genes
    if ($has_operon) {
        foreach my $gid (@$t_operon_unique_gids) {
            my $tfbss = $ctfsa->fetch_by_gene(
                -tf_id              => $tf_id,
                -gene_id            => $gid,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );

            if ($tfbss) {
                $gid_tfbss{$gid} = $tfbss;

                my $gene = $ga->fetch_by_id($gid);
                push @tf_genes, $gene;
            }
        }

        foreach my $gid (@$t_gids) {
            my $fgid = $t_operon_gids->{$gid};

            if ($fgid && $gid_tfbss{$fgid}) {
                my $gene = $ga->fetch_by_id($gid);
                push @{$operon_genes{$fgid}}, $gene;
            }
        }
    } else {
        #
        # Species does not have operons.
        #
        foreach my $gid (@$t_gids) {
            my $tfbss = $ctfsa->fetch_by_gene(
                -tf_id              => $tf_id,
                -gene_id            => $gid,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );

            if ($tfbss) {
                $gid_tfbss{$gid} = $tfbss;

                my $gene = $ga->fetch_by_id($gid);
                push @tf_genes, $gene;
            }
        }
    }

    #printf STDERR "\ntfbs_details tf_genes:\n%s\n\n"
    #   . Data::Dumper::Dumper(@tf_genes);

    #printf STDERR "\ntfbs_details operon_genes:\n%s\n\n"
    #   . Data::Dumper::Dumper(%operon_genes);

    my $result_file = $tf->ID() . ".txt";

    my $abs_result_file = sprintf("%s/%s/%s",
        ABS_HTDOCS_RESULTS_PATH,
        $results_subdir,
        $result_file
    );
    
    $self->write_tfbs_details(
        $abs_result_file,
        $tf,
        \@tf_genes,
        \%gid_tfbss,
        $gene_id_type,
        $gid_gene_ids,
        $has_operon,
        \%operon_genes
    );

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        bg_color_class      => $state->bg_color_class(),
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        rel_htdocs_results_path => REL_HTDOCS_RESULTS_PATH,
        jaspar_url          => JASPAR_URL,
        low_matrix_ic       => LOW_MATRIX_IC,
        high_matrix_ic      => HIGH_MATRIX_IC,
        low_matrix_gc       => LOW_MATRIX_GC,
        high_matrix_gc      => HIGH_MATRIX_GC,
        low_seq_gc          => LOW_SEQ_GC,
        high_seq_gc         => HIGH_SEQ_GC,
        title               => $state->title(),
        heading             => $state->heading(),
        section             => $tf->name() . ' Binding Site Details',
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        dflt_gene_id_type   => DFLT_GENE_ID_TYPE,
        result_retain_days  => REMOVE_TEMPFILES_OLDER_THAN,
        tf                  => $tf,
        collection          => $collection,
        genes               => \@tf_genes,
        gid_tfbss           => \%gid_tfbss,
        gene_id_type        => $gene_id_type,
        gid_gene_ids        => $gid_gene_ids,
        #gene_id_gids        => $gene_id_gids,
        has_operon          => $state->has_operon(),
        operon_genes        => \%operon_genes,
        results_subdir      => $results_subdir,
        result_file         => $result_file,
        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'N/A'
                               }
    };

    my $output = $self->process_template('tfbs_details_gene_ssa.html', $vars);
    #print STDERR "tfbs_details output\n" . Data::Dumper::Dumper($output);
    return $output;
}

sub text_results
{
    my $self = shift;

    my $q = $self->query;

    my $state = $self->state();

    my $species = $state->species();

    my $result_text_file = $q->param('result_text_file');

    my $vars = {
        abs_htdocs_path  => ABS_HTDOCS_PATH,
        rel_htdocs_path  => REL_HTDOCS_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        bg_color_class   => $state->bg_color_class(),
        title            => $state->title(),
        heading          => $state->heading(),
        section          => 'Analysis Results',
        version          => VERSION,
        devel_version    => DEVEL_VERSION,
        text_file        => $result_text_file
    };

    my $output = $self->process_template('blank.html', $vars);

    return $output;
}

sub text_tfbs_details
{
    my $self = shift;

    my $q = $self->query;

    my $state = $self->state();

    my $species = $state->species();

    my $tf_gene_text = $q->param('tf_gene_text');

    my $vars = {
        abs_htdocs_path  => ABS_HTDOCS_PATH,
        rel_htdocs_path  => REL_HTDOCS_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        bg_color_class   => $state->bg_color_class(),
        title            => $state->title(),
        heading          => $state->heading(),
        section          => 'TFBS Details',
        version          => VERSION,
        devel_version    => DEVEL_VERSION,
        text             => $tf_gene_text
    };

    my $output = $self->process_template('blank.txt', $vars);
    #my $template = Template->new();
    #my $input = "[%text%]";

    #my $output;
    #$template->process(\$input, $vars, \$output);
    #print STDERR "text_tfbs_details:\n$string\n";

    $self->header_props(-type => 'text/plain');
    return $output;
}

sub initialize_state
{
    my ($self, $state) = @_;

    my $species = $self->param('species');
    
    unless ($species) {
        return $self->error("Species is undefined");
    }
    $state->species($species);

    my $heading = sprintf "%s Single Site Analysis",
        ucfirst $state->species();
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    $state->errors(undef);
    $state->warnings(undef);
}

1;
