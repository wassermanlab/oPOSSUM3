package oPossumGeneTCAWeb;

use base 'CGI::Application';

use oPossumWebInclude;
use oPossumGeneWebInclude;
use oPossumTCAWebInclude;

use Data::Dumper;    # for debugging only
use File::Temp qw/ tempfile tempdir /;

use Bio::Seq;

use OPOSSUM::Web::State;
use OPOSSUM::Analysis::Cluster::Fisher;
use OPOSSUM::Analysis::Cluster::Zscore;
use OPOSSUM::Analysis::Cluster::CombinedResultSet;

use CGI::Carp qw(carpout);    # fatalsToBrowser;

use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_gene_tca';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_gene_tca";
$log_file .= '_devel' if DEVEL_VERSION;
$log_file .= '_$USER' if $USER;
$log_file .= ".log";

open(LOG, ">>$log_file") || die "Error opening open log file $log_file - $!\n";

carpout(\*LOG);

sub setup
{
    my $self = shift;

    #print STDERR "setup\n";
    
    $self->start_mode('input');
    $self->mode_param('rm');
    $self->run_modes(
        'input'                     => 'input',
        'results'                   => 'results',
        'tfbs_cluster_details'      => 'tfbs_cluster_details',
        'text_results'              => 'text_results',
        'text_tfbs_cluster_details' => 'text_tfbs_cluster_details',
		'tfbs_cluster_info'         => 'tfbs_cluster_info'
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

        $self->state($state);
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
    
    unless ($self->opossum_cluster_db_connect()) {
        return $self->error("Could not connect to oPOSSUM_cluster DB");
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
    
    if ($self->cldba()) {
        if ($self->cldba()->dbc()) {
            $self->cldba()->dbc()->disconnect();
        }
        $self->{-cldba} = undef;
    }

    my $state = $self->state();
    if ($state) {
        $state->dumper->Purity(1);
        $state->dumper->Deepcopy(1);
        $state->commit();
    }

    $self->_clean_tempfiles;
}


sub input
{
    my $self = shift;

    #print STDERR "input\n";
    
    my $q = $self->query;
    #print STDERR "input query:\n"
    #		    . Data::Dumper::Dumper($q);

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

    my @cl_levels  = sort keys %$cl_hash;
    my @thl_levels = sort keys %$thl_hash;
    my @srl_levels = sort keys %$srl_hash;

    my $max_upstream_bp = $db_info->max_upstream_bp();

	$state->has_operon($db_info->has_operon());

    #
    # Connect to oPOSSUM_cluster DB and retrieve TFCluster info
    #
    $self->opossum_cluster_db_connect();
	$self->jaspar_db_connect();

    # instead of showing individual clusters to the user, why not
    # have the user specify the class and family instead?
    my $tf_cluster_set = $self->fetch_tf_cluster_set();
    #print STDERR "tf cluster set size = " . $tf_cluster_set->size . "\n";
    #my $tf_families = $tf_cluster_set->get_tf_families;
    #print STDERR Data::Dumper::Dumper($tf_families) . "\n";
    
    # let's keep it around so that the second db access won't be necessary 
    #$state->tf_cluster_set($tf_cluster_set);

    #printf STDERR "\ntf_cluster_set:\n"
    #    . Data::Dumper::Dumper($tf_cluster_set) . "\n\n";

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
        tf_cluster_set          => $tf_cluster_set,
        in_t_gene_id_type       => $in_t_gene_id_type,
        in_t_gene_ids           => \@in_t_gene_ids,
        var_template            => "input_gene_tca.html"
    };

    my $output = $self->process_template('master.html', $vars);
    #print STDERR "input results:\n"
    #		    . Data::Dumper::Dumper($output);

    return $output;
}


sub results
{
    my $self = shift;

    #print STDERR "results\n";

    my $q = $self->query;
    #print STDERR "results query:\n"
    #		    . Data::Dumper::Dumper($q);

    my $state = $self->state();

    my $analysis_type = 'default';

    my $species = $state->species() || $self->param('species');

    if (!defined $species) {
        return $self->error("Species not specified");
    }

    # OPOSSUM
    my $opdba = $self->opdba();
    my $dbia = $opdba->get_DBInfoAdaptor();
    my $db_info = $dbia->fetch_db_info();
    
    # TFBSCluster
    my $cldba = $self->cldba();
    my $cldbia = $cldba->get_DBInfoAdaptor();
    my $cl_db_info = $cldbia->fetch_db_info();

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

        while (my $line = <$gene_file>) {
            $t_gene_text .= $line;
        }

        return $self->error("No target gene IDs uploaded") unless $t_gene_text;
    }

    my $t_gene_ids = $self->parse_gene_id_text($t_gene_text);

    return $self->error("Error parsing target gene IDs") unless $t_gene_ids;

    $state->t_gene_ids($t_gene_ids);

	# t_operon_gids would be empty for species without any operons
    my ($t_gids,
		$t_included_gene_ids,
		$t_missing_gene_ids,
        $t_gid_gene_ids,
		$t_operon_gids
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

    $state->t_gids($t_gids);
    $state->t_included_gene_ids($t_included_gene_ids);
    $state->t_missing_gene_ids($t_missing_gene_ids);
    #$state->t_gid_ensids($t_gid_ensids);
    $state->t_gid_gene_ids($t_gid_gene_ids);
    #$state->t_gene_id_gids($t_gene_id_gids);
    $state->t_operon_gids($t_operon_gids);

    my $bg_gids;
    my $bg_gene_ids;
    my $bg_included_gene_ids;
    my $bg_missing_gene_ids;
    my $bg_gid_gene_ids;
    my $bg_operon_gids;
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

        ($bg_gids, $bg_operon_gids) = $self->fetch_random_opossum_gene_ids(
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
		 $bg_operon_gids) = $self->fetch_opossum_gene_ids(
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
	
    my $cluster_type = $q->param("cluster_type");
    if (!$cluster_type) {
        return $self->error("TFBS clusters not specified");
    }
    $state->cluster_type($cluster_type);

    my $cl_select_criteria;
    if ($cluster_type =~ /specific/) {
        $cl_select_criteria = 'specific';
    } else {
        $cl_select_criteria = 'all';
    }
    $state->cl_select_criteria($cl_select_criteria);

    my @tf_cluster_families;
    if ($cl_select_criteria eq 'specific')
    {
        
        push @tf_cluster_families, $q->param('tf_cluster_families');
        if (!@tf_cluster_families) {
            return $self->error("No specific TFBS cluster families selected");
        }
    }
    $state->tf_cluster_families(\@tf_cluster_families);

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
        $state->search_region_level($search_region_level);

        $state->upstream_bp($srl_hash->{$search_region_level}->upstream_bp());
        $state->downstream_bp(
            $srl_hash->{$search_region_level}->downstream_bp()
        );
    }
    $state->analysis_type($analysis_type);

    my $result_type = $q->param('result_type');
    $state->result_type($result_type);
    if ($result_type eq 'top') {
        $state->num_display_results($q->param('num_display_results'));
        $state->zscore_cutoff(undef);
        $state->fisher_cutoff(undef);
    } elsif ($result_type eq 'significant') {
        $state->num_display_results(undef);
        $state->zscore_cutoff($q->param('zscore_cutoff'));
        $state->fisher_cutoff($q->param('fisher_cutoff'));
    }

    $state->result_sort_by($q->param('result_sort_by'));


    #
    # Retrieve information from oPOSSUM and TFBSCluster based on user input values
    #

    #print STDERR "before tf check:\n" . Data::Dumper::Dumper($self);

    #
    # Connect to TFBSCluster DB and retrieve TFCluster info
    #
    
    # Here, you are retrieving the clusters for the second time.
    # this could be reworked to minimize db access
    # so in input(), I put tf_cluster_sets into state
    # key = tax group, val = tf set for the tax group (OPOSSUM::TFSet)

    $self->opossum_cluster_db_connect();
    
    my %matrix_args;
    if (@tf_cluster_families) {
        $matrix_args{-families} = \@tf_cluster_families;
    }

    my $tf_cluster_set = $self->fetch_tf_cluster_set(%matrix_args);

    if (!$tf_cluster_set) {
        return $self->error("Error fetching TFBSCluster::TFClusterSet");
    }

	# deleted the code for keepting the tf_cluster_set in state to minimize
	# db access, as I keep getting tainted data error
	# if needed, look at previous backup files
	
    $state->tf_cluster_set($tf_cluster_set);

    my $opdba = $self->opdba();

    my $aca = $opdba->get_AnalysisClusterCountsAdaptor();
    if (!$aca) {
        return $self->error("Could not get AnalysisClusterCountsAdaptor");
    }

    #printf STDERR sprintf("\n\noPOSSUM State:\n%s\n\n",
    #    Data::Dumper::Dumper($self->state())
    #);
    
    #
    # fetch target and background conserved region GC content
    #
    my $t_cr_gc_content = 0;
    #$t_cr_gc_content = $self->fetch_cr_gc_content(
    #    $state->t_gids, $state->conservation_level,
    #    $state->upstream_bp, $state->downstream_bp,
	#	$state->biotype
    #);
    
    my $bg_cr_gc_content = 0;
    #$bg_cr_gc_content = $self->fetch_cr_gc_content(
    #    $state->bg_gids, $state->conservation_level,
    #    $state->upstream_bp, $state->downstream_bp,
	#	$state->biotype
    #);

    # if operon genes present, the retrieved gene counts are all based on
    # the first gene search region, taken care by the CountsAdaptor.
    # no further action necessary.    
    my $t_counts;
    my $bg_counts;
    my $t_cr_length;
    my $bg_cr_length;
    my $crla = $opdba->get_ConservedRegionLengthAdaptor();
    if ($analysis_type eq 'default') {
        $t_counts = $self->fetch_cluster_analysis_counts(
            $analysis_type,
            -gene_ids               => $state->t_gids(),
            -has_operon             => $state->has_operon(), #optional
			-operon_gene_ids		=> $state->t_operon_gids(),
            -cluster_ids            => $state->tf_cluster_set->cluster_ids(),
            -conservation_level     => $state->conservation_level(),
            -threshold_level        => $state->threshold_level(),
            -search_region_level    => $state->search_region_level()
        );

        return $self->error("Error fetching target gene TFBS cluster counts")
            unless $t_counts;

        $bg_counts = $self->fetch_cluster_analysis_counts(
            $analysis_type,
            -gene_ids               => $state->bg_gids(),
            -has_operon             => $state->has_operon(), #not optional, if choosing all genes
			-operon_gene_ids		=> $state->bg_operon_gids(),
            -cluster_ids               => $state->tf_cluster_set->cluster_ids(),
            -conservation_level     => $state->conservation_level(),
            -threshold_level        => $state->threshold_level(),
            -search_region_level    => $state->search_region_level()
        );

        return $self->error("Error fetching background gene TFBS cluster counts")
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
		
        $t_counts = $self->fetch_cluster_analysis_counts(
            $analysis_type,
            -gene_ids           => $state->t_gids(),
            -has_operon         => $state->has_operon(), # optional
			-operon_gene_ids	=> $state->t_operon_gids(),
            -clusters        	=> $state->tf_cluster_set->get_tf_cluster_list(),
            -conservation_level => $state->conservation_level(),
            -threshold          => $state->threshold(),
            -upstream_bp        => $state->upstream_bp(),
            -downstream_bp      => $state->downstream_bp()
        );

        return $self->error("Error fetching target gene TFBS cluster counts")
            unless $t_counts;

		my $t_counts_gids = $t_counts->gene_ids;
		my $t_counts_cids = $t_counts->cluster_ids;
		#print STDERR "\n";
		#print STDERR "\n\nt_counts:\n" . Data::Dumper::Dumper($t_counts) . "\n\n";
		
        $bg_counts = $self->fetch_cluster_analysis_counts(
            $analysis_type,
            -gene_ids           => $state->bg_gids(),
            -has_operon         => $state->has_operon(), # not optional
			-operon_gene_ids	=> $state->t_operon_gids(),
            -clusters        	=> $state->tf_cluster_set->get_tf_cluster_list(),
            -conservation_level => $state->conservation_level(),
            -threshold          => $state->threshold(),
            -upstream_bp        => $state->upstream_bp(),
            -downstream_bp      => $state->downstream_bp()
        );

        return $self->error("Error fetching background gene TFBS cluster counts")
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

    my $fisher = OPOSSUM::Analysis::Cluster::Fisher->new();
    return $self->error("Error initializing Fisher analysis") unless $fisher;

    my $fresult_set = $fisher->calculate_Fisher_probability(
        $bg_counts,
        $t_counts
    );
    return $self->error("Error performing Fisher analysis")
        unless $fresult_set;

    my $zscore = OPOSSUM::Analysis::Cluster::Zscore->new();
    return $self->error("Error initializing z-score analysis") unless $zscore;
    
    # fetch cluster widths for z-score analysis
    #my $cluster_ids = $state->tf_cluster_set->cluster_ids;
    #my %cluster_avg_widths;
    #foreach my $cluster_id (@$cluster_ids) {
    #	my $cluster = $tf_cluster_set->get_tf_cluster($cluster_id);
    #	$cluster_avg_widths{$cluster_id} = $cluster->avg_width;
    #}
    
    my $zresult_set = $zscore->calculate_Zscore(
        $bg_counts,
        $t_counts,
        $bg_cr_length,
        $t_cr_length
        #\%cluster_avg_widths
    );
    return $self->error("Error computing z-score") unless $zresult_set;

    my $result_set = OPOSSUM::Analysis::Cluster::CombinedResultSet->new(
        -fisher_result_set  => $fresult_set,
        -zscore_result_set  => $zresult_set
    );

    return $self->error("Error combining Fisher and z-score result_set")
        unless $result_set;

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
		$self->write_cluster_results($abs_result_file, $results, $tf_cluster_set);
	}
	
    #print STDERR "results self:\n" . Data::Dumper::Dumper($self);

    # passed through to cluster_details routine
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
		rel_htdocs_data_path => REL_HTDOCS_DATA_PATH,
        result_file		    => RESULTS_TEXT_FILENAME,
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
        cl_select_criteria  => $state->cl_select_criteria(),
        #collection          => $state->collection(),
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
        tf_cluster_set      => $tf_cluster_set,
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
        var_template        => "results_gene_tca.html"
    };

    #print STDERR "results vars:\n" . Data::Dumper::Dumper($vars);

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

# for the given tfbs cluster, finds the list of genes with cluster hits and enumerate
# the hits
sub tfbs_cluster_details
{
    my $self = shift;

    #print STDERR "tfbs cluster_details\n";
    #print STDERR "cluster_details self:\n" . Data::Dumper::Dumper($self);

    my $q = $self->query;

    my $cluster_id = $q->param('cluster_id');

    #print STDERR "cluster_details query\n" . Data::Dumper::Dumper($q);

    my $state = $self->state();

    my $species             = $state->species();
    #my $collection          = $state->collection();
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
    my $t_operon_gids       = $state->t_operon_gids();
    my $gid_gene_ids        = $state->t_gid_gene_ids();
    #my $gene_id_gids        = $state->t_gene_id_gids();
    #my $gid_ensids          = $state->t_gid_ensids();
    my $tf_cluster_set      = $state->tf_cluster_set();
	my $results_subdir      = $state->results_subdir();

    my $cl = $tf_cluster_set->get_tf_cluster($cluster_id);

    # get opossum db connectors
    my $opdba = $self->opdba();
    unless ($opdba) {
        $opdba = $self->opossum_db_connect();
    }
	
    my $oa = $opdba->get_OperonAdaptor();
    my $ga = $opdba->get_GeneAdaptor();
    my $ctfsa = $opdba->get_ConservedTFBSAdaptor();
		
    my %gid_cluster_tfbss; # key = gid, val = list of cluster tfbs hits
    my @cl_genes;
	my %operon_genes; # key = first gene id, val = operon genes
	
    #
    # for operon genes, key = gid, val = first gene gid
    # t_gids contain trailing operon genes as well as first genes
    # tf_genes will be passed to write_tfbs_details
    # for operon genes, tf_genes should contain only the first genes
    # operon_genes, keyed by the first gene id, contains all operon genes
    # that are in t_gids. If the first gene is not in t_gids, it will not be
    # included. (but key = still first gene id)
    # 
	
    foreach my $gid (@$t_gids)
    {
        my $gene = $ga->fetch_by_id($gid);
        my $promoters;
		my $fgene;
		
        # if operon gene, get the first gene 		
        if ($$t_operon_gids{$gid}) {
            $fgene = $ga->fetch_by_id($$t_operon_gids{$gid});
        } else {
            $fgene = $gene;
        }
		
        # cannot search by cluster id with ctfsa
        # get sites for each tf in the cluster
        # merge all sites, and remove overlapping ones?
        # or just list all the member tfs?
        # or do both?
        # show the merged sites. otherwise it's too many of the same thing.
        my $tfids = $cl->tf_ids;
        my $tfset = OPOSSUM::ConservedTFBSSet->new();
        my $sites = $ctfsa->fetch_by_gene(
                -gene_id            => $fgene->id,
                -tf_ids             => $tfids,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
        );
        
        if (!$sites or scalar(@$sites) == 0) {
            $gid_cluster_tfbss{$gid} = $sites;
            next;
        }
        
        #print STDERR "# sites for $cluster_id = " . scalar(@$sites) . "\n";
        $tfset->add_tf_site_list($sites);
        #print STDERR "tfset size = " . $tfset->size . "\n";

        my $tfsites = $tfset->tf_sites('start'); # just to be safe
        my $merged_sites = merge_cluster_sites($tfsites, $cluster_id);
        
        push @cl_genes, $fgene; # first gene -> cl_genes
        $gid_cluster_tfbss{$gid} = $merged_sites;
		push @{$operon_genes{$fgene->id}}, $gene;
	#$gid_cluster_tfbss{$gid} = $tfsites;
    }

    #print STDERR "cluster_details tfbs genes:\n"
    #   . Data::Dumper::Dumper(@tf_genes);

    # cl_genes: list of genes with the cluster
    # key = gid, val = individual sites
    # so for this to work, I should have merged the sites within cluster by now
    # so when you call write_cluster_details, it goes through the list of cl_genes
    # and for each gene in the list, it goes through the sites listed and print them

    my $result_file = "Cluster_" . $cl->id . ".txt";
	
	my $abs_result_file = sprintf("%s/%s/%s.txt",
		ABS_HTDOCS_RESULTS_PATH,
		$results_subdir,
		$result_file
	);
	
	$self->write_cluster_details(
		$abs_result_file,
        $cl,
        \@cl_genes,
        \%gid_cluster_tfbss,
        $gene_id_type,
        $gid_gene_ids,
		$state->has_operon(),
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
        section             => 'TFBS Cluster Details',
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        dflt_gene_id_type   => DFLT_GENE_ID_TYPE,
        result_retain_days  => REMOVE_TEMPFILES_OLDER_THAN,
        tf_cluster          => $cl,
        genes               => \@cl_genes,
        gid_cluster_tfbss   => \%gid_cluster_tfbss,
        gene_id_type        => $gene_id_type,
        gid_gene_ids        => $gid_gene_ids,
        #gene_id_gids        => $gene_id_gids,
		has_operon			=> $state->has_operon(),
        operon_genes         => \%operon_genes,
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

    my $output = $self->process_template('tfbs_cluster_details_gene_tca.html', $vars);
    #print STDERR "cluster_details output\n" . Data::Dumper::Dumper($output);
    return $output;
}

# genes would be those that contain the specific cluster hits (cl)
# gid_clusters = ?
sub write_cluster_details
{
    my ($self, $filename, $cl, $genes, $gid_cluster_tfbss, $gene_id_type,
        $gid_gene_ids, $has_operon, $operon_genes) = @_;

    my $text = sprintf("%s\n\n", $cl->name());
    $text .= sprintf("Class:\t%s\n", $cl->class());
    $text .= sprintf("Family:\t%s\n", $cl->family());

    $text .= "\n\nConserved Binding Sites:\n\n";
    
    unless (defined $gene_id_type && $gene_id_type == DFLT_GENE_ID_TYPE) {
        $text .= "Gene ID(s)";
    }

    $text .= "\tEnsembl ID";

    if ($has_operon) {
        $text .= "\tOperon ID";
    }

    $text .= qq{\tChr\tGene Start\tGene End\tStrand\tCluster Start\tCluster End\tNearest TSS\tCluster Rel. Start\tCluster Rel. End\tCluster Strand\tTop Rel. Score\tRep. TFBS Sequence\n};
    
	#
    # Notes for operon genes
    # $genes -> only contain first operon genes. However, not all first genes may
    #   have been included in the input gene set. These can be checked with
    # $operon_genes -> contain all operon genes in the input set
    # So, whenever a gene in $genes maps to operon_genes, must get the corresponding
    # gene array from operon_genes and list them
    #
    foreach my $gene (@$genes)
    {
        my $gid = $gene->id;

        my $ensembl_id  = $gene->ensembl_id();
        my $gene_start  = $gene->start();
        my $gene_end    = $gene->end();
        my $strand      = $gene->strand();
        my $chr         = $gene->chr();
        my $tss         = $gene->tss();

        my $promoters   = $gene->fetch_promoters();

        my $prom_start;
        my $prom_end;
        if ($strand == 1) {
            $prom_start = $tss;
            $prom_end   = $gene_end;
        } else {
            $prom_start = $gene_start;
            $prom_end   = $tss;
        }

        if ($has_operon and defined $$operon_genes{$gid}) {
            unless (defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE)
            {
                my $count = 0;
                foreach my $in_gene (@{$$operon_genes{$gid}}){
                    $text .= ',' if $count > 0;
                    $text .= join(',', @{$gid_gene_ids->{$in_gene->id}});
                    $count++;
                }
                $text .= "\t";
            }
            my @ensembl_ids;
            foreach my $in_gene (@{$$operon_genes{$gid}}){
                push @ensembl_ids, $in_gene->ensembl_id();
            }
            $text .= join(',', @ensembl_ids) . "\t";
		
		} else {
			
            unless (defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE)
            {
                $text .= join(',', @{$gid_gene_ids->{$gid}}) . "\t";
            }        
            $ensembl_id  = $gene->ensembl_id();
            $strand      = $gene->strand();
            $chr         = $gene->chr();
            $tss         = $gene->tss(); # this tss is not used.
                                         # overwritten later.
            $text .= "$ensembl_id\t";
		}
		
        if ($has_operon) {
            my $symbol = "-";
            if (defined $gene->operon()) {
                $symbol = $gene->operon->symbol();
            }
            $text .= sprintf("\t%s", $symbol);
        }

        $text .= sprintf("\t%s\t%d\t%d\t%s",
            $chr,
            $prom_start,
            $prom_end,
            $strand == 1 ? '+' : '-'
        );

        my $tfbs_sites = $$gid_cluster_tfbss{$gid}; # sorted by start		
        if (!$tfbs_sites or scalar(@$tfbs_sites) == 0) {
                #$text .= "Zilch for $gid\n";
                $text .= "\n";
                next;
        }
        #print STDERR Data::Dumper::Dumper($tfbs_sites) ."\n";
        
		# gene_start/gene_end should refer to the first gene here
        # what about tss? which one should I use? show both?
        # actually, no. closest tss refers to the actual promoter tss
        # since the first gene tss is used, should be kept.
        
		my $first = 1;
        my $prev_site = $$tfbs_sites[0];
        foreach my $tfbs (@$tfbs_sites) {
            my $site_start      = $gene_start + $tfbs->start() - 1;
            my $site_end        = $gene_start + $tfbs->end() - 1;
            my $site_seq        = $tfbs->seq();
            my $site_strand     = $tfbs->strand();
            my $site_score      = $tfbs->score();
            my $site_rel_score  = $tfbs->rel_score();

            my $closest_tss;
            my $min_tss_dist = 999999;
            foreach my $promoter (@$promoters) {
                my $tss = $promoter->tss();

                my $start_tss_dist = abs($site_start - $tss);
                my $end_tss_dist   = abs($site_end - $tss);

                if ($start_tss_dist < $min_tss_dist) {
                    $min_tss_dist = $start_tss_dist;
                    $closest_tss = $tss;
                }

                if ($end_tss_dist < $min_tss_dist) {
                    $min_tss_dist = $end_tss_dist;
                    $closest_tss = $tss;
                }
            }

            my ($rel_start, $rel_end);
            if ($strand == 1) {
                $rel_start = $site_start - $closest_tss;
                if ($site_start >= $closest_tss) {
                    $rel_start++;
                }

                $rel_end = $site_end - $closest_tss;
                if ($site_end >= $closest_tss) {
                    $rel_end++;
                }
            } else {
                $rel_start = $closest_tss - $site_start;
                if ($site_start <= $closest_tss) {
                    $rel_start++;
                }

                $rel_end = $closest_tss - $site_end;
                if ($site_end <= $closest_tss) {
                    $rel_end++;
                }

                # swap coords so start is more upstream than end
                ($rel_start, $rel_end) = ($rel_end, $rel_start);
            }

            unless ($first) {
                
                unless (overlap($prev_site, $tfbs)) {
                    $text .= "-------------------------\n";
                }
                
                $text .= "\t\t\t\t\t";

                unless ($gene_id_type == DFLT_GENE_ID_TYPE) {
                    $text .= "\t";
                }
            }
            
            $text .= sprintf("\t%d\t%d\t%d\t%d\t%d\t%s\t%.3f (%.1f%%)\t%s\n",
                $site_start,
                $site_end,
                $closest_tss,
                $rel_start,
                $rel_end,
                $site_strand == 1 ? '+' : '-',
                $site_score,
                $site_rel_score * 100,
                $site_seq
            );

            $first = 0;
            $prev_site = $tfbs;
        } # end foreach tfbs
    } # end foreach gene

    unless (open(FH, ">$filename")) {
        $self->_error("Unable to create TFBS details file $filename - $!");
        return;
    }

    print FH $text;
    close(FH);

    $filename =~ s/.*\///;

    return $filename;
}

sub text_results
{
    my $self = shift;

    my $q = $self->query;

    my $state = $self->state();

    my $species = $state->species();

    my $result_text_file = $q->param('result_cluster_text_file');

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

sub text_tfbs_cluster_details
{
    my $self = shift;

    my $q = $self->query;

    my $state = $self->state();

    my $species = $state->species();

    my $cl_gene_text = $q->param('tf_cluster_gene_text');

    my $vars = {
        abs_htdocs_path    => ABS_HTDOCS_PATH,
        rel_htdocs_path    => REL_HTDOCS_PATH,
        abs_cgi_bin_path   => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path   => REL_CGI_BIN_PATH,
        bg_color_class     => $state->bg_color_class(),
        title              => $state->title(),
        heading            => $state->heading(),
        section            => 'TFBS Cluster Details',
        version            => VERSION,
        devel_version      => DEVEL_VERSION,
        text               => $cl_gene_text
    };

    my $output = $self->process_template('blank.txt', $vars);
    #my $template = Template->new();
    #my $input = "[%text%]";

    #my $output;
    #$template->process(\$input, $vars, \$output);
    #print STDERR "text_cluster_details:\n$string\n";

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

    my $heading = sprintf "%s TFBS Cluster Analysis",
        ucfirst $state->species();
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    $state->errors(undef);
    $state->warnings(undef);
}

1;
