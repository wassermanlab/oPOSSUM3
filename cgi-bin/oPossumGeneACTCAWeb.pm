package oPossumGeneACTCAWeb;

use base 'CGI::Application';

use oPossumWebInclude;
use oPossumGeneWebInclude;
use oPossumTCAWebInclude;
use oPossumACSAWebOpt;
use oPossumTCAWebOpt;

use Data::Dumper;    # for debugging only
use File::Temp qw/ tempdir /;

use OPOSSUM::Web::State;
#use OPOSSUM::Analysis::Cluster::Fisher;
#use OPOSSUM::Analysis::Cluster::Zscore;
#use OPOSSUM::Analysis::Cluster::CombinedResultSet;

use CGI::Carp qw(carpout);    # fatalsToBrowser;

use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_gene_actca';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_gene_actca";
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
        'process'           => 'process'
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
    $self->_clean_resultfiles;
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
    my $biotype = 'protein_coding' if $species eq 'worm';
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

    my $tax_group       = $db_info->tax_group();
    my $min_ic          = $db_info->min_ic();
    my $max_upstream_bp = $db_info->max_upstream_bp();

    # for worms, tax_group include nematodes,vertebrates,insects
    # you could give this option for insects as well
    # in the input menu, you need to be able to choose from 3 tax groups
    # easier to deal with if you separate them
    my @tax_groups      = split /\s*,\s*/, $tax_group;
    my $num_tax_groups  = scalar @tax_groups;

    $state->has_operon($db_info->has_operon());
    $state->tax_groups(\@tax_groups);
    $state->num_tax_groups($num_tax_groups);

    #
    # When these are stored in state and retrieved later in the process()
    # routine, they seem to lose their association with the proper
    # class, and method calls on them fail!
    # e.g. An OPOSSUM::ConservationLevel object retrieved from state
    # is still recognized as such but calling min_conservation on it fails.
    #
    $state->conservation_level_hash($cl_hash);
    $state->threshold_level_hash($thl_hash);
    $state->search_region_level_hash($srl_hash);

    #
    # Connect to oPOSSUM_cluster DB and retrieve TFCluster info
    #
    $self->opossum_cluster_db_connect();
    $self->jaspar_db_connect();

    #printf STDERR "tax_groups: %s\n", join ',', @tax_groups;

    # Have the user select an anchor TF
    # then search using the TFCluster that the anchor TF belongs to
    my $tf_set = $self->fetch_tf_set(
        -collection => 'CORE',
        -tax_group  => \@tax_groups,
        -min_ic     => $min_ic
    );
    my $tf_cluster_set = $self->fetch_tf_cluster_set();

    #printf STDERR "tf_set:\n%s\n", Data::Dumper::Dumper($tf_set);

    #my %tax_tf_sets;
    #foreach my $tax_group (@tax_groups) {
    #    $tax_tf_sets{$tax_group} = $tf_set->subset(-tax_groups => $tax_group);
    #}

    #printf STDERR "tax_tf_sets:\n%s\n", Data::Dumper::Dumper(%tax_tf_sets);

    #printf STDERR "\ntf_set:\n"
    #    . Data::Dumper::Dumper($tf_set) . "\n\n";

    #printf STDERR "input_gene_ids:\n" . Data::Dumper::Dumper(\@input_gene_ids);

    #
    # Format a tax group list for display on the web page. The $tax_group
    # variable is obtained from the db and may be a single value or a comma
    # separated list (no space after comma).
    #
    my @s_tax_groups;
    foreach my $tg (@tax_groups) {
        # remove trailing 's' (un-pluralize)
        my $stg = $tg;
        $stg =~ s/s$//;
        push @s_tax_groups, $stg;
    }

    my $tax_group_list = join ' / ', @s_tax_groups;

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
        dflt_inter_binding_dist => DFLT_INTER_BINDING_DIST,
        max_inter_binding_dist  => MAX_INTER_BINDING_DIST,
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
        tax_group_list          => $tax_group_list,
        tf_set                  => $tf_set,
        tf_cluster_set          => $tf_cluster_set,
        #tax_tf_sets             => \%tax_tf_sets,
        in_t_gene_id_type       => $in_t_gene_id_type,
        in_t_gene_ids           => \@in_t_gene_ids,
        var_template            => "input_gene_actca.html"
    };

    my $output = $self->process_template('master.html', $vars);
    #print STDERR "input results:\n"
    #		    . Data::Dumper::Dumper($output);

    return $output;
}

sub process
{
    my $self = shift;

    my $q = $self->query;

    my $state = $self->state();

    my $tf_db = JASPAR_DB_NAME;
    my $cl_db = TFBS_CLUSTER_DB_NAME;

    #my $tf_collection = 'CORE';
    #$state->collection($tf_collection);
    my $dflt_collections = ['CORE','PBM','PENDING'];
    $state->collections($dflt_collections);

    my $species = $state->species() || $self->param('species');

    if (!defined $species) {
        return $self->error("Species not specified");
    }

    my $opdba = $self->opdba();

    my $dbia = $opdba->get_DBInfoAdaptor();
    my $db_info = $dbia->fetch_db_info();

    #my $tax_group = $db_info->tax_group();
    my $dflt_min_ic = $db_info->min_ic();


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
    if ($bg_id_input_method eq 'paste' || $bg_id_input_method eq 'upload') {
        $bg_gene_id_type = $q->param("bg_gene_id_type");
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

    #printf STDERR "\nt_gene_text:\n$t_gene_text\n\n";

    my $t_gene_ids = $self->parse_gene_id_text($t_gene_text);

    #printf STDERR "\nt_gene_ids:\n%s\n\n", Data::Dumper::Dumper($t_gene_ids);

    return $self->error("Error parsing target gene IDs") unless $t_gene_ids;

    $state->t_gene_ids($t_gene_ids);

    my $bg_gene_ids;
    my $bg_num_rand_genes;
    if ($bg_id_input_method eq 'paste') {
        my $bg_gene_text = $q->param('bg_gene_text');

        return $self->error("No background gene IDs pasted")
            unless $bg_gene_text;

        $bg_gene_ids = $self->parse_gene_id_text($bg_gene_text);

        return $self->error("Error parsing pasted background gene IDs")
            unless $bg_gene_ids;
    } elsif ($bg_id_input_method eq 'upload') {
        my $gene_file = $q->param('bg_gene_file');
        my $gene_fh   = $q->upload('bg_gene_file');

        my $bg_gene_text;
        while (my $line = <$gene_fh>) {
            $bg_gene_text .= $line;
        }

        return $self->error("No background gene IDs uploaded")
            unless $bg_gene_text;

        $bg_gene_ids = $self->parse_gene_id_text($bg_gene_text);

        return $self->error("Error parsing uploaded background gene IDs")
            unless $bg_gene_ids;
    } elsif ($bg_id_input_method eq 'random') {
        $bg_num_rand_genes = $q->param('bg_num_rand_genes');
    }

    $state->bg_gene_ids($bg_gene_ids);
    $state->bg_num_rand_genes($bg_num_rand_genes);

    #printf STDERR "\nbg oPOSSUM gene IDs:\n%s\n\n",
    #    Data::Dumper::Dumper($bg_gids) if $bg_gids;

    my $anchor_tf_id = $q->param('anchor_tf_id');
    if (!$anchor_tf_id) {
        return $self->error("Anchoring TF not specified");
    }
    $state->anchor_tf_id($anchor_tf_id);

    my $tf_cluster_select_method = $q->param("tf_cluster_select_method");
    if (!$tf_cluster_select_method) {
        return $self->error("TF families not specified");
    }
    $state->tf_cluster_select_method($tf_cluster_select_method);

    my $tf_cluster_select_criteria;
    if ($tf_cluster_select_method =~ /specific/) {
        $tf_cluster_select_criteria = 'specific';
    } else {
        $tf_cluster_select_criteria = 'all';
    }
    $state->tf_cluster_select_criteria($tf_cluster_select_criteria);

    #my $min_ic;
    #my @tax_groups;
    my @tf_families;
    my $tf_families_str;
    #my $tax_group_str;
    #my @tf_ids;
    if ($tf_cluster_select_criteria eq 'specific') {
        foreach my $fam ($q->param('tf_families')) {
            push @tf_families, $fam;
        }

        if (!@tf_families) {
            return $self->error("No specific TF families selected");
        }
    }
    $tf_families_str = join(',', @tf_families);

    my $conservation_level = $q->param('conservation_level');
    printf STDERR "conservation level: $conservation_level\n";
    $state->conservation_level($conservation_level);

    #
    # For some reason methods retrieved from state lose their method
    # association!
    #
    #my $cl_hash = $state->conservation_level_hash();
    my $cl_hash = $self->fetch_conservation_levels();
    #printf STDERR "cl hash:\n%s\n\n", Data::Dumper::Dumper($cl_hash);

    #my $cl = $cl_hash->{$conservation_level};
    #printf STDERR "cl:\n%s\n\n", Data::Dumper::Dumper($cl);
    #printf STDERR "cl is an OPOSSUM::ConservationLevel: %s\n\n",
    #    $cl->isa("OPOSSUM::ConservationLevel");
    #printf STDERR "cl can do min_conservation: %s\n\n",
    #    $cl->can("min_conservation");

    my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();
    $state->min_conservation($min_conservation);

    my $threshold_level = $q->param('threshold_level');
    my $threshold       = $q->param('threshold');

    #my $thl_hash = $state->threshold_level_hash();
    my $thl_hash = $self->fetch_threshold_levels();

    if ($threshold) {
        $threshold_level = undef;

        if ($threshold =~ /(.+)%$/) {
            $threshold = $1;
        }
        $threshold /= 100;
    } else {
        $threshold = $thl_hash->{$threshold_level}->threshold();
    }
    $state->threshold_level($threshold_level);
    $state->threshold($threshold);

    my $search_region_level = $q->param('search_region_level');
    my $upstream_bp         = $q->param('upstream_bp');
    my $downstream_bp       = $q->param('downstream_bp');

    #my $srl_hash = $state->search_region_level_hash();
    my $srl_hash = $self->fetch_search_region_levels();

    if (   (defined $upstream_bp && $upstream_bp ne '' && $upstream_bp > 0)
        || (defined $downstream_bp && $downstream_bp ne ''
            && $downstream_bp > 0)
    ) {
        #
        # Custom search region entered
        #
        $search_region_level = undef;

        #
        # Check if entered search region corresponds to a default search
        # region level
        #
        #foreach my $level (keys %$srl_hash) {
        #    if ($species eq 'yeast') {
        #        if ($upstream_bp == $srl_hash->{$level}->upstream_bp()) {
        #            $search_region_level = $level;
        #            last;
        #        }
        #    } else {
        #        if (   $upstream_bp == $srl_hash->{$level}->upstream_bp()
        #            && $downstream_bp == $srl_hash->{$level}->downstream_bp()
        #        ) {
        #            $search_region_level = $level;
        #            last;
        #        }
        #    }
        #}
    } else {
        $upstream_bp    = $srl_hash->{$search_region_level}->upstream_bp();
        $downstream_bp  = $srl_hash->{$search_region_level}->downstream_bp();
    }
    $state->search_region_level($search_region_level);
    $state->upstream_bp($upstream_bp);
    $state->downstream_bp($downstream_bp);

    my $inter_binding_dist = $q->param('inter_binding_dist');
    if ($inter_binding_dist > MAX_INTER_BINDING_DIST) {
        return $self->error(
            sprintf(
                "Specified inter-binding distance exceeds maximum of %d allowed",
                MAX_INTER_BINDING_DIST
            )
        );
    }

    my $num_display_results;
    my $zscore_cutoff;
    my $fisher_cutoff;
    my $result_type = $q->param('result_type');
    $state->result_type($result_type);
    if ($result_type eq 'top_x_results') {
        $num_display_results = $q->param('num_display_results');
    } elsif ($result_type eq 'significant_hits') {
        $zscore_cutoff = $q->param('zscore_cutoff');
        $fisher_cutoff = $q->param('fisher_cutoff');
    }

    $state->num_display_results($num_display_results);
    $state->zscore_cutoff($zscore_cutoff);
    $state->fisher_cutoff($fisher_cutoff);

    my $result_sort_by = $q->param('result_sort_by');
    $state->result_sort_by($result_sort_by);

    my $email = $q->param('email');
    $state->email($email);

    my $user_t_gene_file;
    if ($t_id_input_method eq 'upload') {
        $user_t_gene_file = $q->param('t_gene_file');

        $user_t_gene_file =~ s/.*\///;
    }

    my $user_bg_gene_file;
    if ($bg_id_input_method eq 'upload') {
        $user_bg_gene_file = $q->param('bg_gene_file');

        $user_bg_gene_file =~ s/.*\///;
    }

    #
    # Create a temporary working directory for all the temp. input files and
    # output results file as a sub-directory under the defined temp. dir.
    #
    my $tempdir = tempdir(DIR => ABS_HTDOCS_RESULTS_PATH);

    my $job_id = $tempdir;
    $job_id =~ s/.*\///;

    my $t_gene_ids_file = $self->create_local_working_file(
        $tempdir, 't_gene_ids', $t_gene_ids
    );

    unless ($t_gene_ids_file) {
        return $self->error("Could not create target gene IDs file");
    }

    my $bg_gene_ids_file;
    if ($bg_gene_ids) {
        $bg_gene_ids_file = $self->create_local_working_file(
            $tempdir, 'bg_gene_ids', $bg_gene_ids
        );

        unless ($bg_gene_ids_file) {
            return $self->error("Could not create background gene IDs file");
        }
    }

    #my $tf_ids_file;
    #if ($tf_select_criteria eq 'specific') {
    #    $tf_ids_file = $self->create_local_working_file(
    #        $tempdir, 'tf_ids', \@tf_ids
    #    )
    #}

    my $command = "./opossum_gene_actca.pl -j $job_id -d $tempdir -s $species"
        . " -g $t_gene_ids_file -tdb $tf_db -cdb $cl_db -id $anchor_tf_id"
        . " -cl $conservation_level -dist $inter_binding_dist";

    if ($bg_gene_ids_file) {
        $command .= " -b $bg_gene_ids_file";
    } elsif ($bg_num_rand_genes) {
        $command .= " -bnr $bg_num_rand_genes";
    }

    if (defined $t_gene_id_type) {
        $command .= " -gt $t_gene_id_type";
    }

    if (defined $bg_gene_id_type) {
        $command .= " -bt $bg_gene_id_type";
    }

    if ($tf_families_str) {
        $command .= " -fam '$tf_families_str'";
    }

    if (defined $threshold) {
        $command .= " -th $threshold";
    }

    if (defined $upstream_bp) {
        $command .= " -up $upstream_bp";
    }

    if (defined $downstream_bp) {
        $command .= " -dn $downstream_bp";
    }

    if (defined $num_display_results) {
        $command .= " -n $num_display_results";
    } else {
        $command .= " -zcutoff $zscore_cutoff";
        $command .= " -fcutoff $fisher_cutoff";
    }

    if ($result_sort_by) {
        $command .= " -sr $result_sort_by";
    }

    if ($email) {
        $command .= " -m $email";
    }

    if ($state->has_operon()) {
        $command .= " -has_operon";
    }

    my $biotype = $state->biotype();
    if ($biotype) {
        if (ref $biotype eq 'ARRAY') {
            $command .= sprintf " -biotype '%s'", join ',', @$biotype;
        } else {
            $command .= sprintf " -biotype '%s'", $biotype;
        }
    }

    printf LOG "Starting ACTCA analysis at %s:\n$command\n\n",
        scalar localtime(time);

    system("$command >/dev/null 2>&1 &");

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        bg_color_class      => BG_COLOR_CLASS,
        rel_htdocs_results_path => REL_HTDOCS_RESULTS_PATH,
        rel_htdocs_data_path => REL_HTDOCS_DATA_PATH,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        heading             => $state->heading(),
        title               => $state->title(),
        section             => 'Analysis Submitted',
        jaspar_url          => JASPAR_URL,
        result_retain_days  => REMOVE_TEMPFILES_OLDER_THAN,
        sid                 => $state->sid(),
        job_id              => $job_id,
        user_t_gene_file    => $user_t_gene_file,
        user_bg_gene_file   => $user_bg_gene_file,
        species             => $species,
        tf_cluster_select_criteria  => $state->tf_cluster_select_criteria(),
        t_gene_id_type      => $state->t_gene_id_type(),
        bg_gene_id_type     => $state->bg_gene_id_type(),
        tf_cluster_select_criteria => $tf_cluster_select_criteria,
        #collection          => $state->collection(),
        #tf_families_str     => $tf_families_str,
        tf_families         => \@tf_families,
        #min_ic              => $state->min_ic(),
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
        gene_ids            => $state->t_gene_ids(),
        num_genes           => scalar @{$state->t_gene_ids()},
        email               => $email,
        var_template        => "analysis_summary_gene_actca.html"
    };

    #print STDERR "results vars:\n" . Data::Dumper::Dumper($vars);

    my $output = $self->process_template('master.html', $vars);

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

    my $heading = sprintf "%s Anchored Combination TFBS Cluster Analysis",
        ucfirst $state->species();
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    $state->errors(undef);
    $state->warnings(undef);
}

1;
