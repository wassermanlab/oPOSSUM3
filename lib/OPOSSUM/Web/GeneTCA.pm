package GeneTCA;

use base 'CGI::Application';

use OPOSSUM::Web::Include::BaseInclude;
use OPOSSUM::Web::Include::GeneInclude;
use OPOSSUM::Web::Include::TCAInclude;

use Data::Dumper;    # for debugging only
use Template;
use CGI;
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
$log_file .= "_$USER" if $USER;
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
        'input'          => 'input',
        'process'        => 'process',
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

    $self->errors(undef);
    $self->warnings(undef);

    #printf STDERR sprintf("\n\noPOSSUM State:\n%s\n\n",
    #    Data::Dumper::Dumper($self->state())
    #);

    #printf STDERR "\n\nrun mode = %s\n\n", $q->param('rm');
}


sub teardown
{
    my $self = shift;

    #print STDERR "teardown\n";

    my $state = $self->state();
    if ($state) {
        $state->dumper->Purity(1);
        $state->dumper->Deepcopy(1);
        $state->commit();
    }

    #$self->_clean_tempfiles;
    #$self->_clean_resultfiles;
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
    #my $in_t_gene_id_type = $state->t_gene_id_type();
    my @in_t_gene_ids = @{$state->t_gene_ids()} if $state->t_gene_ids();

    #
    # If gene IDs passed in via the query string from an external
    # application, these take precedent. These are always assumed
    # to be Ensembl IDs (gene ID type = 0).
    #
    if ($q->param("id")) {
        foreach my $id ($q->param("id")) {
            push @in_t_gene_ids, $id;
        }
        #$in_t_gene_id_type = DFLT_GENE_ID_TYPE;
    }

    # this should be later modified to give the user a choice
    my $biotype = DFLT_BIOTYPE;
    #my $biotype = 'protein_coding' if $species eq 'worm';
    #$biotype = _parse_biotype($biotype);
    #$state->biotype($biotype);
	
    my $total_genes;
    if ($biotype) {
        $total_genes = $self->fetch_gene_count("biotype = '$biotype'");
    } else {
        $total_genes = $self->fetch_gene_count();
    }
    my $db_info = $self->fetch_db_info();
    #my $xgid_types = $self->fetch_external_gene_id_types();
    my $cl_hash = $self->fetch_conservation_levels();
    my $thl_hash = $self->fetch_threshold_levels();
    my $srl_hash = $self->fetch_search_region_levels();

    my @cl_levels  = sort keys %$cl_hash;
    my @thl_levels = sort keys %$thl_hash;
    my @srl_levels = sort keys %$srl_hash;

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
	
    my $tax_group_str   = $db_info->tax_group();
    my $min_ic          = $db_info->min_ic();
    my $max_upstream_bp = $db_info->max_upstream_bp();
	
    # for worms, tax_group include nematodes,vertebrates,insects
    # you could give this option for insects as well
    # in the input menu, you need to be able to choose from 3 tax groups
    # easier to deal with if you separate them
    my @tax_groups      = split /\s*,\s*/, $tax_group_str;
    my $num_tax_groups  = scalar @tax_groups;


    $state->has_operon($db_info->has_operon());
    $state->tax_groups(\@tax_groups);
    $state->num_tax_groups($num_tax_groups);

    #
    # Format a tax group list for display on the web page. The $tax_group
    # variable is obtained from the db and may be a single value or a comma
    # separated list (no space after comma).
    #
    #my @s_tax_groups;
    #foreach my $tg (@tax_groups) {
        # remove trailing 's' (un-pluralize)
    #    my $stg = $tg;
    #    $stg =~ s/s$//;
    #    push @s_tax_groups, $stg;
    #}

    #my $tax_group_list = join ' / ', @s_tax_groups;

    #
    # Connect to oPOSSUM_cluster DB and retrieve TFCluster info
    # Connect to JASPAR DB and retrieve TF info
    #
	$self->jaspar_db_connect();
    $self->opossum_cluster_db_connect();

    my $tf_cluster_set = $self->fetch_tf_cluster_set();


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
        max_target_genes        => MAX_TARGET_GENES,
        dflt_bg_num_rand_genes  => DFLT_BG_NUM_RAND_GENES,
        sid                     => $state->sid(),
        species                 => $species,
        db_info                 => $db_info,
        total_genes             => $total_genes,
        #xgid_types              => $xgid_types,
        cl_hash                 => $cl_hash,
        thl_hash                => $thl_hash,
        srl_hash                => $srl_hash,
        cl_levels               => \@cl_levels,
        thl_levels              => \@thl_levels,
        srl_levels              => \@srl_levels,
		collections				=> JASPAR_COLLECTIONS_USED,
        tax_groups              => \@tax_groups,
        num_tax_groups          => $num_tax_groups,
        #tax_group_list          => $tax_group_list,
        tf_cluster_set          => $tf_cluster_set,
        #in_t_gene_id_type       => $in_t_gene_id_type,
        in_t_gene_ids           => \@in_t_gene_ids,
		min_ic					=> $min_ic,
        var_template            => "input_gene_tca.html"
    };

    my $output = $self->process_template('master.html', $vars);

    return $output;
}


sub process
{
    my $self = shift;

    my $state = $self->state();

    my $q = $self->query();

    my $species = $state->species() || $self->param('species');
    if (!defined $species) {
        return $self->error("Species not specified");
    }

    my $dflt_min_ic     = DFLT_CORE_MIN_IC;
    my $dflt_collections = JASPAR_COLLECTIONS_USED;
    my $dflt_tax_groups = JASPAR_CORE_TAX_GROUPS;

    my $t_id_input_method = $q->param('t_id_input_method');
    if (!$t_id_input_method) {
        return $self->error("No target gene IDs specified");
    }

	my $t_gene_id_type		= $q->param('t_gene_id_type');

    my $bg_id_input_method = $q->param('bg_id_input_method');
    if (!$bg_id_input_method) {
        return $self->error("Background gene ID input method not specified");
    }

    my $bg_gene_id_type;
    my $bg_num_rand_genes;
    if ($bg_id_input_method eq 'paste' || $bg_id_input_method eq 'upload') {
        $bg_gene_id_type = $q->param("bg_gene_id_type");
        #unless (defined $bg_gene_id_type) {
        #    return $self->error("Background gene ID type not specified");
        #}
    } elsif ($bg_id_input_method eq 'random') {
        $bg_gene_id_type = DFLT_GENE_ID_TYPE;
        $bg_num_rand_genes = $q->param('bg_num_rand_genes');
        unless ($bg_num_rand_genes) {
            return $self->error(
                "Number of random background genes not specified"
            );
        }
    }

    my $conservation_level;
    if ($species eq 'yeast') {
        $conservation_level = 1;
    } else {
        $conservation_level = $q->param('conservation_level');
        if (!$conservation_level) {
            $conservation_level = DFLT_CONSERVATION_LEVEL;
        }
    }
	my $threshold_level		= $q->param('threshold_level');
    my $threshold           = $q->param('threshold');
	my $search_region_level = $q->param('search_region_level');
	my $upstream_bp			= $q->param('upstream_bp');
	my $downstream_bp		= $q->param('downstream_bp');
    my $result_type         = $q->param('result_type');
    my $num_display_results = $q->param('num_display_results');
    my $zscore_cutoff       = $q->param('zscore_cutoff');
    my $fisher_cutoff       = $q->param('fisher_cutoff');
    my $result_sort_by      = $q->param('result_sort_by');
    my $email               = $q->param('email');    

    unless ($email) {
        return $self->error(
          "Please provide an e-mail address where your results will be mailed."
        );
    }

    #
    # For some reason methods retrieved from state lose their method
    # association!
	# Andrew: well, since we are going to be calling a separate script to
	# run the analysis anyways, don't really need to save in state.
    #
    my $cl_hash = $self->fetch_conservation_levels();
    my $thl_hash = $self->fetch_threshold_levels();
    my $srl_hash = $self->fetch_search_region_levels();

	#my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();
	#$state->min_conservation($min_conservation);
	#$state->conservation_level($conservation_level);

    #
    # If threshold, upstream/downstream bp are defined strip any blanks and
    # if result is empty string, set to undef as future code needs to
    # distinguish undefined values from empty strings. Note: values of 0 are OK.
    # Should really check that if these values are set, they are actually
    # numeric.
    #
    # DJA 2012/03/23
    #
    if (defined $upstream_bp) {
        $upstream_bp =~ s/\s+//g;

        if ($upstream_bp eq '') {
            $upstream_bp = undef;
        }
    }

    if (defined $downstream_bp) {
        $downstream_bp =~ s/\s+//g;

        if ($downstream_bp eq '') {
            $downstream_bp = undef;
        }
    }

    if (defined $threshold) {
        $threshold =~ s/\s+//g;

        if ($threshold eq '') {
            $threshold = undef;
        }
    }

    #
    # Convert to number in range 0-1. Rather than '%' as analysis script
    # will interpret correctly. This is also consistent with the way this
    # value would be retrieved from the threshold level hash.
    #
    # DJA 2012/03/07
    #
    if (defined $threshold) {
        $threshold /= 100;
    }

    #
    # If neither threshold or threshold level is set, then set threshold level
    # to default but DON'T set actual threshold as it messes up the logic
    # further down with regard to  which arguments are passed to analysis
    # script, e.g. if threshold is set then it is assumed we are doing a custom
    # analysis and threshold level is NOT passed to analysis script.
    #
    # NOTE: This should never happen as the threshold level will ALWAYS have
    # some value since it is pull down menu.
    #
    # DJA 2012/03/07
    #
    if (!defined $threshold && !$threshold_level) {
        $threshold_level = DFLT_THRESHOLD_LEVEL;
        #$threshold = $thl_hash->{$threshold_level}->threshold();
    }

	#$state->threshold($threshold);
	#$state->threshold_level($threshold_level);

    if ($species eq 'yeast') {
        if (!defined($upstream_bp) && !$search_region_level) {
            $search_region_level = DFLT_SEARCH_REGION_LEVEL;
        }
    } else {
        if (defined($upstream_bp) != defined($downstream_bp)) {
            return $self->error(
                "Both upstream and downstream search regions must be specified"
            );
        } elsif (!$upstream_bp && !$downstream_bp) {
            if (!$search_region_level) {
                $search_region_level = DFLT_SEARCH_REGION_LEVEL;
            }

            #
            # DON'T set upstream/downstream bp for the same reasons as outlined
            # for not setting the threshold when setting default threshold level
            # above.
            #
            #$upstream_bp = $srl_hash->{$search_region_level}->upstream_bp();
            #$downstream_bp = $srl_hash->{$search_region_level}->downstream_bp();
        }
    }

	#printf LOG ("Upstream_bp " . $upstream_bp);
	#printf LOG ("Downstream_bp " . $downstream_bp);
	#printf LOG ("SRL: " . $search_region_level);

	#$state->search_region_level($search_region_level);
	#$state->upstream_bp($upstream_bp);
	#$state->downstream_bp($downstream_bp);

    #
    # For display on analysis summary page
    #
    my $disp_min_conservation;
    my $disp_upstream_bp;
    my $disp_downstream_bp;
    my $disp_threshold;

	$disp_min_conservation
        = $cl_hash->{$conservation_level}->min_conservation();

    if (defined $upstream_bp) {
        $disp_upstream_bp = $upstream_bp;
    } else {
        $disp_upstream_bp
            = $srl_hash->{$search_region_level}->upstream_bp();
    }

    if (defined $downstream_bp) {
        $disp_downstream_bp = $downstream_bp;
    } else {
        $disp_downstream_bp
            = $srl_hash->{$search_region_level}->downstream_bp();
    }

    if (defined $threshold) {
        $disp_threshold = $threshold;
    } else {
        $disp_threshold
            = $thl_hash->{$threshold_level}->threshold();
    }
	
    my $tf_db = JASPAR_DB_NAME;
    my $cl_db = TFBS_CLUSTER_DB_NAME;
	
    my $min_ic;
    $min_ic = $q->param('tf_min_ic');
    
	#my @tax_groups;
    #my $tax_groups_str;
    #push @tax_groups, $q->param('tf_tax_groups');
    #if (!@tax_groups or scalar(@tax_groups == 0)) {
    #    push @tax_groups, @$dflt_tax_groups;
    #}
    #$tax_groups_str = join(',', @tax_groups);
	#my $num_tax_groups = scalar @tax_groups;
    
    #my @collections;
    #my $collections_str;   
    #push @collections, $q->param('tf_collections');
    #if (!@collections or scalar(@collections == 0)) {
        #push @collections, @$dflt_collections;
        #return $self->error("No JASPAR collection selected");
    #}
    #$collections_str = join(',', @collections);
    	
    my $tf_family_select_method = $q->param("tf_family_select_method");
    if (!$tf_family_select_method) {
        return $self->error("No TFBS family selected");
    }
    
    my $user_tf_family_file;
    if ($tf_family_select_method eq 'upload') {
        $user_tf_family_file = $q->param('tf_family_upload_file');
        $user_tf_family_file =~ s/.*\///;
    }
    
    my $user_t_gene_file;
    if ($t_id_input_method eq 'upload') {
        $user_t_gene_file = $q->param('t_gene_file');
        $user_t_gene_file =~ s/.*\///;
        $user_t_gene_file =~ s/\s+//g;
    }

    my $user_bg_gene_file;
    if ($bg_id_input_method eq 'upload') {
        $user_bg_gene_file = $q->param('bg_gene_file');
        $user_bg_gene_file =~ s/.*\///;
        $user_bg_gene_file =~ s/\s+//g;
    }

    #
    # Create a temporary working directory for all the input genes and
    # output results file as a sub-directory under the defined temp. dir.
    #
    my $tempdir = tempdir(DIR => ABS_HTDOCS_RESULTS_PATH);
    
    # for testing
    my $job_id = $tempdir;
    $job_id =~ s/.*\///;

    my $tf_family_filename;
    if ($tf_family_select_method !~ /all/)
    {        
        $tf_family_filename = $self->get_tf_family_file(
            $tf_family_select_method,
            $tempdir
        );
        
        return $self->error() if !$tf_family_filename;
    }

    my $t_gene_filename  = $self->get_t_gene_file(
        $t_id_input_method,
        $tempdir
    );

    return $self->error("Could not create target gene IDs file") 
		if !$t_gene_filename;

	my $bg_gene_filename;
    if ($bg_id_input_method eq 'upload' || $bg_id_input_method eq 'paste') {
		$bg_gene_filename = $self->get_bg_gene_file(
			$bg_id_input_method,
			$tempdir
		);

	    return $self->error('Could not create background gene IDs file')
            if !$bg_gene_filename;
	}


    # Close the I/O handles
    #close(STDERR);

	#
	# Call the analysis script
	#
    my $command = OPOSSUM_SCRIPTS_PATH . "/opossum_gene_tca.pl"
        . " -j $job_id"
        . " -tgf $t_gene_filename"
		. " -s $species"
        . " -d $tempdir"
        . " -web"
        . " -m $email";
    
	if ($t_gene_id_type) {
		$command .= " -tgidt $t_gene_id_type";
	}
	
    if ($bg_id_input_method eq 'upload' || $bg_id_input_method eq 'paste') {
        $command .= " -bgf $bg_gene_filename";
    } elsif ($bg_id_input_method eq 'random') {
        $command .= " -bnr $bg_num_rand_genes";
    }
	
	if ($bg_gene_id_type) {
		$command .= " -bgidt $bg_gene_id_type";
	}
	
	if ($state->has_operon()) {
		$command .= " -has_operon";
	}
	
    if ($tf_db) {
        $command .= " -tfdb $tf_db";
    }
    
    if ($cl_db) {
        $command .= " -cldb $cl_db";
    }

    #if ($collections_str) {
    #    $command .= " -co $collections_str";
    #}

    #if ($tax_groups_str) {
    #    $command .= " -tax '$tax_groups_str'";
    #}

    if ($tf_family_filename) {
        $command .= " -famf $tf_family_filename";
    }

    if ($min_ic) {
        $command .= " -ic $min_ic";
    }

    if ($conservation_level) {
        $command .= " -cl $conservation_level";
    }

    #
    # User defined threshold overrides threshold_level
    #
    if (defined $threshold) {
        $command .= " -th $threshold";
    } else {
        $command .= " -thl $threshold_level";
    }
	
    #
    # User defined search region overrides search_region_level
    #
    if (defined $upstream_bp) {
        $command .= " -up $upstream_bp";
    }

    if (defined $downstream_bp) {
        $command .= " -dn $downstream_bp";
    }

    if (!defined $upstream_bp && !defined $downstream_bp) {
        $command .= " -srl $search_region_level";
    }

    if ($result_type eq 'top_x_results') {
        $command .= " -n $num_display_results";
    } elsif ($result_type eq 'significant_hits') {
        $command .= " -zcutoff $zscore_cutoff";
        $command .= " -fcutoff $fisher_cutoff";
    }

    if ($result_sort_by) {
        $command .= " -sr $result_sort_by";
    }

    if ($user_tf_family_file) {
        $command .= " -ufamf '$user_tf_family_file'";
    }

    if ($user_t_gene_file) {
        $command .= " -utgf '$user_t_gene_file'";
    }

    if ($user_bg_gene_file) {
        $command .= " -ubgf '$user_bg_gene_file'";
    }
    
    my $submitted_time = scalar localtime(time);
    printf LOG
        "\nStarting gene-based TCA analysis at $submitted_time:\n$command\n\n";

    my $rv = system("$command >/dev/null 2>&1 &");
    if ($rv) {
        printf LOG "Analysis script returned $rv - $!\n";
    }

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        rel_htdocs_results_path => REL_HTDOCS_RESULTS_PATH,
        rel_htdocs_data_path => REL_HTDOCS_DATA_PATH,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        bg_color_class      => $state->bg_color_class(),
        heading             => $state->heading(),
        title               => $state->title(),
        section             => 'Analysis Submitted',
        job_id              => $job_id,
        submitted_time      => $submitted_time,
		species				=> $species,
        user_tf_family_file => $user_tf_family_file,
        user_t_gene_file    => $user_t_gene_file,
        user_bg_gene_file   => $user_bg_gene_file,
        tf_db               => $tf_db,
        cl_db               => $cl_db,
        #tf_cluster_set      => $tf_cluster_set,
        #collections         => \@collections,
        #tax_groups          => \@tax_groups,
        min_ic              => $min_ic,
		min_conservation	=> $disp_min_conservation,
		upstream_bp			=> $disp_upstream_bp,
		downstream_bp		=> $disp_downstream_bp,
        threshold           => $disp_threshold,
        result_type         => $result_type,
        num_display_results => $num_display_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $result_sort_by,
        email               => $email,
        var_template        => "analysis_summary_gene_tca.html"
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

    my $heading = sprintf "%s TFBS Cluster Analysis",
        ucfirst $state->species();
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    #$state->errors(undef);
    #$state->warnings(undef);
}

1;
