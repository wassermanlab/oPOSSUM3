package SeqACTCA;

use base 'CGI::Application';

use OPOSSUM::Opt::ACSAOpt;
#use oPossumOpt::oPossumTCAOpt;

use OPOSSUM::Web::Include::BaseInclude;
use OPOSSUM::Web::Include::SeqInclude;
use OPOSSUM::Web::Include::TCAInclude;

use Data::Dumper;    # for debugging only
use Template;
use CGI::Carp qw(carpout);    # fatalsToBrowser;
use File::Temp qw/ tempdir /;

use OPOSSUM::Web::State;


use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_seq_actca';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_seq_actca";
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

    $self->_clean_tempfiles;
    $self->_clean_resultfiles;
}

sub input
{
    my $self = shift;


    my $q = $self->query;

    my $state = $self->state();

    #my $collections = JASPAR_COLLECTIONS;
    #my $dflt_collections = ['CORE','PBM','PENDING'];
    my $dflt_collections = JASPAR_COLLECTIONS_USED;

    #my $dflt_tax_groups = JASPAR_CORE_TAX_GROUPS;
    my $dflt_tax_groups = ['vertebrates','insects','nematodes'];
    
    my $min_ic     = DFLT_CORE_MIN_IC;
    #my $tax_groups = JASPAR_CORE_TAX_GROUPS;
    my $tax_groups = $dflt_tax_groups;

    my $num_tax_groups = scalar @$tax_groups;

    #
    # Connect to oPOSSUM_cluster DB and retrieve TFCluster info
    # Connect to JASPAR DB and retrieve TF info
    #
    $self->jaspar_db_connect();
    $self->opossum_cluster_db_connect();
    
    my $tf_cluster_set = $self->fetch_tf_cluster_set();

    # Have the user select an anchor TF
    # then search using the TFCluster that the anchor TF belongs to
    my $tf_set = $self->fetch_tf_set(
        -collection => 'CORE',
        -tax_group  => $tax_groups,
        -min_ic     => $min_ic
    );

    #printf STDERR "tf_set:\n%s\n", Data::Dumper::Dumper($tf_set);

    #
    # Format a tax group list for display on the web page. The $tax_group
    # variable is obtained from the db and may be a single value or a comma
    # separated list (no space after comma).
    #
    my @s_tax_groups;
    foreach my $tg (@$tax_groups) {
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
        dflt_threshold          => DFLT_THRESHOLD,
        min_threshold           => MIN_THRESHOLD,
        dflt_inter_binding_dist => DFLT_INTER_BINDING_DIST,
        max_inter_binding_dist  => MAX_INTER_BINDING_DIST,
        collections             => $dflt_collections,
        bg_seq_set_keys         => BG_SEQ_SET_KEYS,
        bg_seq_set_names        => BG_SEQ_SET_NAMES,
        sid                     => $state->sid(),
        tax_groups              => $tax_groups,
        tax_group_list          => $tax_group_list,
        min_ic                  => $min_ic,
        tf_cluster_set          => $tf_cluster_set,
        tf_set                  => $tf_set,
        var_template            => "input_seq_actca.html"
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

    #my $collections = JASPAR_COLLECTIONS;
    #my $dflt_collections = ['CORE','PBM','PENDING'];
    my $dflt_collections = JASPAR_COLLECTIONS_USED;

    my $dflt_min_ic     = DFLT_CORE_MIN_IC;
    #my $dflt_tax_groups = JASPAR_CORE_TAX_GROUPS;
    my $dflt_tax_groups = ['vertebrates','insects','nematodes'];
    my $num_tax_groups  = scalar @$dflt_tax_groups;

    my $inter_binding_dist = $q->param('inter_binding_dist');
    if ($inter_binding_dist > MAX_INTER_BINDING_DIST) {
        return $self->error(
            sprintf(
                "Specified inter-binding distance exceeds maximum of %d allowed",
                MAX_INTER_BINDING_DIST
            )
        );
    }

    my $threshold           = $q->param('threshold');
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

    unless ($threshold) {
        return $self->error(
            "No TFBS profile matrix score threshold provided."
        );
    }
    
	#printf LOG "Email and threshold set\n";

    my $tf_db = JASPAR_DB_NAME;
    my $cl_db = TFBS_CLUSTER_DB_NAME;


    #my @tf_cluster_ids;
    #my $tf_cluster_ids_str;
    my @collections;
    my $collections_str;
    my $min_ic;
    my @tax_groups;
    my $tax_groups_str;
    #my $matrix_set;
    my $tf_cluster_set;
    my @tf_families;
    my $tf_families_str;
    
    $min_ic = $q->param('tf_min_ic');
    
    push @tax_groups, $q->param('tf_tax_groups');
    if (!@tax_groups or scalar(@tax_groups == 0)) {
        push @tax_groups, @$dflt_tax_groups;
    }
    $tax_groups_str = join(',', @tax_groups);
    
    
    push @collections, $q->param('tf_collections');
    if (!@collections or scalar(@collections == 0)) {
        #push @collections, @$dflt_collections;
        return $self->error("No JASPAR collection selected");
    }
    $collections_str = join(',', @collections);
    
    my $tf_family_select_method = $q->param("tf_family_select_method");
    if (!$tf_family_select_method) {
        return $self->error("No TFBS cluster selected");
    }
    
    #my $tf_family_select_criteria;
    #if ($tf_family_select_method =~ /specific/) {
    #    $tf_family_select_criteria = 'specific';
    #} else {
    #    $tf_family_select_criteria = 'all';
    #}
    
    my $user_tf_family_file;
    if ($tf_family_select_method =~ /specific/) {
        push @tf_families, $q->param('tf_families');
        if (!@tf_families or scalar (@tf_families == 0)) {
            return $self->error("No specific TFBS cluster families selected");
        }
        $tf_families_str = join(',', @tf_families);
    } elsif ($tf_family_select_method eq 'upload') {
        $user_tf_family_file = $q->param('tf_family_upload_file');
        $user_tf_family_file =~ s/.*\///;
    }

    my $anchor_tf_id = $q->param('anchor_tf_id');
    if (!$anchor_tf_id) {
        return $self->error("Anchoring TF not specified");
    }
    $state->anchor_tf_id($anchor_tf_id);

    my $seq_input_method = $q->param('seq_input_method');
    if (!$seq_input_method) {
        return $self->error("No target sequences specified");
    }

    my $bg_seq_input_method = $q->param('bg_seq_input_method');
    if (!$bg_seq_input_method) {
        return $self->error("No background sequence specified");
    }

	#printf LOG "Input parameters stored\n";

    my $user_seq_file;
    if ($seq_input_method eq 'upload') {
        $user_seq_file = $q->param('seq_file');

        $user_seq_file =~ s/.*\///;
    }

    my $user_bg_seq_file;
    my $bg_seq_set_name;
    if ($bg_seq_input_method eq 'upload') {
        $user_bg_seq_file = $q->param('bg_seq_file');

        $user_bg_seq_file =~ s/.*\///;
    } elsif ($bg_seq_input_method eq 'default') {
        my $bg_seq_set_names = BG_SEQ_SET_NAMES;
        my $bg_seq_set_key = $q->param('bg_seq_set_key');

        $bg_seq_set_name = $bg_seq_set_names->{$bg_seq_set_key};
    }

	#printf LOG "User sequence files stored\n";

    #
    # Create a temporary working directory for all the input sequence and
    # output results file as a sub-directory under the defined temp. dir.
    #
    my $tempdir = tempdir(DIR => ABS_HTDOCS_RESULTS_PATH);
    
    # for testing
    #my $job_id = 'd9yqr8V7WJ';
    #my $tempdir = $job_id;
    my $job_id = $tempdir;
    $job_id =~ s/.*\///;

    my $tf_family_filename;
    if ($tf_family_select_method eq 'upload') {
        $tf_family_filename = $self->get_tf_family_file(
            $tf_family_select_method,
            $tempdir
        );

        return $self->error() if !$tf_family_filename;
    }

    my $seq_filename  = $self->get_seq_file(
        $seq_input_method,
        $tempdir
    );

    unless ($seq_filename) {
        if ($seq_input_method eq 'upload') {
            return $self->error(
                "There was a problem with the provided target sequence file"
            );
        } elsif ($seq_input_method eq 'paste') {
            return $self->error(
                "There was a problem with the pasted target sequence text"
            );
        }
    }

    my $bg_seq_filename = $self->get_back_seq_file(
        $bg_seq_input_method,
        $tempdir
    );

    unless ($bg_seq_filename) {
        if ($bg_seq_input_method eq 'upload') {
            return $self->error(
                "There was a problem with the provided background sequence file"
            );
        } elsif ($bg_seq_input_method eq 'paste') {
            return $self->error(
                "There was a problem with the pasted background sequence text"
            );
        } else {
            return $self->error();
        }
    }

	#printf LOG "\n\ntempdir = $tempdir\n";


	#
	# Call the analysis script
	#
    my $command = OPOSSUM_SCRIPTS_PATH . "/opossum_seq_actca.pl"
   		. " -j $job_id"
        . " -tsf $seq_filename"
	    . " -bsf $bg_seq_filename"
	    . " -d $tempdir"
		. " -atfid $anchor_tf_id"
        . " -dist $inter_binding_dist"
		. " -web"
		. " -m $email";

		#printf LOG "1:  $command\n\n";

    if ($tf_db) {
        $command .= " -tfdb $tf_db";
    }
    
    if ($cl_db) {
        $command .= " -cldb $cl_db";
    }

	#printf LOG "2:  $command\n\n";
    if ($collections_str) {
        $command .= " -co '$collections_str'";
    }

    if ($tax_groups_str) {
        $command .= " -tax '$tax_groups_str'";
    }

    if ($tf_families_str) {
        $command .= " -fam '$tf_families_str'";
    } elsif ($tf_family_filename) {
        $command .= " -famf $tf_family_filename";
    }

	#printf LOG "3:  $command\n\n";
    if ($min_ic) {
        $command .= " -ic $min_ic";
    }

    if (defined $threshold) {
        $command .= " -th '$threshold%'";
    }

	#printf LOG "4:  $command\n\n";
    if ($result_type eq 'top_x_results') {
        $command .= " -n $num_display_results";
    } elsif ($result_type eq 'significant_hits') {
        $command .= " -zcutoff $zscore_cutoff";
        $command .= " -fcutoff $fisher_cutoff";
    }

	#printf LOG "5:  $command\n\n";
    if ($result_sort_by) {
        $command .= " -sr $result_sort_by";
    }

	#printf LOG "6:  $command\n\n";
    if ($user_seq_file) {
        $command .= " -utsf $user_seq_file";
    }

    if ($user_bg_seq_file) {
        $command .= " -ubsf $user_bg_seq_file";
    }

    if ($bg_seq_set_name) {
        $command .= " -bss '$bg_seq_set_name'";
    }

    if ($user_tf_family_file) {
        $command .= " -ufamf $user_tf_family_file";
    }

	#printf LOG "Done\n";
	my $submitted_time = scalar localtime(time);
	#printf LOG "Prepared script command at $submitted_time";
	#printf LOG "\n";
	#printf LOG "Final:  $command\n\n";

	printf LOG "\nStarting aCTCA analysis at $submitted_time:\n$command\n\n";

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
        bg_color_class      => BG_COLOR_CLASS,
        heading             => $state->heading(),
        title               => $state->title(),
        section             => 'Analysis Submitted',
        job_id              => $job_id,
        submitted_time      => $submitted_time,
		user_seq_file		=> $user_seq_file,
		user_bg_seq_file	=> $user_bg_seq_file,
        user_tf_family_file => $user_tf_family_file,
        tf_db               => $tf_db,
        cl_db               => $cl_db,
        tf_cluster_set      => $tf_cluster_set,
        anchor_tf_id        => $anchor_tf_id,
        collections         => \@collections,
        tax_groups          => \@tax_groups,
		tf_families			=> \@tf_families,
        min_ic              => $min_ic,
        threshold           => $threshold,
        result_type         => $result_type,
        num_display_results => $num_display_results,
        result_sort_by      => $result_sort_by,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        email               => $email,
        var_template        => "analysis_summary_seq_actca.html"
    };

    #print STDERR "results vars:\n" . Data::Dumper::Dumper($vars);

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

sub initialize_state
{
    my ($self, $state) = @_;

    my $heading = "Sequence-based Anchored Combination TFBS Cluster Analysis";
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    #$state->errors(undef);
    #$state->warnings(undef);
}

1;
