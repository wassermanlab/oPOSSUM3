package oPossumSeqTCAWeb;

use base 'CGI::Application';

use oPossumWebInclude;
use oPossumSeqWebInclude;
use oPossumTCAWebInclude;

#use lib OPOSSUM_LIB_PATH;   # defined in oPossumWebInclude

use Data::Dumper;    # for debugging only
use Template;
use CGI;
use File::Temp qw/ tempfile tempdir /;

use OPOSSUM::Web::State;
use TFBS::DB::JASPAR5;

use CGI::Carp qw(carpout);    # fatalsToBrowser;
#use CGI::Debug( report => 'everything', on => 'anything' );

use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_seq_tca';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_seq_tca";
$log_file .= '_devel' if DEVEL_VERSION;
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

    my $state = $self->state();

    my $q = $self->query;
    
    #my $collections = JASPAR_COLLECTIONS;
    my $dflt_collections = ['CORE','PBM','PENDING'];

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
    
=head3
    my %core_tf_sets;
    my %pending_tf_sets;
    my %pbm_tf_sets;
    foreach my $tax_group (@$tax_groups) {
        my $core_tf_set = $self->fetch_tf_set(
            -collection => 'CORE',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $core_tf_sets{$tax_group} = $core_tf_set;
        
        my $pbm_tf_set = $self->fetch_tf_set(
            -collection => 'PBM',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $pbm_tf_sets{$tax_group} = $pbm_tf_set;
        
        my $pending_tf_set = $self->fetch_tf_set(
            -collection => 'PENDING',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $pending_tf_sets{$tax_group} = $pending_tf_set;
    }
=cut

    my $sample_seqs = qq{>seq1
AAGAAAATAAAGGTAACTCAATTTAGCTCCCTGGGATCCTGACTTACATACAACTGCATA
TTTGTGATGACTGTTTTCCATCCAATAATTATATACTTTTTTTCTAGGAAACAGAAGGCA
ACTTATTGCGCTTTTCTTTTTTCCTTGTAAGTGTTTTCAGGCAGAGATTTTTTTTTTGGA
AGGTATTTCCCTTCCCCTGTG
>seq2
CATATAGTCAATACAATGTGTTTAATAAAATCAATTTATTTCCTTTTAGTCAACTAATAA
TCATGCCCTCACACAAATTATTTGCAACATTTGATTGAAAATATTCATATTCCAAGCATC
TCTGTTCTCTGACATATGTTTAAAGTGGGTTCTAAGCTATAAAACGAAACTGCTAACATA
TTAGCTTATTACAATCTGACA
>seq3
GTAGCCAAGTTGAGTGCTGGGCCAGCACTGAGTCATATGCTGATTGCATCAAAAGTAACT
TCACAAATGCTTTTCTTTCCCTGTGTGTGCTCAGTCATGCTTACATTACTGCTTTCTGAT
AATCCTGGCTAAAAAGATATAAGAGACACTCATGGAAATTAAAAGGTGATTTACAGTGCA
CCAAAGATGATGAATAGAATA
>seq4
GAAATCATATCTCTGTGCCTATTTCCTCATCTCTAGAACAAATGAATGCGTTGCCAGGCA
GCACTTAATGAGAGTACAACTGCCTTGTAAGGAAAGCATGACCAGGCAGTTTCACTGTTG
GGTGAACCTCACAGAGGCATTGAGGCATTGGCAGATACAGCCTCCTAATACTGCCGAAAG
ATACAATAAGAATATATATAT
>seq5
CATTTGGGACAAAGACTGCAGAGGAAGAAGATAATGGAGATTGGAGCCTGACTCAGCACA
GACTCTTTGGAGCATCCCTGTGAGACATTGATGCTGACTCATTGCTCGGGGCTAGTTTCC
AAGCATATAAACAGAAACAGTAAGGTCGGCTGAAAATCCACAGAGCTTTAGTGTACAAAG
CACTAAAGAAAAACAAGGGCA};

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
        #collections             => JASPAR_COLLECTIONS,
        collections             => $dflt_collections,
        bg_seq_set_keys         => BG_SEQ_SET_KEYS,
        bg_seq_set_names        => BG_SEQ_SET_NAMES,
        sid                     => $state->sid(),
        tax_groups              => $tax_groups,
        #num_tax_groups          => $num_tax_groups,
        min_ic                  => $min_ic,
        tf_cluster_set          => $tf_cluster_set,
        #core_tf_sets            => \%core_tf_sets,
        #pbm_tf_sets             => \%pbm_tf_sets,
        #pending_tf_sets         => \%pending_tf_sets,
        sample_sequences        => $sample_seqs,
        var_template            => "input_seq_tca.html"
    };

    my $output = $self->process_template('master.html', $vars);

    return $output;
}


sub process
{
    my $self = shift;

    my $state = $self->state();

    my $q = $self->query();

    #my $collections = JASPAR_COLLECTIONS;
    my $dflt_collections = ['CORE','PBM','PENDING'];

    my $dflt_min_ic     = DFLT_CORE_MIN_IC;
    #my $dflt_tax_groups = JASPAR_CORE_TAX_GROUPS;
    my $dflt_tax_groups = ['vertebrates','insects','nematodes'];
    my $num_tax_groups  = scalar @$dflt_tax_groups;

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
    my @tf_cluster_families;
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
    
    my $tf_cluster_select_method = $q->param("tf_cluster_select_method");
    if (!$tf_cluster_select_method) {
        return $self->error("No TFBS cluster selected");
    }
    
    #my $tf_cluster_select_criteria;
    #if ($tf_cluster_select_method =~ /specific/) {
    #    $tf_cluster_select_criteria = 'specific';
    #} else {
    #    $tf_cluster_select_criteria = 'all';
    #}
    
    if ($tf_cluster_select_method =~ /specific/) {
        push @tf_cluster_families, $q->param('tf_cluster_families');
        if (!@tf_cluster_families) {
            return $self->error("No specific TFBS cluster families selected");
        }
    }
    $tf_families_str = join(',', @tf_cluster_families);
    
    #my %matrix_args;
    #if (@tf_cluster_families) {
    #    $matrix_args{-families} = \@tf_cluster_families;
    #}

    #my $tf_cluster_set = $self->fetch_tf_cluster_set(%matrix_args);

    #if (!$tf_cluster_set) {
    #    return $self->error("Error fetching TFBSCluster::TFClusterSet");
    #}
    
    #@tf_cluster_ids = @{$tf_cluster_set->ids()};
    #$tf_cluster_ids_str = join ",", @tf_cluster_ids;
    
    my $seq_input_method = $q->param('seq_input_method');
    if (!$seq_input_method) {
        $self->error("No target sequences specified");
    }

    my $bg_seq_input_method = $q->param('bg_seq_input_method');
    if (!$bg_seq_input_method) {
        $self->error("No background sequence specified");
    }

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

    my $seq_filename  = $self->get_seq_file(
        $seq_input_method,
        $tempdir
    );

    return $self->error() if !$seq_filename;

    my $bg_seq_filename = $self->get_back_seq_file(
        $bg_seq_input_method,
        $tempdir
    );

    return $self->error('No background sequence file') if !$bg_seq_filename;

    # Close the I/O handles
    #close(STDERR);

    my $command = "./opossum_seq_tca.pl"
        . " -j $job_id"
        . " -s $seq_filename"
        . " -b $bg_seq_filename"
        . " -d $tempdir"
        . " -th '$threshold%'"
        . " -m $email";
    
    if ($tf_db) {
        $command .= " -tdb $tf_db";
    }
    
    if ($cl_db) {
        $command .= " -cdb $cl_db";
    }

    if ($collections_str) {
        $command .= " -co $collections_str";
    }

    if ($tax_groups_str) {
        $command .= " -tax '$tax_groups_str'";
    }

    if ($min_ic) {
        $command .= " -ic $min_ic";
    }

    #if ($tf_cluster_ids_str) {
    #    $command .= " -ids '$tf_cluster_ids_str'";
    #}
    
    if ($tf_families_str) {
        $command .= " -fam '$tf_families_str'";
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

    if ($user_seq_file) {
        $command .= " -usf $user_seq_file";
    }

    if ($user_bg_seq_file) {
        $command .= " -ubsf $user_bg_seq_file";
    }

    if ($bg_seq_set_name) {
        $command .= " -bss $bg_seq_set_name";
    }

    printf LOG "\nStarting sequence-based analysis at %s:\n$command\n\n",
        scalar localtime(time);

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
        user_seq_file       => $user_seq_file,
        user_bg_seq_file    => $user_bg_seq_file,
        bg_seq_set_name     => $bg_seq_set_name,
        #tf_select_criteria  => $tf_cluster_select_criteria,
        tf_db               => $tf_db,
        cl_db               => $cl_db,
        tf_cluster_set      => $tf_cluster_set,
        collections         => \@collections,
        tax_groups          => \@tax_groups,
        min_ic              => $min_ic,
        threshold           => $threshold,
        result_type         => $result_type,
        num_display_results => $num_display_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $result_sort_by,
        email               => $email,
        var_template        => "analysis_summary_seq_tca.html"
    };

    #print STDERR "results vars:\n" . Data::Dumper::Dumper($vars);

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

sub db
{
    return$_[0]->{db};
}

sub jaspar_core_tf_list
{
    return $_[0]->{jaspar_core_tf_list};
}

sub _get_core_tf_list
{
    my $self = shift;

    my $db = $self->db;
    if (!$db) {
         die "no database connection\n";
    }
 
    my $core_matrixset = $db->get_MatrixSet(-collection => "CORE");
    if (!$core_matrixset) {
        die "error getting CORE matrixset\n";
    }

    my $core_matrix_iterator = $core_matrixset->Iterator();
    my @jaspar_core_tfs;
    while (my $core_matrix = $core_matrix_iterator->next) {
        push @jaspar_core_tfs, $core_matrix;
    }

    @jaspar_core_tfs = sort {uc $a->name cmp uc $b->name} @jaspar_core_tfs;

    $self->{jaspar_core_tf_list} = \@jaspar_core_tfs;
}

sub initialize_state
{
    my ($self, $state) = @_;

    my $heading = "Sequence-based TFBS Cluster Analysis";
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    $state->errors(undef);
    $state->warnings(undef);
}

1;
