#!/usr/bin/env perl

=head1 NAME

opossum_seq_tca.pl

=head1 SYNOPSIS
  opossum_seq_ssa.pl
      -d results_dir
      -s t_seq_file
      -b bg_seq_file
      [-spf t_peak_pos_file]
      [-bpf bg_peak_pos_file]
      [-db tf_database]
      [-cdb tf_cluster_database]
      [-famf tf_family_file] | [-fam tf_families]
      [-co collections] [-tax tax_groups] [-ic min_ic]
      [-th threshold]
      [-n num_results | -zcutoff cutoff -fcutoff cutoff -kscutoff cutoff]
      [-sr sort_by]
      [-nh]
      [-web]
      [-j job_id]
      [-m email]
      [-usf user_seq_file]
      [-ubsf user_bg_seq_file]
      [-uspf user_peak_pos_file]
      [-ubpf user_bg_peak_pos_file]
      [-ufamf user_tf_family_file]
      [-bss bg_seq_set_name]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

Some switches can take multiple values. To specify multiple values either
use multiple instances of the switch or a single switch followed by a comma
separated string of values (or some combination thereof).
e.g: -tax vertebrates -tax "insects, nematodes"

    -d, -dir dir
            Name of directory used for output results files. If the
            directory does not already exist it will be created.

    -s, -tsf t_seq_file
            Input target sequences file (fasta format).

    -b, -bsf bg_seq_file
            Input background sequences file (fasta format).

    -spf, -tpf t_peak_pos_file
            Max confidence (peak max) position file for input target
            sequences.

    -bpf, -bpf bg_peak_pos_file
            Max confidence (peak max) position file for input background
            sequences.

    -db, -tfdb db_name
            Specifies which TF database to use; default = JASPAR_2010.

    -fam tf_families
            Specify one or more JASPAR TF families to use.

    -famf, -tff tf_family_file
            File containing the list of JASPAR TF families to search for.

    -co collections
            Specify one or more JASPAR TFBS profile collections;
            default = CORE.

    -tax tax_groups
            Limit the analysis to use only JASPAR TFBS profile matrices
            which belong to one or more of these tax groups, e.g:
              -tax vertebrates -tax insects -tax nematodes
              -tax "vertebrates,insects,nematodes"

    -ic min_ic
            Specify minimum information content (specificity) of JASPAR
            TFBS profile matrices to use.

    -th threshold
            Minimum relative TFBS position weight matrix (PWM) score to
            report in the analysis. The thresold may be spesified as a
            percentage string, e.g. '85%', or as a decimal number, e.g.
            0.85
            Default = '80%' (min. = '75%')
    -n num_results
            The number of results to output. Numeric or string 'All'.
            Default = 'All'

    -zcutoff score
            Z-score cutoff of results to display. Only output results with
            at least this Z-score.

    -fcutoff score
            Fisher score cutoff of results to display. Only output results
            with at least this Fisher score.

    -kscutoff score
            KS score cutoff of results to display. Only output results
            with at least this KS score.

    -nh
            Do NOT write TF binding site details. Unless specified, in
            addition to the main results ranking of TFBS
            over-representation, files detailing the target sequence binding
            site positions are written for each of the input TFs (one file
            per TF).

    -sr sort_by
            Sort results by this score ('zscore', 'fisher'), highest to
            lowest.
            Default = 'zscore'.

    -help, -h, -?
            Help. Print usage message and exit.

=head2 Web Server Specific Options

    The following options are passed to the script by web-based oPOSSUM.
    These are not required when running the scripts directly on the command
    line and can generally be ignored.

    -web
            Web server switch. Indicates that the script caller is the web
            server, and HTML results files should also be created.

    -j, -job_id job_ID
            The oPOSSUM job ID.

    -m email
            E-mail address of user. An e-mail is sent to the user to notify
            him/her when the analysis has completed with a URL to the
            HTML results page.

    -usf, -utsf user_seq_file
            User supplied target sequences upload file name.

    -ubsf user_bg_seq_file
            User supplied background sequences upload file name.

    -uspf, -utpf user_peak_pos_file
            User supplied max confidence (peak max) position file name for
            input target sequences.

    -ubpf user_bg_peak_pos_file
            User supplied max confidence (peak max) position file name for
            input background sequences.

    -ufamf, -utff user_tf_family_file
            User supplied TF family upload file name.

    -bss bg_seq_set_name
            Background sequence set name.

=head1 DESCRIPTION

Take a set of peak sequences as provided in the input peaks file and 
a set of control peak sequences. Optionally take: 1) maximum confidence
position files, 2) a subset of transcription factors (TFs) limited by external 
(JASPAR) database name and information content or taxonomic supergroup or all 
TFs in the oPOSSUM database, 3) TFBS structural families to include in analysis, and 4) PWM score threshold.

Count the number of TFBS clusters for each cluster type which was found at the 
given PWM score threshold for both the test and control set of peaks. Perform
Fisher exact test, z-score and KS analyses, and output these results to the
output file. Optionally write details of TFBS clusters found in test set to 
detailed TFBS cluster hits file.

=head1 AUTHOR

  David Arenillas, Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use warnings;

use lib '/apps/oPOSSUM3/lib';

use Getopt::Long;
use Pod::Usage;
use Carp;

use File::Temp;
use File::Temp qw/ tempfile tempdir /;
use Log::Log4perl qw(get_logger :levels);
#use Data::Dumper;

use Bio::SeqIO;

#use TFBS::DB::JASPAR5;

#use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::Include::SeqTCAInclude;
use OPOSSUM::TFSet;
use OPOSSUM::Analysis::Cluster::Zscore;
use OPOSSUM::Analysis::Cluster::Fisher;
use OPOSSUM::Analysis::Cluster::KS;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::Values;
use OPOSSUM::Analysis::Cluster::CombinedResultSet;

use Statistics::Distributions;

#use constant DEBUG          => 0;

my $help;
my $job_id;
my $t_seq_file;
my $bg_seq_file;
my $t_peak_pos_file;
my $bg_peak_pos_file;
my $results_dir;
my $web;
my $tf_db;
my $tfcl_db;
my @collections;
my @tax_groups;
my @tf_families;
my $tf_family_file;
my $min_ic;
my $threshold;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $ks_cutoff;
my $sort_by;
my $nh;
my $email;
my $user_tf_family_file;
my $user_seq_file;
my $user_bg_seq_file;
my $user_peak_pos_file;
my $user_bg_peak_pos_file;
my $bg_seq_set_name;
GetOptions(
    'dir|d=s'       => \$results_dir,
    'tsf|s=s'       => \$t_seq_file,
    'bsf|b=s'       => \$bg_seq_file,
    'tpf|spf=s'     => \$t_peak_pos_file,
    'bpf=s'         => \$bg_peak_pos_file,
    'tfdb|tdb|db=s' => \$tf_db,
    'cldb|cdb=s'    => \$tfcl_db,
    'co=s'          => \@collections,
    'tax=s'         => \@tax_groups,
    'ic=s'          => \$min_ic,
    'fam=s'         => \@tf_families,
    'famf|tff=s'    => \$tf_family_file,
    'th=s'          => \$threshold,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'kscutoff=f'    => \$ks_cutoff,
    'sr=s'          => \$sort_by,
    'nh'            => \$nh,
    'm=s'           => \$email,
    'web'           => \$web,
    'job_id|j=s'    => \$job_id,
    'utsf|usf=s'    => \$user_seq_file,
    'ubsf=s'        => \$user_bg_seq_file,
    'utpf|uspf=s'   => \$user_peak_pos_file,
    'ubpf|ubpf=s'   => \$user_bg_peak_pos_file,
    'ufamf|utff=s'  => \$user_tf_family_file,
    'bss=s'         => \$bg_seq_set_name,
    'help|h|?'      => \$help
);

if ($help) {    
    pod2usage(
        -verbose    => 1
    );
}

unless ($results_dir) {
    pod2usage(
        -msg        => "No results directory specified."
                     . " Use -help switch for detailed option list.\n",
        -verbose    => 1
    );
}

my $message = "";
my $ok = 1;

my $abs_results_dir = $results_dir;
my $rel_results_dir = $results_dir;

# called by web server? then create web pages
if ($web) {
    # Remove absolute path
    $rel_results_dir =~ s/.*\///;

    # Add relative path
    $rel_results_dir = REL_HTDOCS_RESULTS_PATH . "/$rel_results_dir";
} else {
    unless (-d $results_dir) {
        mkdir $results_dir
            || die "Error creating results directory $results_dir - $!\n";;
    }

    #$abs_results_dir = $ENV{PWD} . "/$rel_results_dir";
}

#
# Initialize logging
#
my $log_file = get_log_filename('opossum_seq_tca', $results_dir);

my $logger = get_logger();
if (DEBUG) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}

#my $layout = Log::Log4perl::Layout::PatternLayout->new("%M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("[%d] %L %p\t%m%n");

my $appender = Log::Log4perl::Appender->new(
    "Log::Dispatch::File",
    filename    => $log_file,
    mode        => "append"
);

$appender->layout($layout);
$logger->add_appender($appender);

#open(ERR, ">>$log_file") || die "Could not open log file $log_file\n";

#carpout(\*ERR);

$logger->info("Starting analysis");

my $heading = "Sequence-Based TFBS Cluster Analysis";

my %job_args = (
    -job_id => $job_id,
    -heading => $heading,
    -email => $email,
    -logger => $logger
);

#
# check for required arguments
#

unless ($t_seq_file) {
    fatal("No sequences file specified", %job_args);
}
unless (-e $t_seq_file) {
    fatal("Sequences file $t_seq_file does not exist", %job_args);
}
unless ($bg_seq_file) {
    fatal("No background sequences file specified", %job_args);
}
unless (-e $bg_seq_file) {
    fatal("Background sequences file $bg_seq_file does not exist", %job_args);
}

$logger->info("Abs results dir: $abs_results_dir");
$logger->info("Rel results dir: $rel_results_dir");

#
# read in target sequences
#
$logger->info("Reading target sequences from $t_seq_file");
my $t_seqs = read_seqs($t_seq_file)
    || fatal("Could not read sequences", %job_args);

check_seqs($t_seqs, "target", %job_args);

print "t_seqs: " . scalar @$t_seqs . "\n";

my ($t_seq_ids, $t_seq_id_seqs, $t_seq_id_display_ids) = id_seqs($t_seqs);

my $t_peak_pos;
if ($t_peak_pos_file) {
    $logger->info("Reading target sequence peak max pos file from $t_peak_pos_file");
    $t_peak_pos = read_peak_pos($t_peak_pos_file, $t_seq_ids)
    || fatal("Could not read peak max positions", %job_args);
}

my ($t_seq_gc_content, $t_seq_len) = calc_seqs_total_gc_content($t_seqs);

$logger->debug("Target seq GC content: $t_seq_gc_content");
$logger->debug("Target seq total length: $t_seq_len");

#
# read in background sequences
#
$logger->info("Reading background sequences from $bg_seq_file");
my $bg_seqs = read_seqs($bg_seq_file)
    || fatal("Could not read background sequences", %job_args);

check_seqs($bg_seqs, "background", %job_args);

my ($bg_seq_ids, $bg_seq_id_seqs, $bg_seq_id_display_ids) = id_seqs($bg_seqs);

my $bg_peak_pos;
if ($bg_peak_pos_file) {
    $logger->info("Reading background sequence peak max pos file from $bg_peak_pos_file");
    $bg_peak_pos = read_peak_pos($bg_peak_pos_file, $bg_seq_ids)
    || fatal("Could not read peak max positions", %job_args);
}

my ($bg_seq_gc_content, $bg_seq_len) = calc_seqs_total_gc_content($bg_seqs);

$logger->debug("Background seq GC content: $bg_seq_gc_content");
$logger->debug("Background seq total length: $bg_seq_len");


#
# JASPAR / TF parameter settings
# Also, TFBSCluster settings
#

$tf_db          = JASPAR_DB_NAME if !$tf_db;

# if this is default, JASPAR has to be default as well
$tfcl_db = TFBS_CLUSTER_DB_NAME if !$tfcl_db;

if ($tfcl_db eq TFBS_CLUSTER_DB_NAME and $tf_db ne JASPAR_DB_NAME) {
    fatal("\nDefault TFBS cluster DB requires the default JASPAR DB",
        %job_args);
}

my $collections_str;
if (@collections) {
    $collections_str = join(',', @collections);
    push @collections, split(/\s*,\s*/, $collections_str);
} else {
    $collections_str = 'CORE';
    push @collections, 'CORE';
}

my $tax_groups_str;
if (@tax_groups) {
    $tax_groups_str = join(',', @tax_groups);
    @tax_groups  = split(/\s*,\s*/, $tax_groups_str);

    unless (@tax_groups) {
        fatal("Error parsing tax groups", %job_args);
    }
}

my $tf_families_str;
if (@tf_families) {
    $tf_families_str = join(',', @tf_families);
    @tf_families  = split(/\s*,\s*/, $tf_families_str);
} elsif ($tf_family_file) {    
    $logger->info("Reading TFBS family file");

    @tf_families = @{read_tf_families_from_file($tf_family_file, %job_args)};

    unless (@tf_families) {
        fatal("No TFBS families read from file $tf_family_file", %job_args);
    }

    $tf_families_str = join(',', @tf_families);
}


# AK: No longer either/or: you can have specific families AND min_ic
# The question is: what do you call tf_select_criteria now?
# Do I actually need this parameter?
# Hmm
# even if you specify the families, there's no reason it shouldn't also
# filter on min_ic and tax groups
# something that might be worthwhie porting over to gene based tca
# This should work in JASPAR5

#
# Connect to JASPAR database and retrieve the matrices
#

# this is different from gene-based...
# need to go over in detail over different possibilities
$threshold = DFLT_THRESHOLD . "%" if !$threshold;

#
# Connect to JASPAR database and retrieve the matrices
#

my $jdb = jaspar_db_connect($tf_db)
    || fatal("Could not connect to JASPAR database $tf_db", %job_args);

my %get_matrix_args = (
    -matrixtype => 'PFM'
);

$get_matrix_args{-collection} = \@collections if @collections;
$get_matrix_args{-tax_group} = \@tax_groups if @tax_groups;
$get_matrix_args{-min_ic} = $min_ic || DFLT_CORE_MIN_IC;
$get_matrix_args{-family} = \@tf_families if @tf_families;


#$logger->debug("Fetching TF matrix arguments" .
#Data::Dumper::Dumper(%get_matrix_args));

$logger->debug("Fetch TF matrices");
my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

unless ($matrix_set && $matrix_set->size > 0) {
    fatal("Error fetching TF profiles from JASPAR DB", %job_args);
}

#$logger->debug("Matrix set:\n" . Data::Dumper::Dumper($matrix_set));

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

#
# retrieve the TFBS cluster set
#

my $cdb = tfbs_cluster_db_connect($tfcl_db);
if (!$cdb) {
    fatal("Could not connect to $tfcl_db", %job_args);
}

# Note, cluster ids only consist of numbers
my $tf_cluster_set = fetch_tf_cluster_set($cdb, %get_matrix_args);

unless ($tf_cluster_set and $tf_cluster_set->size > 0) {
    fatal("Error fetching TFBS clusters", %job_args);
}
#$logger->debug("TFBSClusterSet:\n" . Data::Dumper::Dumper($tf_cluster_set));

#my $tf_ids = $tf_set->ids();
my $tfcl_ids = $tf_cluster_set->ids();

#
# Search sequences with the matrix set
#

$logger->info("Searching target sequences for TFBS clusters");
#my ($t_tf_seq_siteset, $t_tf_seqs)
my $t_cluster_seq_siteset
    = tf_cluster_set_search_seqs($tf_cluster_set, $tf_set,
                                 $t_seq_id_seqs, $threshold);

fatal("No TFBS clusters found in test sequences", %job_args)
    if !$t_cluster_seq_siteset;

$logger->info("Searching background sequences for TFBSs");
#my ($bg_tf_seq_siteset, $bg_tf_seqs)
my $bg_cluster_seq_siteset
    = tf_cluster_set_search_seqs($tf_cluster_set, $tf_set,
                                 $bg_seq_id_seqs, $threshold);

fatal("No TFBS clusters found in control sequences", %job_args)
    if !$bg_cluster_seq_siteset;

#
# Compute TFBS cluster counts
#

$logger->info("Computing target sequence TFBS counts");
my $t_counts = compute_cluster_counts(
    $tf_cluster_set, $t_seq_ids, $t_cluster_seq_siteset);

#$logger->debug("Target TFBS counts:\n" . Data::Dumper::Dumper($t_counts));

$logger->info("Computing background sequence TFBS counts");
my $bg_counts = compute_cluster_counts(
    $tf_cluster_set, $bg_seq_ids, $bg_cluster_seq_siteset);

unless ($t_counts) {
    fatal("Error fetching target sequence TFBS cluster counts", %job_args);
}

unless ($bg_counts) {
    fatal("Error fetching background sequence TFBS cluster counts", %job_args);
}

#
# Compute peak to site distances
#

$logger->info("Computing sequence peak to TFBS distances");

my $t_dists = compute_cluster_peak_distances(
    $tf_cluster_set, $t_seq_id_seqs, $t_seq_id_display_ids, $t_cluster_seq_siteset, $t_peak_pos
);

my $bg_dists = compute_cluster_peak_distances(
    $tf_cluster_set, $bg_seq_id_seqs, $bg_seq_id_display_ids, $bg_cluster_seq_siteset, $bg_peak_pos
);

#
# Score calculations
#

my $fisher = OPOSSUM::Analysis::Cluster::Fisher->new();
fatal("Error initializing Fisher analysis", %job_args)
    unless $fisher;

$logger->info("Computing Fisher scores");
my $fresult_set = $fisher->calculate_Fisher_probability(
    $bg_counts,
    $t_counts
);
fatal("Error performing Fisher analysis", %job_args)
    unless $fresult_set;

#$logger->debug("Fisher results:\n" . Data::Dumper::Dumper($fresults));

my $zscore = OPOSSUM::Analysis::Cluster::Zscore->new();
fatal("Error initializing z-score analysis", %job_args)
    unless $zscore;

$logger->info("Computing z-scores");
my $zresult_set = $zscore->calculate_Zscore(
    $bg_counts,
    $t_counts,
    $bg_seq_len,
    $t_seq_len,
);
fatal("Error computing z-score", %job_args)
    unless $zresult_set;

#$logger->debug("Z-score results:\n" . Data::Dumper::Dumper($zresults));

################

my $ks = OPOSSUM::Analysis::Cluster::KS->new();
fatal("Error initializing KS-test", %job_args)
    unless $ks;

$logger->info("Computing KS for site to peak distances");

#my $ksresult_set = $ks->calculate_KS_with_bg_distribution(
#    DFLT_PEAK_DIST_DISTRIBUTION,
#    $t_dists
#);
my $ksresult_set = $ks->calculate_KS_with_bg_values(
    $bg_dists,
    $t_dists
);
fatal("Error computing KS-test", %job_args)
    unless $ksresult_set;

#$logger->debug("KS test results:\n" . Data::Dumper::Dumper($ksresults));


#
# Use new OPOSSUM::Analysis::CombinedResultSet to combine Fisher and
# Z-score result sets.
#
$logger->info("Combining Fisher and Z-scores");
my $cresult_set = OPOSSUM::Analysis::Cluster::CombinedResultSet->new(
    -fisher_result_set  => $fresult_set,
    -zscore_result_set  => $zresult_set,
    -ks_result_set      => $ksresult_set
);
fatal("Error combining Fisher and z-score result_set", %job_args)
    unless $cresult_set;

#
# Get results as a list
#
my %result_params;
$result_params{-num_results} = $num_results if defined $num_results;
$result_params{-zscore_cutoff} = $zscore_cutoff if defined $zscore_cutoff;
$result_params{-fisher_cutoff} = $fisher_cutoff if defined $fisher_cutoff;
$result_params{-ks_cutoff} = $ks_cutoff if defined $ks_cutoff;

if (defined $sort_by) {
    if ($sort_by =~ /^fisher/) {
        $sort_by = 'fisher_p_value';
    } elsif ($sort_by =~ /^z_score/ || $sort_by =~ /^z-score/) {
        $sort_by = 'zscore';
    } elsif ($sort_by =~ /^ks/) {
        $sort_by = 'ks_p_value';
    }

    $result_params{-sort_by} = $sort_by;

    # should check if Fisher is in p-val form here or not
    #if ($sort_by eq 'zscore') {
        # Sort z-score from highest to lowest
        $result_params{-reverse} = 1;
    #}
}

$logger->info("Getting filtered/sorted result list");
my $cresults = $cresult_set->get_list(%result_params);

unless ($cresults) {
    $message = "No TFBSs scored above the selected Z-score/Fisher/KS"
        . " thresholds";
    $logger->info($message);
    $ok = 0;
    #
    # XXX
    # The ok/message stuff is not handled properly later, so just fatal it
    # for now.
    # DJA 2012/05/03
    #
    fatal($message, %job_args);
}


if ($web) {
    $logger->info("Writing HTML results");
    write_results_html($tf_cluster_set, $cresults); 
}

if ($ok)
{
    $logger->info("Writing text results");
    my $out_file = "$abs_results_dir/" . RESULTS_TEXT_FILENAME;
    write_results_text($out_file, $cresults, $tf_cluster_set, %job_args);
    
    if (!$nh) {
        $logger->info("Writing TFBS details");
        write_tfbs_cluster_details();
    }

}

$logger->info("Sending notification email to $email");

$job_args{-web} = $web;
$job_args{-abs_results_dir} = $abs_results_dir;
$job_args{-rel_results_dir} = $rel_results_dir;
$job_args{-user_t_file} = $user_seq_file;
$job_args{-user_bg_file} = $user_bg_seq_file;
$job_args{-user_t_peak_file} = $user_peak_pos_file;
$job_args{-user_bg_peak_file} = $user_bg_peak_pos_file;
$job_args{-t_num} = scalar @$t_seqs;
$job_args{-bg_num} = scalar @$bg_seqs;
$job_args{-tf_db} = $tf_db;
$job_args{-collections} = $collections_str;
$job_args{-tax_groups} = $tax_groups_str;
$job_args{-families} = $tf_families_str;
$job_args{-min_ic} = $min_ic;
$job_args{threshold} = $threshold;
$job_args{-z_cutoff} = $zscore_cutoff;
$job_args{-f_cutoff} = $fisher_cutoff;
$job_args{-ks_cutoff} = $ks_cutoff;
$job_args{-num_results} = $num_results;
$job_args{-sort_by} = $sort_by;

send_email(%job_args) if $email;

$logger->info("Finished analysis");

exit;


##########################################

#
# Ouput combined z-score/Fisher/KS-test results as HTML
#
sub write_results_html
{
    my ($tf_cluster_set, $results) = @_;
    
    my $result_type;
    if (defined $zscore_cutoff || defined $fisher_cutoff || defined $ks_cutoff)
    {
        $result_type = 'significant_hits';
    } else {
        $result_type = 'top_x_results';
    }

    my $warn_zero_bg_hits;
    foreach my $result (@$results) {
        if ($result->bg_gene_hits() == 0 && $result->t_gene_hits() > 0) {
            $warn_zero_bg_hits = 1;
            last;
        }
    }
    
    #my $tfcl_select_criteria = 'specific' if $tf_families;

    my $title = "oPOSSUM $heading";

    my $vars = {
        abs_htdocs_path         => ABS_HTDOCS_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_htdocs_path         => REL_HTDOCS_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        rel_htdocs_data_path    => REL_HTDOCS_DATA_PATH,
        version                 => VERSION,
        devel_version           => DEVEL_VERSION,
        jaspar_url              => JASPAR_URL,
        bg_color_class          => BG_COLOR_CLASS,
        heading                 => $heading,
        title                   => $title,
        section                 => 'Analysis Results',
        result_retain_days      => REMOVE_RESULTFILES_OLDER_THAN,
        low_matrix_ic           => LOW_MATRIX_IC,
        high_matrix_ic          => HIGH_MATRIX_IC,
        low_matrix_gc           => LOW_MATRIX_GC,
        high_matrix_gc          => HIGH_MATRIX_GC,
        low_seq_gc              => LOW_SEQ_GC,
        high_seq_gc             => HIGH_SEQ_GC,
        job_id                  => $job_id,
        num_t_seqs              => scalar @$t_seqs,
        num_bg_seqs             => scalar @$bg_seqs,
        t_seq_len               => $t_seq_len,
        bg_seq_len              => $bg_seq_len,
        t_seq_gc_content        => sprintf("%.3f", $t_seq_gc_content),
        bg_seq_gc_content       => sprintf("%.3f", $bg_seq_gc_content),
        tf_db                   => $tf_db,
        cl_db                   => $tfcl_db,
        #cl_select_criteria      => $tfcl_select_criteria,
        tf_cluster_set          => $tf_cluster_set,
        collection              => \@collections,
        tax_groups              => \@tax_groups,
        tf_families             => \@tf_families,
        min_ic                  => $min_ic,
        threshold               => $threshold,
        email                   => $email,
        results                 => $results,
        rel_results_dir         => $rel_results_dir,
        result_type             => $result_type,
        num_display_results     => $num_results,
        zscore_cutoff           => $zscore_cutoff,
        fisher_cutoff           => $fisher_cutoff,
        ks_cutoff               => $ks_cutoff,
        result_sort_by          => $sort_by,
        warn_zero_bg_hits       => $warn_zero_bg_hits,
        results_file            => RESULTS_TEXT_FILENAME,
        user_tf_family_file     => $user_tf_family_file,
        user_seq_file           => $user_seq_file,
        user_bg_seq_file        => $user_bg_seq_file,
        user_peak_pos_file      => $user_peak_pos_file,
        user_bg_peak_pos_file   => $user_bg_peak_pos_file,
        bg_seq_set_name         => $bg_seq_set_name,

        formatf                 => sub {
                                        my $dec = shift;
                                        my $f = shift;
                                        return ($f || $f eq '0')
                                            ? sprintf("%.*f", $dec, $f)
                                            : 'NA'
                                   },

        formatg                 => sub {
                                        my $dec = shift;
                                        my $f = shift;
                                        return ($f || $f eq '0')
                                            ? sprintf("%.*g", $dec, $f)
                                            : 'NA'
                                   },

        var_template            => "results_seq_tca.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    my $html_filename = "$results_dir/" . RESULTS_HTDOCS_FILENAME;
    
    open(OUT, ">$html_filename")
        || fatal("Could not create HTML results file $html_filename", %job_args);

    print OUT $output;

    close(OUT);

    return $html_filename;
}

#
# For each TFCluster/seq, write the details of the putative TFBS clusters
# out to text and html files.
#
sub write_tfbs_cluster_details
{
    my $logger = $job_args{-logger};
	$logger->info("Writing TFBS cluster details");

    # if results are truncated, show the relevant tf cluster ids only
    my $tfcl_ids;
    #my $t_gids = $ac->gene_ids;
    foreach my $result (@$cresults) {
        push @$tfcl_ids, $result->id;
    }

    foreach my $tfcl_id (@$tfcl_ids)
    {
        my $tfcl = $tf_cluster_set->get_tf_cluster($tfcl_id);
        my $tfcl_id = $tfcl->id();
        my $tf_ids = $tfcl->tf_ids;
        
        my $text_filename = sprintf "$abs_results_dir/c$tfcl_id.txt";
        my $html_filename = sprintf "$abs_results_dir/c$tfcl_id.html";

        my $t_seq_siteset = $t_cluster_seq_siteset->{$tfcl_id};

        if ($t_seq_siteset) {
            my @seq_ids = keys %$t_seq_siteset;
        
            write_tfbs_cluster_details_text(
                $text_filename,
                $tfcl, $tfcl_db,
                \@seq_ids,
                $t_seq_id_display_ids,
                $t_seq_siteset,
                %job_args
            );
            
            write_tfbs_cluster_details_html(
                $html_filename, $rel_results_dir,
                $tfcl, $tfcl_db,
                \@seq_ids,
                $t_seq_id_display_ids,
                $t_seq_siteset,
                %job_args
            ) if $web;
        }
    }
}
