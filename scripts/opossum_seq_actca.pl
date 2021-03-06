#!/usr/bin/env perl

=head1 NAME

opossum_seq_actca.pl

=head1 SYNOPSIS

  opossum_seq_actca.pl
      -d results_dir
      -s t_seq_file
      -b bg_seq_file
      -aid anchor_tf_id
      -dist site_distance
      [-tdb tf_database]
      [-cdb tf_cluster_database]
      [-fam tf_families] | [-famf tf_family_file]
      [-co collections] [-tax tax_groups] [-ic min_ic]
      [-th threshold]
      [-n num_results | -zcutoff cutoff -fcutoff cutoff]
      [-nh]
      [-sr sort_by]
      [-web]
      [-j job_id]
      [-m email]
	  [-usf user_t_seq_file]
	  [-ubsf user_bg_seq_file]
      [-ufamf user_tf_family_file]
      [-bss bg_seq_seq_name]

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

    -aid, -atfid tf_id
            Anchoring JASPAR TF ID.

    -dist site_distance
            Inter-binding site distance. This is the maximum distance
            between the binding sites of anchoring TF clusters and the
            binding sites of the other TF clusters of interest.

    -db, -tfdb db_name
            Specifies which TF database to use; default = JASPAR_2010.

    -cdb, -cldb cluster_db_name
            Specifies which TFBS cluster database to use;
            default = TFBS_clusters.

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

    -usf, -utsf user_t_seq_file
            User supplied target sequences upload file name.

    -ubsf user_bg_seq_file
            User supplied background sequences upload file name.

    -ufamf, -utff user_tf_family_file
            User supplied TF family upload file name.

    -bss bg_seq_set_name
            Background sequence set name.

=head1 DESCRIPTION

Take a set of peak sequences as provided in the input peaks file and a set of
control peak sequences. Specify an anchoring TF ID and the maximum inter-binding
site distance. Optionally take: 1) a subset of transcription factors (TFs) 
limited by external (JASPAR) database name and information content or taxonomic
supergroup or all TFs in the oPOSSUM database, 2) TFBS structural families to 
include in analysis, and 3) PWM score threshold.

Count the number of TFBS clusters located within the specified site distance 
from each of the anchor TFBS cluster found at the given PWM score threshold for
both the testcontrol set of peaks. Perform Fisher exact test and z-score 
analysis and output these results to the output file. Optionally write the
details of TFBS clusters found in test set to detailed TFBS cluster hits file.

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

use OPOSSUM::Include::SeqACTCAInclude;
use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::Analysis::Cluster::Zscore;
use OPOSSUM::Analysis::Cluster::Fisher;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::CombinedResultSet;

use Statistics::Distributions;

#use constant DEBUG          => 0;

my $help;
my $job_id;
my $results_dir;
my $web;
my $species;
my $t_seq_file;
my $bg_seq_file;
my $tf_db;
my $tfcl_db;
my @collections;
my @tax_groups;
my $tf_family_file;
my @tf_families;
my $anchor_tf_id;
my $max_site_dist;
my $min_ic;
my $threshold;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $nh;
my $email;
my $user_t_seq_file;
my $user_bg_seq_file;
my $user_tf_family_file;
my $bg_seq_set_name;
GetOptions(
    'dir|d=s'       => \$results_dir,
    'tsf|s=s'       => \$t_seq_file,
    'bsf|b=s'       => \$bg_seq_file,
    'tfdb|tdb|db=s' => \$tf_db,
    'cldb|cdb=s'    => \$tfcl_db,
    'atfid|aid=s'   => \$anchor_tf_id,
    'dist=i'        => \$max_site_dist,
    'co=s'          => \@collections,
    'tax=s'         => \@tax_groups,
    'ic=s'          => \$min_ic,
    'fam=s'         => \@tf_families,
    'famf|tff=s'    => \$tf_family_file,
    'th=s'          => \$threshold,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'sr=s'          => \$sort_by,
    'nh'            => \$nh,
    'm=s'           => \$email,
    'web'           => \$web,
    'job_id|j=s'    => \$job_id,
    'utsf|usf=s'    => \$user_t_seq_file,
    'ubsf=s'        => \$user_bg_seq_file,
    'ufamf|utff=s'  => \$user_tf_family_file,
    'bss=s'         => \$bg_seq_set_name,
    'help|h|?'      => \$help
);

if ($help) {
	pod2usage(
		-verbose	=> 1
	);
}

unless ($results_dir) {
    pod2usage(
        -msg        => "No results directory specified."
                     . " Use -help switch for detailed option list.\n",
        -verbose    => 1
    );
}

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
my $log_file = get_log_filename('opossum_seq_actca', $results_dir);

my $logger = get_logger();
if (DEBUG) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}

#my $layout = Log::Log4perl::Layout::PatternLayout->new("%M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("[%d] %p\t%m%n");

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

#
# Begin setting parameters
#
if (!$job_id) {
	$job_id = 'opossum_seq_actca';
}

my $message = "";
my $ok = 1;

my %job_args = (
	-job_id => $job_id,
	-email => $email,
	-logger => $logger
);

my $heading = sprintf(
	"Sequence-Based Anchored Combination TFBS Cluster Analysis"
);
$job_args{-heading} = $heading;

unless ($anchor_tf_id) {
    fatal("No anchor TF ID specified", %job_args);
}

unless ($t_seq_file) {
    fatal("No sequences file specified", %job_args);
}
unless (-e $t_seq_file) {
	fatal("Sequence file $t_seq_file does not exist", %job_args);
}

unless ($bg_seq_file) {
	fatal("No background sequences file specified", %job_args);
}
unless (-e $bg_seq_file) {
	fatal("Background sequence file $bg_seq_file does not exist", %job_args);
}

#
# Read in target sequences
#
$logger->info("Reading target sequences from $t_seq_file");
my $t_seqs = read_seqs($t_seq_file)
	|| fatal("Could not read sequences", %job_args);

check_seqs($t_seqs, "target", %job_args);

my ($t_seq_ids, $t_seq_id_seqs, $t_seq_id_display_ids) = id_seqs($t_seqs);
my ($t_seq_gc_content, $t_seq_len) = calc_seqs_total_gc_content($t_seqs);

#
# Read in background sequences
#
$logger->info("Reading background sequences from $bg_seq_file");
my $bg_seqs = read_seqs($bg_seq_file)
	|| fatal("Could not read sequences", %job_args);

check_seqs($bg_seqs, "background", %job_args);

my ($bg_seq_ids, $bg_seq_id_seqs, $bg_seq_id_display_ids) = id_seqs($bg_seqs);
my ($bg_seq_gc_content, $bg_seq_len) = calc_seqs_total_gc_content($bg_seqs);


#
# set optional parameters to default values if not provided by the user
#

my $collections_str;
if (@collections) {
    $collections_str = join(',', @collections);
	@collections = split /\s*,\s*/, $collections_str;

	unless (@collections) {
		fatal("Error parsing collections", %job_args);
	}
} else {
    $collections_str = 'CORE';
    push @collections, 'CORE';
}

my $tf_families_str;
if (@tf_families) {
    $tf_families_str = join(',', @tf_families);
	@tf_families = split /\s*,\s*/, $tf_families_str;

	unless (@tf_families) {
		fatal("Error parsing tf_families", %job_args);
	}
} elsif ($tf_family_file) {
    $logger->info("Reading TFBS family file");

   @tf_families = @{read_tf_families_from_file($tf_family_file, %job_args)};

    unless (@tf_families) {
        fatal("No TFBS families read from file $tf_family_file", %job_args);
    }

    $tf_families_str = join(',', @tf_families);
}

my $tax_groups_str;
if (@tax_groups) {
    $tax_groups_str = join(',', @tax_groups);
	@tax_groups = split /\s*,\s*/, $tax_groups_str;

	unless (@tax_groups) {
		fatal("Error parsing tax groups", %job_args);
	}
}

$threshold = DFLT_THRESHOLD . "%" if !$threshold;
$max_site_dist = DFLT_INTER_BINDING_DIST if !$max_site_dist;


#
# JASPAR / TF parameter settings
# Also, TFBSCluster settings
# No TF profile file possible here
#

$tf_db = JASPAR_DB_NAME if !$tf_db;

# if this is default, JASPAR has to be default as well
$tfcl_db = TFBS_CLUSTER_DB_NAME if !$tfcl_db;

if ($tfcl_db eq TFBS_CLUSTER_DB_NAME and $tf_db ne JASPAR_DB_NAME) {
	fatal("\nDefault TFBS cluster DB requires the default JASPAR DB",
		%job_args);
}

#
# Connect to JASPAR database
#
$logger->info("Fetching matrices from JASPAR");

my $jdb = jaspar_db_connect($tf_db)
	|| fatal("Could not connect to JASPAR database $tf_db", %job_args);

my %get_matrix_args = (
	-matrixtype => 'PFM'
);

$get_matrix_args{-collection} = \@collections if @collections;
$get_matrix_args{-tax_group} = \@tax_groups if @tax_groups;
$get_matrix_args{-min_ic} = $min_ic || DFLT_CORE_MIN_IC;
$get_matrix_args{-family} = \@tf_families if @tf_families;

$logger->info("Fetching TF matrices");
my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

unless ($matrix_set && $matrix_set->size > 0) {
    fatal("Error fetching TF profiles from JASPAR DB", %job_args);
}

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

$logger->info("Fetching anchoring TF matrix");
my $anchor_tf = $jdb->get_Matrix_by_ID($anchor_tf_id, 'PFM');
unless ($anchor_tf) {
    fatal("Error fetching anchoring TF profile from JASPAR DB", %job_args);
}

#
# retrieve the TFBS cluster set
#

my $cdb = tfbs_cluster_db_connect($tfcl_db);
if (!$cdb) {
	fatal("Could not connect to $tfcl_db", %job_args);
}

my $tfca = $cdb->get_TFClusterAdaptor;
fatal("Error fetching TFClusterAdaptor") if !$tfca;

my $tf_cluster_set = fetch_tf_cluster_set($cdb, %get_matrix_args);
unless ($tf_cluster_set and $tf_cluster_set->size > 0) {
	fatal("Error fetching TFBS clusters", %job_args);
}

$logger->info("Fetching anchoring TFBS cluster");
my $anchor_cluster = $tfca->fetch_by_tf_id($anchor_tf_id);

#my $tf_ids = $tf_set->ids();
#my $tfcl_ids = $tf_cluster_set->ids();


#
# Search for anchor cluster - TFBS cluster pairs
#

$logger->info("Fetching target sequence anchored TFBS clusters");

my ($t_cluster_seq_sites, $t_cluster_seq_sitepairs) = 
	anchored_tf_cluster_set_search_seqs(
		$tf_cluster_set, $tf_set, $anchor_cluster, $max_site_dist,
		$t_seq_id_seqs, $threshold, %job_args
	);

unless ($t_cluster_seq_sites) {
	$message = "No anchored TFBS clusters found in target sequences";
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

my $cresults;
if ($ok) {
	$logger->info("Searching background sequences for TFBS clusters");
	my ($bg_cluster_seq_sites) = 
		anchored_tf_cluster_set_search_seqs(
			$tf_cluster_set, $tf_set, $anchor_cluster, $max_site_dist,
			$bg_seq_id_seqs, $threshold, %job_args
	);

	$logger->info("No anchored TFBS clusters found in control sequences")
		unless $bg_cluster_seq_sites;

	$logger->info("Computing target sequence TFBS cluster counts");
	my $t_counts = compute_site_cluster_counts(
		$tf_cluster_set, $t_seq_ids, $t_cluster_seq_sites
	);

	$logger->info("Computing background sequence TFBS cluster counts");
	my $bg_counts = compute_site_cluster_counts(
		$tf_cluster_set, $bg_seq_ids, $bg_cluster_seq_sites
	);

	#
	# Score calculations
	#

	my $fisher = OPOSSUM::Analysis::Cluster::Fisher->new();
	fatal("Error initializing Fisher analysis", %job_args) unless $fisher;

	$logger->info("Computing Fisher scores");
	my $fresult_set = $fisher->calculate_Fisher_probability(
    	$bg_counts,
    	$t_counts
	);
	fatal("Error performing Fisher analysis", %job_args) unless $fresult_set;

	my $zscore = OPOSSUM::Analysis::Cluster::Zscore->new();
	fatal("Error initializing z-score analysis", %job_args) unless $zscore;

	$logger->info("Computing Z-scores");
	my $zresult_set = $zscore->calculate_Zscore(
    	$bg_counts,	$t_counts, $bg_seq_len, $t_seq_len
	);
	fatal("Error computing z-score", %job_args) unless $zresult_set;

	#
	# Use new OPOSSUM::Analysis::Cluster::CombinedResultSet to combine 
	# Fisher and Z-score result sets
	#
	$logger->info("Combining Fisher and Z-scores");
	my $cresult_set = OPOSSUM::Analysis::Cluster::CombinedResultSet->new(
    	-fisher_result_set  => $fresult_set,
    	-zscore_result_set  => $zresult_set
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

	if (defined $sort_by) {
    	if ($sort_by =~ /^fisher/) {
        	$sort_by = 'fisher_p_value';
    	} elsif ($sort_by =~ /^z_score/ || $sort_by =~ /^z-score/) {
        	$sort_by = 'zscore';
    	}

    	$result_params{-sort_by} = $sort_by;

        #
        # Fisher scores are now -ln() values so all results are reverse sorted
        # (largest value first).
        # DJA 14/03/2011
        #
        #if ($sort_by eq 'zscore') {
        #    # Sort z-score from highest to lowest
        #    $result_params{-reverse} = 1;
        #}
    	$result_params{-reverse} = 1;
	}

	$logger->info("Getting filtered/sorted result list");
	$cresults = $cresult_set->get_list(%result_params);

	unless ($cresults) {
    	$message = "No anchored TFBS clusters scored above the selected "
			. "Z-score/Fisher thresholds";
    	$logger->info($message);
    	$ok = 0;
        #
        # XXX
        # The ok/message stuff is not handled properly later, so just fatal it
        # for now.
        # DJA 2012/05/03
        #
        fatal($message, %job_args);
	} else {
		if ($web) {
			$logger->info("Writing HTML results");
			write_results_html();
		}

    	$logger->info("Writing text results");
		my $out_file = "$abs_results_dir/" . RESULTS_TEXT_FILENAME;
    	write_results_text($out_file, $cresults, $tf_cluster_set, %job_args);

		if (!$nh) {
    		$logger->info("Writing TFBS cluster details");
    		write_tfbs_cluster_details();
		}
	}
}

$job_args{-web} = $web;
$job_args{-abs_results_dir} = $abs_results_dir;
$job_args{-rel_results_dir} = $rel_results_dir;
$job_args{-user_t_file} = $user_t_seq_file;
$job_args{-user_bg_file} = $user_bg_seq_file;
$job_args{-t_num} = scalar @$t_seqs if $t_seqs;
$job_args{-bg_num} = scalar @$bg_seqs if $bg_seqs;
$job_args{-tf_db} = $tf_db;
$job_args{-collections} = $collections_str if $collections_str;
$job_args{-tax_groups} = $tax_groups_str if $tax_groups_str;
$job_args{-families} = $tf_families_str if $tf_families_str;
$job_args{-min_ic} = $min_ic;
$job_args{-threshold} = $threshold;
$job_args{-z_cutoff} = $zscore_cutoff;
$job_args{-f_cutoff} = $fisher_cutoff;
$job_args{-num_results} = $num_results;
$job_args{-sort_by} = $sort_by;

send_email(%job_args) if $email;

$logger->info("Finished analysis");

exit;

#################################################

#
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
	my ($filename, $results, $tf_cluster_set, %job_args) = @_;

    return unless $results && $results->[0];

    open(FH, ">$filename")
        || fatal("Could not create analysis results file $filename", %job_args);

    $logger->info("Writing analysis results to $filename\n");

    printf FH "TFBS cluster name\tClass\tFamily\tTarget sequence hits\tTarget sequence non-hits\tBackground sequence hits\tBackground sequence non-hits\tTarget TFBS cluster hits\tTarget TFBS cluster nucleotide rate\tBackground TFBS cluster hits\tBackground TFBS cluster nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tfcl = $tf_cluster_set->get_tf_cluster($result->id());

        printf FH 
            "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",
            $tfcl->name(),
            $tfcl->class() || 'NA',
            $tfcl->family() || 'NA',
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_cluster_hits() || 0,
            defined $result->t_cluster_rate()
                ? sprintf("%.3f", $result->t_cluster_rate()) : 'NA',
            $result->bg_cluster_hits() || 0,
            defined $result->bg_cluster_rate()
                ? sprintf("%.3f", $result->bg_cluster_rate()) : 'NA',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'NA',
            defined $result->fisher_p_value()
                ? sprintf("%.3g", $result->fisher_p_value()) : 'NA';
    }
    close(FH);
}


#
# Ouput combined z-score/Fisher results as HTML
#
sub write_results_html
{
	my $tfcl_ids = $tf_cluster_set->ids();

    my $warn_zero_bg_hits = 0;
    foreach my $result (@$cresults) {
        if ($result->bg_gene_hits() == 0) {
            $warn_zero_bg_hits = 1;
            last;
        }
    }

    my $result_type;
    if (defined $zscore_cutoff || defined $fisher_cutoff) {
        $result_type = 'significant_hits';
    } else {
        $result_type = 'top_x_results';
    }

	my $tfcl_select_criteria = 'specific' if @tf_families;

    my $title = "oPOSSUM $heading";

	if ($threshold =~ /%$/) {
		$threshold =~ s/%//;
		$threshold /= 100;
	}

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
		rel_htdocs_data_path=> REL_HTDOCS_DATA_PATH,
        jaspar_url          => JASPAR_URL,
        title               => $title,
        heading             => $heading,
        section             => 'Analysis Results',
        bg_color_class      => BG_COLOR_CLASS,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        result_retain_days  => REMOVE_RESULTFILES_OLDER_THAN,
		#low_matrix_ic       => LOW_MATRIX_IC,
		#high_matrix_ic      => HIGH_MATRIX_IC,
		#low_matrix_gc       => LOW_MATRIX_GC,
		#high_matrix_gc      => HIGH_MATRIX_GC,
		low_seq_gc          => LOW_SEQ_GC,
		high_seq_gc         => HIGH_SEQ_GC,
        job_id              => $job_id,
        num_t_seqs          => $t_seqs ? scalar @$t_seqs : 0,
        num_bg_seqs         => $bg_seqs ? scalar @$bg_seqs : 0,
		t_seq_len			=> $t_seq_len,
		bg_seq_len			=> $bg_seq_len,
		t_seq_gc_content	=> sprintf("%.3f", $t_seq_gc_content),
		bg_seq_gc_content	=> sprintf("%.3f", $bg_seq_gc_content),
        tf_db               => $tf_db,
		cl_db				=> $tfcl_db,
        cl_select_criteria  => $tfcl_select_criteria,
		tf_set				=> $tf_set,
        tf_cluster_set      => $tf_cluster_set,
        anchor_tf_name      => $anchor_tf->name,
		anchor_cluster		=> $anchor_cluster,
        max_site_distance   => $max_site_dist,
		collections         => \@collections,
		tax_groups          => \@tax_groups,
		min_ic              => $min_ic,
		tf_families			=> \@tf_families,
        threshold           => $threshold,
        results             => $cresults,
        rel_results_dir     => $rel_results_dir,
        result_type         => $result_type,
        num_display_results => $num_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $sort_by,
        warn_zero_bg_hits   => $warn_zero_bg_hits,
        results_file        => RESULTS_TEXT_FILENAME,
        message             => $message,
		user_t_seq_file		=> $user_t_seq_file,
		user_bg_seq_file	=> $user_bg_seq_file,
        user_tf_family_file => $user_tf_family_file,
		email				=> $email,

        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'NA'
                               },

        formatg             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*g", $dec, $f)
                                        : 'NA'
                               },

        var_template        => "results_seq_actca.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    my $html_filename = "$results_dir/" . RESULTS_HTDOCS_FILENAME;

    open(OUT, ">$html_filename")
        || fatal("Could not create HTML results file $html_filename",
				%job_args);

    print OUT $output;

    close(OUT);

    $logger->info("Wrote HTML formatted results to $html_filename");

    return $html_filename;
}

#
# Write the details of the putative TFBS clusters for each TFCluster/sequence. 
# Create a html and text file for each TFCluster.
#
sub write_tfbs_cluster_details
{
	my $anchor_cluster_tf_ids = $anchor_cluster->tf_ids();

	my $cluster_ids = $tf_cluster_set->ids();

	foreach my $cl_id (@$cluster_ids) 
	{
		my $seq_sitepairs = $t_cluster_seq_sitepairs->{$cl_id};

		if ($seq_sitepairs) {
			my $tfcl = $tf_cluster_set->get_tf_cluster($cl_id);
			my @seq_ids = keys %$seq_sitepairs;

			my $text_filename = "$abs_results_dir/c$cl_id.txt";
			my $html_filename = "$abs_results_dir/c$cl_id.html";

	        write_tfbs_cluster_details_text(
    	        $text_filename, 
				$tfcl, $anchor_cluster, 
				\@seq_ids, $t_seq_id_display_ids, $seq_sitepairs,
				%job_args
        	);

        	write_tfbs_cluster_details_html(
            	$html_filename, $rel_results_dir,
				$tfcl, $anchor_cluster,
				\@seq_ids, $t_seq_id_display_ids, $seq_sitepairs,
				%job_args
        	) if $web;
		}
    }
}

#
# Write the details of the putative TFBS clusters for each TF/sequence. Create a
# text file for each TFBS cluster.
#
sub write_tfbs_cluster_details_text
{
    my ($filename, 
		$tfcl, $anchor_cluster,
		$seq_ids, $seq_id_display_ids, $seq_sitepairs,
		%job_args
	) = @_;


	my $logger = $job_args{-logger};

    open(FH, ">$filename") || fatal(
        "Could not create TFBS cluster details text file $filename", %job_args
    );

	my $tfcl_id = $tfcl->id();
	my $tfcl_name = $tfcl->name();
	my $anchor_id = $anchor_cluster->id();
	my $anchor_name = $anchor_cluster->name();

    $logger->info("Writing '$tfcl_name' TFBS cluster details to $filename");

	printf FH "$anchor_name\n\n";

    printf FH "TFBS Cluster ID:\t%s\n", $anchor_cluster->id();
    printf FH "Class:    \t%s\n", $anchor_cluster->class() || 'NA';
    printf FH "Family:   \t%s\n", $anchor_cluster->family() || 'NA';

    printf FH "\n\n";

    printf FH "%s\n\n", $tfcl->name();

    printf FH "TFBS Cluster ID:\t%s\n", $tfcl->id();
    printf FH "Class:    \t%s\n", $tfcl->class() || 'NA';
    printf FH "Family:   \t%s\n", $tfcl->family() || 'NA';

    printf FH "\n\nc%s - c%s Binding Site Cluster Combinations\n\n",
        $anchor_cluster->name(), $tfcl->name();

    printf FH "Sequence ID\tAnchoring TFBS Cluster\tStart\tEnd\tStrand\tScore\tTFBS Cluster Sequence\tAnchored TFBS Cluster\tStart\tEnd\tStrand\tScore\tTFBS Cluster Sequence\tDistance\n";

    foreach my $seq_id (@$seq_ids) 
	{
		my $sitepairs = $seq_sitepairs->{$seq_id};
		next unless $sitepairs && $sitepairs->[0];

		my $display_id = $seq_id_display_ids->{$seq_id};
        
        my $first = 1;
		printf FH "%-31s", $display_id;

        foreach my $sitepair (@$sitepairs) 
		{
            my $anchor_site	= $sitepair->{anchor_site};
            my $cl_site     = $sitepair->{cluster_site};
            my $distance    = $sitepair->{distance};

			printf FH "\t%s\t%d\t%d\t%s\t%.1f%%\t%s\t%s\t%d\t%d\t%s\t%.1f%%\t%s\t%d\n",
				$anchor_name,
				$anchor_site->start(),
				$anchor_site->end(),
				$anchor_site->strand() == 1 ? '+' : '-',
				$anchor_site->rel_score() * 100,
				$anchor_site->seq,
				$tfcl->name(),
				$cl_site->start(),
				$cl_site->end(),
				$cl_site->strand() == 1 ? '+' : '-',
				$cl_site->rel_score() * 100,
				$cl_site->seq,
				$distance;
		}

		print FH "\n";
	}

	close(FH);
}


#
# Write the details of the putative TFBSs for each TF/sequence. Create an
# HTML file for each TF.
#
sub write_tfbs_cluster_details_html
{
    my ($filename, $rel_results_dir, 
		$tfcl, $anchor_cluster, 
		$seq_ids, $seq_id_display_ids, $seq_sitepairs,
		%job_args
	) = @_;

    my $tfcl_name = $tfcl->name();
    my $tfcl_id   = $tfcl->id();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS cluster details html file $filename", %job_args
    );

    $logger->info("Writing '$tfcl_name' TFBS cluster details to $filename");

    my $heading = "Sequence-Based Anchored Combination TFBS Cluster Analysis";

    my $title = "oPOSSUM $heading";

    my $section = sprintf(
        "c%s - c%s Binding Site Cluster Combinations",
        $anchor_cluster->id(), $tfcl_id
    );

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        bg_color_class      => BG_COLOR_CLASS,
        title               => $title,
        heading             => $heading,
        section             => $section,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
		#low_matrix_ic       => LOW_MATRIX_IC,
		#high_matrix_ic      => HIGH_MATRIX_IC,
		#low_matrix_gc       => LOW_MATRIX_GC,
		#high_matrix_gc      => HIGH_MATRIX_GC,
        low_seq_gc          => LOW_SEQ_GC,
        high_seq_gc         => HIGH_SEQ_GC,

        jaspar_url          => JASPAR_URL,

        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'NA'
                               },

		tf_db               => $tf_db,
		cl_db				=> $tfcl_db,
        anchor_tf_name      => $anchor_tf->name,
		anchor_cluster		=> $anchor_cluster,
        tf_cluster          => $tfcl,
		seq_ids				=> $seq_ids,
		seq_id_display_ids	=> $seq_id_display_ids,
		seq_sitepairs		=> $seq_sitepairs,
        rel_results_dir     => $rel_results_dir,
        tfbs_cluster_details_file   => "c$tfcl_id.txt",
        var_template        => "tfbs_cluster_details_seq_actca.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}

