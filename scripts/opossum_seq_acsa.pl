#!/usr/bin/env perl

=head1 NAME

opossum_seq_acsa.pl

=head1 SYNOPSIS

  opossum_seq_acsa.pl
      -d results_dir
      -s t_seq_file
      -b bg_seq_file
      -dist site_distance
      [-atf anchor_matrix_file] | [-aid anchor_tf_id]
      ([-mf matrix_file]
          | [-ids tf_ids]
          | [-co collections] [-tax tax_groups] [-ic min_ic])
      [-th threshold]
      [-n num_results] | [-zcutoff cutoff -fcutoff cutoff]
      [-sr sort_by]
      [-web]
      [-j job_id]
      [-m email]
      [-usf user_seq_file]
      [-ubsf user_bg_seq_file]
      [-uamf user_anchor_matrix_file]
      [-umf user_matrix_file]
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

    -dist site_distance
            Inter-binding site distance. This is the maximum distance
            between binding sites of the anchoring TF and binding sites of
            the other TFs of interest.

    -aid, -atfid tf_id
            Anchoring TF ID.

    -amf, -atf, -atff anchor_matrix_file
            File containing the anchoring TFBS profile matrix. May be a
            position frequency or position weight matrix.

    -mf, -tff matrix_file
            File containing one or more TFBS profile matrices. May be
            position frequency or position weight matrices. If specified,
            this options takes presedence over any of the -tfids, -c, -tax
            and -ic options below.

    -ids, -tfids tf_ids
            Specify one of more JASPAR TFBS profile matrix IDs to include
            use in the analysis. If this option is given it overrides any
            combination of the -co, -tax and -ic options below.

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

    -uamf, -uatff user_anchor_matrix_file
            User supplied anchor matrix file name.

    -umf, -utff user_matrix_file
            User supplied matrix file name.

    -usf, -utsf user_seq_file
            User supplied target sequences upload file name.

    -ubsf user_bg_seq_file
            User supplied background sequences upload file name.

    -bss bg_seq_set_name
            Background sequence set name.

=head1 DESCRIPTION

Take a set of peak sequences as provided in the input peaks file and a set of 
control peak sequences. Specify an anchoring TF ID or a file containing the 
anchoring TFBS matrix and the maximum inter-binding site distance. Optionally, 
take: 1) a subset of transcription factors (TFs) either specified in an input 
file, or limited by external (JASPAR) database name and information content or
taxonomic supergroup or all TFs in the oPOSSUM database, and 2) PWM score 
threshold.

Count the number of TFBSs located within the specified site distance from each
of the anchor TF binding site found at the given PWM score threshold for both
the test and control set of peaks. Perform Fisher exact test and z-score 
analysis and output these results to the output file. Optionally write the 
details of TFBSs found in test set to detailed TFBS hits file.

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

use OPOSSUM::Include::SeqACSAInclude;
use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::CombinedResultSet;
use OPOSSUM::Plot::ScoreVsGC;

use Statistics::Distributions;

#use constant DEBUG          => 0;

my $help;
my $job_id;
my $results_dir;
my $web;
my $t_seq_file;
my $bg_seq_file;
my @collections;
my @tax_groups;
my @tf_ids;
my $matrix_file;
my $anchor_tf_id;
my $anchor_matrix_file;
my $max_site_dist;
my $min_ic;
my $threshold;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $nh;
my $email;
my $user_seq_file;
my $user_bg_seq_file;
my $user_anchor_matrix_file;
my $user_matrix_file;
my $bg_seq_set_name;
GetOptions(
    'dir|d=s'       => \$results_dir,
    'tsf|s=s'       => \$t_seq_file,
    'bsf|b=s'       => \$bg_seq_file,
    'atfid|aid=s'   => \$anchor_tf_id,
	'amf|atff|atf=s'=> \$anchor_matrix_file,
    'dist=i'        => \$max_site_dist,
    'co=s'          => \@collections,
    'tax=s'         => \@tax_groups,
    'ic=s'          => \$min_ic,
    'tfids|ids=s'   => \@tf_ids,
    'mf|tff|tf=s'   => \$matrix_file,
    'th=s'          => \$threshold,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'sr=s'          => \$sort_by,
    'nh'            => \$nh,
    'm=s'           => \$email,
    'web'           => \$web,
    'job_id|j=s'    => \$job_id,
    'utsf|usf=s'    => \$user_seq_file,
    'ubsf=s'        => \$user_bg_seq_file,
    'uamf|uatff=s'  => \$user_anchor_matrix_file,
    'umf|utff=s'    => \$user_matrix_file,
    'bss=s'         => \$bg_seq_set_name,
    'help|h|?'      => \$help
);

#my $exclude_single_hits = $no_exclude_single_hits ? 0 : 1;

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

$job_id = 'opossum_seq_acsa' unless $job_id;
$num_results = 'All' unless $num_results;

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

unless (-d $abs_results_dir) {
    die "Results directory $abs_results_dir does not exist\n";
}

#
# Initialize logging
#
my $log_file = get_log_filename("opossum_seq_acsa", $results_dir);

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

my $heading = "Sequence-Based Anchored Combination Site Analysis";

my %job_args = (
	-job_id => $job_id,
	-heading => $heading,
	-email => $email,
	-logger => $logger
);

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
	fatal("Background sequences file $bg_seq_file does not exist", %job_args);
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
	|| fatal("Could not read background sequences", %job_args);

check_seqs($bg_seqs, "background", %job_args);

my ($bg_seq_ids, $bg_seq_id_seqs, $bg_seq_id_display_ids) = id_seqs($bg_seqs);
my ($bg_seq_gc_content, $bg_seq_len) = calc_seqs_total_gc_content($bg_seqs);

#
# JASPAR / TF parameter settings
# 

#
# A provided TFBS matrix file overrides everything else
#
my $tf_ids_str;
my $collections_str;
my $tax_groups_str;
unless ($matrix_file) {
    if (@tf_ids) {
        $tf_ids_str = join(',', @tf_ids);
        @tf_ids = split /\s*,\s*/, $tf_ids_str;

        unless (@tf_ids) {
            fatal("Error parsing TF IDs", %job_args);
        }
    } else {
        if (@collections) {
            $collections_str = join(',', @collections);
            @collections = split /\s*,\s*/, $collections_str;

            unless (@collections) {
                fatal("Error parsing collections", %job_args);
            }
        } else {
            push @collections, 'CORE';
            $collections_str = 'CORE';
        }

        if (@tax_groups) {
            $tax_groups_str = join(',', @tax_groups);
            @tax_groups = split /\s*,\s*/, $tax_groups_str;

            unless (@tax_groups) {
                fatal("Error parsing tax groups", %job_args);
            }
        }
    }
}

$threshold = DFLT_THRESHOLD . "%" if !$threshold;

$max_site_dist = DFLT_INTER_BINDING_DIST if !$max_site_dist;

#
# this is different from gene-based...
# need to go over in detail over different possibilities
#
my $tf_db;
my $jdb;
unless ($matrix_file && $anchor_matrix_file) {
    $tf_db = JASPAR_DB_NAME;

    $jdb = jaspar_db_connect($tf_db)
        || fatal("Could not connect to JASPAR database $tf_db", %job_args);
}

my $matrix_set;
my $tf_select_criteria;
if ($matrix_file) {
	#
	# 1. User-supplied TF profile file
	# 
	$tf_select_criteria = 'specific';

    $logger->info("Reading matrices from $matrix_file");

	unless (-e $matrix_file) {
		fatal("Specified input TF profile file $matrix_file does not exist",
			%job_args);
	}

	$matrix_set = read_matrices($matrix_file);
	unless ($matrix_set && $matrix_set->size > 0) {
		fatal("Error reading TF profiles from $matrix_file");
	}
} else {
	#
	# 2. Connect to JASPAR database and retrieve the matrices
	#
	$logger->info("Fetching matrices from JASPAR");

	my %get_matrix_args = (
		-matrixtype => 'PFM'
	);

	if (@tf_ids) {
		# This takes precedence over tax groups and min IC
		$tf_select_criteria = 'specific';
		$get_matrix_args{-ID} = \@tf_ids;
	} else {
		$tf_select_criteria = 'min_ic';
		$get_matrix_args{-collection} = \@collections if @collections;
		$get_matrix_args{-tax_group} = \@tax_groups if @tax_groups;
		$get_matrix_args{-min_ic} = $min_ic || DFLT_CORE_MIN_IC;
	}

	$logger->info("Fetching TF matrices");
	$matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

	unless ($matrix_set && $matrix_set->size > 0) {
		fatal("Error fetching TF profiles from JASPAR DB", %job_args);
	}
}

#
# Retrieve the anchor matrix
#
my $anchor_matrix;
if ($anchor_matrix_file) {
	$logger->info("Reading anchor matrix from $anchor_matrix_file");
	unless (-e $anchor_matrix_file) {
		fatal("Specified input anchor TFBS profile matrix file"
            . " $anchor_matrix_file does not exist", %job_args);
	}

	my $anchor_matrix_set = read_matrices($anchor_matrix_file);
	unless ($anchor_matrix_set && $anchor_matrix_set->size > 0) {
		fatal("Error reading anchor TFBS profile matrix from"
            . " $anchor_matrix_file");
	}
	
	if ($anchor_matrix_set->size > 1) {
		$logger->warn(
			"More than one matrix in anchor TFBS profile matrix file "
			. "$anchor_matrix_file - using first matrix found");
	}
	
	my $iter = $anchor_matrix_set->Iterator();
	$anchor_matrix = $iter->next();
} else {
	$logger->info("Fetching anchor matrix from JASPAR");
	
	$anchor_matrix = $jdb->get_Matrix_by_ID($anchor_tf_id, 'PFM');
	unless ($anchor_matrix) {
    	fatal("Error fetching anchoring TF profile from JASPAR DB", %job_args);
	}
}

#
# Search for anchor matrix - TF sitepairs
#

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

# t_tf_seq_sitepairs will be used for writing sitepair details later on
$logger->info("Searching target sequences for TFBSs");
my ($t_tf_seq_sites, $t_tf_seq_sitepairs) = 
	anchored_tf_set_search_seqs(
		$tf_set, $anchor_matrix, $max_site_dist, 
		$t_seq_id_seqs, $threshold, %job_args
	);

unless ($t_tf_seq_sites) {
	$message = "No anchored TFBSs found in target sequences";
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
if ($ok)
{
	$logger->info("Searching background sequences for TFBSs");
	my ($bg_tf_seq_sites) = 
		anchored_tf_set_search_seqs(
			$tf_set, $anchor_matrix, $max_site_dist,
			$bg_seq_id_seqs, $threshold, %job_args
		);

	$logger->info("No anchored TFBSs fofund in control sequences")
		unless $bg_tf_seq_sites;

	$logger->info("Computing target sequence TFBS counts");
	my $t_counts = compute_site_counts($tf_set, $t_seq_ids, $t_tf_seq_sites);

	$logger->info("Computing background sequence TFBS counts");
	my $bg_counts = compute_site_counts(
		$tf_set, $bg_seq_ids, $bg_tf_seq_sites
	);


	#
	# Score calculations
	#

	my $fisher = OPOSSUM::Analysis::Fisher->new();
	fatal("Error initializing Fisher analysis", %job_args) unless $fisher;

	$logger->info("Computing Fisher scores");
	my $fresult_set = $fisher->calculate_Fisher_probability(
    	$bg_counts,
    	$t_counts
	);
	fatal("Error performing Fisher analysis", %job_args) unless $fresult_set;

	my $zscore = OPOSSUM::Analysis::Zscore->new();
	fatal("Error initializing z-score analysis", %job_args) unless $zscore;

	$logger->info("Computing Z-scores");
	my $zresult_set = $zscore->calculate_Zscore(
    	$bg_counts,
    	$t_counts,
   	 	$bg_seq_len,
    	$t_seq_len,
    	$tf_set
	);
	fatal("Error computing z-score", %job_args) unless $zresult_set;

	#
	# Use new OPOSSUM::Analysis::CombinedResultSet to combine Fisher and 
	# Z-score result sets
	#
	my $cresult_set = OPOSSUM::Analysis::CombinedResultSet->new(
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
    	$message = "No anchored TFBSs scored above the selected Z-score/Fisher"
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
	} else {

		if ($web) {
			$logger->info("Writing HTML results");
			write_results_html();
		}

    	$logger->info("Writing text results");
		my $out_file = "$abs_results_dir/" . RESULTS_TEXT_FILENAME;
    	write_results_text($out_file, $cresults, $tf_set, $tf_db, %job_args);

		if (!$nh) {
    		$logger->info("Writing TFBS details");
    		write_tfbs_details($t_tf_seq_sitepairs, $t_seq_id_display_ids, 
				$tf_set, $anchor_matrix, $tf_db, 
				$abs_results_dir, $rel_results_dir, $web, %job_args);
		}
	}
}

$logger->info("Plotting scores vs. profile \%GC content");

my $plotter = OPOSSUM::Plot::ScoreVsGC->new();
unless ($plotter) {
    $logger->error("Could not initialize plotting");
} else {
    my $plot_err;

    my $z_plot_file = "$abs_results_dir/" . ZSCORE_PLOT_FILENAME;
    unless($plotter->plot(
        $cresults, $tf_set, 'Z', ZSCORE_PLOT_SD_FOLD, $z_plot_file,
        \$plot_err
    )) {
        $logger->error("Could not plot Z-scores vs. GC content. $plot_err");
    }

    my $fisher_plot_file = "$abs_results_dir/" . FISHER_PLOT_FILENAME;
    unless($plotter->plot(
        $cresults, $tf_set, 'Fisher', FISHER_PLOT_SD_FOLD, $fisher_plot_file,
        \$plot_err
    )) {
        $logger->error(
            "Could not plot Fisher scores vs. GC content. $plot_err"
        );
    }
}

$job_args{-web} = $web;
$job_args{-abs_results_dir} = $abs_results_dir;
$job_args{-rel_results_dir} = $rel_results_dir;
$job_args{-user_t_file} = $user_seq_file;
$job_args{-user_bg_file} = $user_bg_seq_file;
$job_args{-t_num} = scalar @$t_seqs if $t_seqs;
$job_args{-bg_num} = scalar @$bg_seqs if $bg_seqs;
$job_args{-tf_db} = $tf_db;
$job_args{-collections} = $collections_str if $collections_str;
$job_args{-tax_groups} = $tax_groups_str if $tax_groups_str;
$job_args{-tf_ids} = $tf_ids_str if $tf_ids_str;
$job_args{-tf_ids} = $matrix_file if $matrix_file;
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
# Ouput combined z-score/Fisher results as HTML
#
sub write_results_html
{
	my $tf_ids = $tf_set->ids();

    my $warn_zero_bg_gene_hits = 0;
    foreach my $result (@$cresults) {
        if ($result->bg_gene_hits() == 0 && $result->t_gene_hits() > 0) {
            $warn_zero_bg_gene_hits = 1;
            last;
        }
    }

    my $result_type;
    if (defined $zscore_cutoff || defined $fisher_cutoff) {
        $result_type = 'significant_hits';
    } else {
        $result_type = 'top_x_results';
    }

    my $title = "oPOSSUM $heading";

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        jaspar_url          => JASPAR_URL,
        title               => $title,
        heading             => $heading,
        section             => 'Analysis Results',
        bg_color_class      => BG_COLOR_CLASS,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        result_retain_days  => REMOVE_RESULTFILES_OLDER_THAN,
        low_matrix_ic       => LOW_MATRIX_IC,
        high_matrix_ic      => HIGH_MATRIX_IC,
        low_matrix_gc       => LOW_MATRIX_GC,
        high_matrix_gc      => HIGH_MATRIX_GC,
        low_seq_gc          => LOW_SEQ_GC,
        high_seq_gc         => HIGH_SEQ_GC,
		job_id				=> $job_id,
		num_t_seqs			=> scalar @$t_seqs,
		num_bg_seqs			=> scalar @$bg_seqs,
		t_seq_len			=> $t_seq_len,
		bg_seq_len			=> $bg_seq_len,
		t_seq_gc_content	=> sprintf("%.3f", $t_seq_gc_content),
		bg_seq_gc_content	=> sprintf("%.3f", $bg_seq_gc_content),
        tf_db               => $tf_db,
		matrix_file         => $matrix_file,
		tf_ids				=> \@tf_ids,
        tf_set              => $tf_set,
		tf_select_criteria	=> $tf_select_criteria,
        collections         => \@collections,
        collection          => $collections_str,
        tax_groups          => \@tax_groups,
        anchor_matrix       => $anchor_matrix,
        max_site_distance   => $max_site_dist,
        min_ic              => $min_ic,
        threshold           => $threshold,
        results             => $cresults,
        rel_results_dir     => $rel_results_dir,
        result_type         => $result_type,
        num_display_results => $num_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $sort_by,
        warn_zero_bg_gene_hits  => $warn_zero_bg_gene_hits,
        results_file        => RESULTS_TEXT_FILENAME,
        zscore_plot_file    => ZSCORE_PLOT_FILENAME,
        fisher_plot_file    => FISHER_PLOT_FILENAME,
        message             => $message,
		user_seq_file		=> $user_seq_file,
		user_bg_seq_file	=> $user_bg_seq_file,
		user_anchor_matrix_file => $user_anchor_matrix_file,
		user_matrix_file    => $user_matrix_file,
		bg_seq_set_name		=> $bg_seq_set_name,
		email				=> $email,

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

        var_template        => "results_seq_acsa.html"
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


