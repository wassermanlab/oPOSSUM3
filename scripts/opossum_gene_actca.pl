#!/usr/bin/env perl

=head1 NAME

opossum_gene_actca.pl

=head1 SYNOPSIS

  opossum_gene_actca.pl
      -s species
      -d results_dir
      -aid anchor_tf_id
      -dist site_distance
      -g t_gene_file [-gt t_gene_id_type]
      [-b bg_gene_file] [-bt bg_gene_id_type]
          | [-bnr bg_num_rand_genes]
      [-has_operon]
      [-biotype biotype]
      [-db tf_database]
      [-cldb tf_cluster_database]
      ([-fam tf_families] | [-famf tf_family_file])
      [-cl conservation_level]
      [-thl threshold_level] | [-th threshold]
      [-srl search_region_level] | (-up upstream_bp [-dn downstream_bp])
      [-n num_results | -zcutoff cutoff -fcutoff cutoff]
      [-sr sort_by]
      [-nh]
      [-web]
      [-j job_id]
      [-m email]
      [-utgf user_t_gene_file]
      [-ubgf user_bg_gene_file]
      [-ufamf user_tf_family_file]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

Some switches can take multiple values. To specify multiple values either
use multiple instances of the switch or a single switch followed by a comma
separated string of values (or some combination thereof).
e.g: -tax vertebrates -tax "insects, nematodes"

    -s, -species species
            The common species name for which the analysis is being
            performed, e.g.: human.

    -d, -dir dir
            Name of directory used for output results files. If the
            directory does not already exist it will be created.

    -aid, -atfid anchor_tf_id
            Anchoring TF JASPAR ID

    -dist site_distance
            Inter-binding site distance. This is the maximum distance
            between the binding sites of anchoring TF clusters and the
            binding sites of the other TF clusters of interest.

    -g, -tgf FILE
            Input file containing a list of target gene IDs with one ID per
            line. These IDs must match the ID type specified by the -gt
            option.

    -gt, -tgidt id_type
            Target gene ID type; default = 0 (Ensembl ID).

    -b, -bgf FILE
            Input file containing a list of background gene IDs with one ID
            per line. These IDs must match the ID type specified by the -bt
            option. If neither the -b nor the -bnr option is provided, ALL
            genes in the oPOSSUM database are used as background.

    -bt, -bgidt id_type
            Background gene ID type; default = 0 (Ensembl ID).

    -bnr num
            Number of genes selected randomly from the oPOSSUM database to
            use as background genes. If neither the -b nor the -bnr option
            is provided, ALL genes in the oPOSSUM database are used as
            background.

    -has_operon
            If specified, indicates this species has operons.

    -biotype type
            String containing biotype name(s). If specified, only genes of
            this type are used as background. NOT FULLY IMPLEMENTED.

    -db, -tfdb db_name
            Specifies which TF database to use; default = JASPAR_2010.

    -cldb tf_cluster_database
          Specify which TFBSCluster database to use.
          Default = TFBS_cluster.

    -fam tf_families
          Specify one or more JASPAR TF families to use.

    -famf, -tff tf_family_file
          File containing the JASPAR TF families to use.
          Format: one JASPAR TF family per line.

    -cl conservation_level
            Use this pre-defined conservation level in the analysis. The
            conservation level corresponds to average phastCons scores of
            the conserved regions as follows:
            Level    Average phastCons score
              1              0.4
              2              0.6

            Only TFBSs which fall into conserved regions with an average
            conservation score of at least the specified level are included
            in the analysis.
            Default = 1 (min. = 1)

    -th threshold
            Minimum relative TFBS position weight matrix (PWM) score to
            report in the analysis. The thresold may be spesified as a
            percentage string, e.g. '85%', or as a decimal number, e.g.
            0.85
            Default = '80%' (min. = '75%')

    -thl level
            Pre-defined threshold level for PWM score

    -srl level
            Pre-defined amount of sequence upstream/downstream of gene TSSs
            to include in analysis. Range is from 1-6. The levels correspond
            to differing amounts of upstream/downstream sequence depending
            on species. See the web-based species analysis input pages.
            Default = 3

    -up upstream_bp
            Amount of sequence upstream of genes' TSSs to include in
            analysis.

    -dn downstream_bp
            Amount of sequence downstream of genes' TSSs to include in
            analysis.

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
            over-representation, files detailing the target gene binding

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

    -utgf, -ugf FILE
            Original name of the user supplied target gene file for
            informational display purposes only.

    -ubgf, -ubf FILE
            Original name of the user supplied background gene file for
            informational display purposes only.

    -ufamf, -utff user_tf_family_file
            User supplied TF family file name.

=head1 DESCRIPTION

Take a list of gene IDs, optional background gene IDs and optional subset
of transcription factor (TF) families (TFs) specified in an input file.
Specify an anchoring TF ID and inter-binding site distance. Also optionally 
specify PWM score threshold.

Count the number of TFBS clusters located within the specified site distance 
from each of the anchor TFBS cluster found at the given PWM score threshold for
both the test and background set of genes. Perform Fisher exact test and 
z-score analysis and output these results to the output file. Optionally write 
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

#use Bio::SeqIO;

#use TFBS::DB::JASPAR5;

use OPOSSUM::Include::GeneACTCAInclude;
use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;
use OPOSSUM::Analysis::Cluster::Zscore;
use OPOSSUM::Analysis::Cluster::Fisher;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::CombinedResultSet;
use OPOSSUM::Plot::ScoreVsGC;

use Statistics::Distributions;

#use constant DEBUG          => 0;

my $help;
my $job_id;
my $results_dir;
my $web;
my $species;
my $t_gene_file;
my $t_gene_id_type;
my $bg_gene_file;
my $bg_gene_id_type;
my $bg_num_rand_genes;
my $has_operon;
my $biotype;
my $tf_db;
my $tfcl_db;
#my @collections;
#my @tax_groups;
#my $min_ic;
my @tf_families;
my $tf_family_file;
my $anchor_tf_id;
my $max_site_dist;
my $threshold;
my $threshold_level;
my $upstream_bp;
my $downstream_bp;
my $search_region_level;
my $conservation_level;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $nh;
my $email;
my $user_t_gene_file;
my $user_bg_gene_file;
my $user_tf_family_file;
GetOptions(
    'species|s=s'   => \$species,
    'dir|d=s'       => \$results_dir,
    'atfid|aid=s'   => \$anchor_tf_id,
    'dist=i'        => \$max_site_dist,
    'tgf|g=s'       => \$t_gene_file,
    'bgf|b=s'       => \$bg_gene_file,
    'cl=i'          => \$conservation_level,
    'srl=i'         => \$search_region_level,
    'tgidt|gt=s'    => \$t_gene_id_type,
    'bgidt|bt=s'    => \$bg_gene_id_type,
    'bnr=i'         => \$bg_num_rand_genes,
    'has_operon'    => \$has_operon,
    'biotype=s'     => \$biotype,
    'tfdb|db=s'     => \$tf_db,
    'cldb|cdb=s'    => \$tfcl_db,
    'fam=s'         => \@tf_families,
    'famf|tff=s'    => \$tf_family_file,
    #'co=s'          => \@collections,
    #'tax=s'         => \@tax_groups,
    #'ic=s'          => \$min_ic,
    'th=s'          => \$threshold,
    'thl=i'         => \$threshold_level,
    'up=i'          => \$upstream_bp,
    'dn=i'          => \$downstream_bp,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'sr=s'          => \$sort_by,
    'nh'            => \$nh,
    'web'           => \$web,
    'job_id|j=s'    => \$job_id,
    'm=s'           => \$email,
    'utgf|ugf=s'    => \$user_t_gene_file,
    'ubgf|ubf=s'    => \$user_bg_gene_file,
    'ufamf|utff=s'  => \$user_tf_family_file,
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

	#$abs_results_dir = $ENV{PWD} . "/$results_dir";
}

unless (-d $abs_results_dir) {
    die "Results directory $abs_results_dir does not exist\n";
}

#
# Initialize logging
#
my $log_file = get_log_filename('opossum_gene_actca', $results_dir);

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
	$job_id = 'opossum_gene_actca';
}

my $message = "";
my $ok = 1;

my %job_args = (
	-job_id => $job_id,
	-email => $email,
	-logger => $logger
);

unless ($species) {
	$job_args{-heading} = "Gene aCTCA";
	fatal("No species specified", %job_args);
}

my $heading = sprintf(
	"%s Anchored Combination TFBS Cluster Analysis", ucfirst $species
);
$job_args{-heading} = $heading;

unless ($anchor_tf_id) {
    fatal("No anchor TF ID specified", %job_args);
}

unless ($t_gene_file) {
    fatal("No target gene IDs file specified", %job_args);
}

#
# set optional parameters to default values if not provided by the user
#

unless (defined $t_gene_id_type) {
    $t_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $bg_gene_id_type) {
    $bg_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $conservation_level) {
	$conservation_level = DFLT_CONSERVATION_LEVEL;
}

#
# XXX
# NO, in this case use ALL genes in the oPOSSUM DB!
#
#if (!$bg_gene_file and !$bg_num_rand_genes) {
#    $bg_num_rand_genes = DFLT_BG_NUM_RAND_GENES;
#}

if ($bg_gene_file && $bg_num_rand_genes) {
	fatal("Either provide the background gene file OR provide the number " .
		" of background genes to use", %job_args);
}

if (!$search_region_level and !$upstream_bp and !$downstream_bp) {
	$search_region_level = DFLT_SEARCH_REGION_LEVEL;
}

if (!$threshold and !$threshold_level) {
	$threshold_level = DFLT_THRESHOLD_LEVEL;
}

#
# Read gene IDs
#
$logger->info("Reading target gene IDs file");
my $t_gene_ids = read_gene_ids_from_file($t_gene_file, %job_args);
unless ($t_gene_ids) {
    fatal("No target gene IDs read from file $t_gene_file", %job_args);
}

if (scalar @$t_gene_ids > MAX_TARGET_GENES) {
    fatal(
          "Number of target genes input exceeds maximum of " . MAX_TARGET_GENES
        . " allowed", %job_args
   );
}

my $bg_gene_ids;
if ($bg_gene_file) {
    $logger->info("Reading background gene IDs file");
	$bg_gene_ids = read_gene_ids_from_file($bg_gene_file, %job_args);
    unless ($bg_gene_ids) {
        fatal(
            "No background gene IDs read from file $bg_gene_file", %job_args
        );
    }
}

#
# get the necessary oPOSSUM adaptors
#

my $db_name = sprintf("%s_%s", OPOSSUM_DB_NAME, $species);

my $opdba = opossum_db_connect($species)
    || fatal("Could not connect to oPOSSUM database $db_name", %job_args);

my $ga = $opdba->get_GeneAdaptor
    || fatal("Could not get GeneAdaptor", %job_args);

my $oa = $opdba->get_OperonAdaptor
    || fatal("Could not get OperonAdaptor", %job_args);

my $cla = $opdba->get_ConservationLevelAdaptor
    || fatal("Could not get ConservationLevelAdaptor", %job_args);

my $crla = $opdba->get_ConservedRegionLengthAdaptor
    || fatal("Could not get ConservedRegionLengthAdaptor", %job_args);

my $aca = $opdba->get_AnalysisClusterCountsAdaptor
	|| fatal("Could not get AnalysisClusterCountsAdaptor", %job_args);

my $ctfbsa = $opdba->get_ConservedTFBSAdaptor()
    || fatal("Could not get ConservedTFBSAdaptor", %job_args);

my $thla = $opdba->get_ThresholdLevelAdaptor()
    || fatal("Could not get ThresholdLevelAdaptor", %job_args);

my $srla = $opdba->get_SearchRegionLevelAdaptor()
    || fatal("Could not get SearchRegionLevelAdaptor", %job_args);


#
# Fetch the target and background gene information for oPOSSUM analysis
#
$logger->info("Fetching target oPOSSUM gene information");
my (
    $t_gids,
    $t_included_gene_ids,
    $t_missing_gene_ids,
    $t_gid_gene_ids,
    $t_operon_first_gids,
    $t_operon_unique_gids,
	$bg_gids,
	$bg_included_gene_ids,
	$bg_missing_gene_ids,
	$bg_gid_gene_ids,
	$bg_operon_first_gids,
	$bg_operon_unique_gids
) = fetch_gene_data_for_opossum_analysis(
    $ga, $oa, $cla, $has_operon, $biotype,
	$t_gene_id_type, $t_gene_ids,
	$bg_gene_id_type, $bg_gene_ids, $bg_num_rand_genes,
	%job_args
);


#
# JASPAR / TF parameter settings
# Also, TFBSCluster settings
# 

$tf_db = JASPAR_DB_NAME if !$tf_db;

# if this is default, JASPAR has to be default as well
$tfcl_db = TFBS_CLUSTER_DB_NAME if !$tfcl_db;

if ($tfcl_db eq TFBS_CLUSTER_DB_NAME and $tf_db ne JASPAR_DB_NAME) {
	fatal("\nDefault TFBS cluster DB requires the default JASPAR DB",
		%job_args);
}

# if no collection specified, use the default CORE
#my $collections_str;
#if (@collections) {
#   $collections_str = join(',', @collections);
#	@collections = split (/\s*,\s*/, $collections_str);
#
#	unless (@collections) {
#		fatal("Error parsing collections", %job_args);
#	}
#} else {
#	$collections_str = 'CORE';
#	push @collections, 'CORE';
#}
#
#my $tax_groups_str;
#if (@tax_groups) {
#   $tax_groups_str = join(',', @tax_groups);
#   @tax_groups  = split /\s*,\s*/, $tax_groups_str;
#
#   unless (@tax_groups) {
#       fatal("Error parsing tax groups", %job_args);
#   }
#}

my $tf_families_str;
if (@tf_families) {
    $tf_families_str = join(',', @tf_families);
    @tf_families  = split(/\s*,\s*/, $tf_families_str);

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

#
# Connect to JASPAR database
#
$logger->info("Fetching matrices from JASPAR");

my $jdb = jaspar_db_connect($tf_db)
	|| fatal("Could not connect to JASPAR database $tf_db", %job_args);

$logger->info("Fetching anchoring TF matrix");

my $anchor_tf = $jdb->get_Matrix_by_ID($anchor_tf_id, 'PFM');

unless ($anchor_tf) {
    fatal("Error fetching anchoring TF profile from JASPAR DB", %job_args);
}

my %get_matrix_args = (
	-matrixtype => 'PFM'
);

$get_matrix_args{-family} = \@tf_families if @tf_families;
#$get_matrix_args{-min_ic} = $min_ic || DFLT_CORE_MIN_IC;
#$get_matrix_args{-tax_group} = \@tax_groups if @tax_groups;
#$get_matrix_args{-collection} = \@collections if @collections;
#$get_matrix_args{-collection} = 'CORE';

#$logger->info("Fetching TF matrices");
#my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

#unless ($matrix_set && $matrix_set->size > 0) {
#    fatal("Error fetching TF profiles from JASPAR DB", %job_args);
#}

#my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

#
# retrieve the TFBS cluster set
#

my $cldb = tfbs_cluster_db_connect($tfcl_db);
if (!$cldb) {
	fatal("Could not connect to $tfcl_db", %job_args);
}

my $tfca = $cldb->get_TFClusterAdaptor;
fatal("Error fetching TFClusterAdaptor") if !$tfca;

$logger->info("Fetching anchoring TFBS cluster");
my $anchor_cluster = $tfca->fetch_by_tf_id($anchor_tf_id);

my $tf_cluster_set = fetch_tf_cluster_set($cldb, %get_matrix_args);

unless ($tf_cluster_set and $tf_cluster_set->size > 0) {
	fatal("Error fetching TFBS clusters", %job_args);
}

#my $tf_ids = $tf_set->ids();
#my $tfcl_ids = $tf_cluster_set->ids();

#
# set default or custom parameter settings
# maybe this function should be changed to take %args instead
#

#if ($search_region_level) {
#if ($search_region_level and !$upstream_bp and !$downstream_bp) {
#	my $srl_hash = $srla->fetch_search_region_level_hash();
#	$upstream_bp = $srl_hash->{$search_region_level}->upstream_bp;
#	$downstream_bp = $srl_hash->{$search_region_level}->downstream_bp;
#}

#if ($upstream_bp and !$downstream_bp) {
#	$downstream_bp = 0;
#}

# cannot have only 1 of up/down specified
#if ((!$upstream_bp and $downstream_bp) or ($upstream_bp and !$downstream_bp)) {
#if (defined $upstream_bp != defined $downstream_bp) {
#	fatal("For custom analysis, both the upstream and downstream " .
#		"search range must be specified", %job_args);
#}

#
# set default or custom parameter settings
# maybe this function should be changed to take %args instead
#
# XXX
# Fixed fetch_analysis_parameters so that is also returns threshold and
# upstream/downstream bp.
# DJA 2012/03/05
# XXX
#
my $analysis_type;
($analysis_type,
    $threshold_level,
    $search_region_level,
    $threshold,
    $upstream_bp,
    $downstream_bp
) = fetch_analysis_parameters(
        $thla,
        $srla,
        $species,
        $threshold_level,
        $search_region_level,
        $threshold,
        $upstream_bp,
        $downstream_bp
);

$logger->info("Fetching conservation levels");
my $cl_hash = $cla->fetch_conservation_level_hash();
my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();

#
# Analysis Counts
# Andrew's Q - no separate custom counts routines?
#

$logger->info("Fetching target gene anchored TFBS cluster counts");
my $t_counts = $aca->fetch_anchored_counts(
   	    -gene_ids               => $t_gids,
       	-anchor_cluster         => $anchor_cluster,
       	-cluster_set            => $tf_cluster_set,
       	-distance               => $max_site_dist,
       	-conservation_level     => $conservation_level,
       	-threshold              => $threshold,
       	-upstream_bp            => $upstream_bp,
       	-downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
       	-has_operon             => $has_operon,
       	-operon_gene_ids        => $t_operon_first_gids
);

unless ($t_counts) {
    fatal("Error fetching target gene anchored TFBS cluster counts", %job_args);
}

$logger->info("Fetching background gene anchored TFBS counts");
my $bg_counts = $aca->fetch_anchored_counts(
   	    -gene_ids               => $bg_gids,
       	-anchor_cluster         => $anchor_cluster,
       	-cluster_set            => $tf_cluster_set,
   	    -distance               => $max_site_dist,
       	-conservation_level     => $conservation_level,
       	-threshold              => $threshold,
       	-upstream_bp            => $upstream_bp,
   		-downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
       	-has_operon             => $has_operon,
       	-operon_gene_ids        => $bg_operon_first_gids
);

unless ($bg_counts) {
    fatal("Error fetching background gene anchored TFBS cluster counts",
        %job_args);
}

$logger->info("Fetching target total conserved region length");
my $t_cr_length = $crla->fetch_total_length(
    -conservation_level     => $conservation_level,
    -upstream_bp            => $upstream_bp,
    -downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
    -gene_ids               => $t_gids
);

unless ($t_cr_length) {
    fatal("Error fetching target gene total conserved region length",%job_args);
}

$logger->info("Fetching background total conserved region length");
my $bg_cr_length = $crla->fetch_total_length(
    -conservation_level     => $conservation_level,
    -upstream_bp            => $upstream_bp,
    -downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
    -gene_ids               => $bg_gids
);

unless ($bg_cr_length) {
    fatal(
        "Error fetching background gene total conserved region length",
		%job_args
    );
}

# Skipping GC content calculation
=head3
my $t_cr_gc_content = 0;
#$t_cr_gc_content = fetch_cr_gc_content(
#        $t_gids, $conservation_level,
#        $upstream_bp, $downstream_bp,
#        $biotype
#);

unless (defined $t_cr_gc_content) {
    $logger->warn("Could not fetch target gene conserved region GC content");
}

my $bg_cr_gc_content = 0;
#$bg_cr_gc_content = fetch_cr_gc_content(
#        $bg_gids, $conservation_level,
#        $upstream_bp, $downstream_bp,
#        $biotype
#);

unless (defined $bg_cr_gc_content) {
    $logger->warn(
        "Could not fetch background gene conserved region GC content"
    );
}
=cut

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
    $bg_counts,
    $t_counts,
    $bg_cr_length,
    $t_cr_length,
);
fatal("Error computing z-score", %job_args) unless $zresult_set;

#
# Use new OPOSSUM::Analysis::Cluster::CombinedResultSet to combine Fisher and 
# Z-score result sets
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
my $cresults = $cresult_set->get_list(%result_params);

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
}

if ($web) {
	$logger->info("Writing HTML results");
	write_results_html();
}

if ($ok) {
    $logger->info("Writing text results");
	my $out_file = "$abs_results_dir/" . RESULTS_TEXT_FILENAME;
    write_results_text($out_file, $cresults, $tf_cluster_set, %job_args);

	if (!$nh) {
    	$logger->info("Writing TFBS cluster details");
    	write_tfbs_cluster_details();
	}
}

$job_args{-web} = $web;
$job_args{-abs_results_dir} = $abs_results_dir;
$job_args{-rel_results_dir} = $rel_results_dir;
$job_args{-user_t_file} = $user_t_gene_file;
$job_args{-user_bg_file} = $user_bg_gene_file;
$job_args{-t_num} = scalar @$t_gids;
$job_args{-bg_num} = scalar @$bg_gids;
$job_args{-tf_db} = $tf_db;
#$job_args{-collections} = $collections_str;
#$job_args{-tax_groups} = $tax_groups_str;
$job_args{-families} = $tf_families_str;
#$job_args{-min_ic} = $min_ic;
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
=head3
sub write_results_text
{
	my ($filename, $results, $tf_cluster_set, %job_args) = @_;

    return unless $results && $results->[0];

    open(FH, ">$filename")
        || fatal("Could not create analysis results file $filename", %job_args);

    $logger->info("Writing analysis results to $filename\n");

    printf FH "TFBS cluster name\tClass\tFamily\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS cluster hits\tTarget TFBS cluster nucleotide rate\tBackground TFBS cluster hits\tBackground TFBS cluster nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tfcl = $tf_cluster_set->get_tf_cluster($result->id());

        printf FH 
            "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $tfcl->name(),
            $tfcl->class() || 'NA',
            $tfcl->family() || 'NA',
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_cluster_hits() || 0,
            $result->bg_cluster_hits() || 0,
            defined $result->t_cluster_rate()
                ? sprintf("%.3f", $result->t_cluster_rate()) : 'NA',
            defined $result->bg_cluster_rate()
                ? sprintf("%.3f", $result->bg_cluster_rate()) : 'NA',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'NA',
            defined $result->fisher_p_value()
                ? sprintf("%.3g", $result->fisher_p_value()) : 'NA';
    }
    close(FH);
}

=cut

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
		#low_seq_gc          => LOW_SEQ_GC,
		#high_seq_gc         => HIGH_SEQ_GC,
        species             => $species,
        job_id              => $job_id,
        gene_id_type        => $t_gene_id_type,
        t_gids              => $t_gids,
        t_gene_ids          => $t_gene_ids,
        t_included_gene_ids => $t_included_gene_ids,
        t_missing_gene_ids  => $t_missing_gene_ids,
        bg_gids             => $bg_gids,
        bg_gene_ids         => $bg_gene_ids,
        bg_included_gene_ids=> $bg_included_gene_ids,
        bg_missing_gene_ids => $bg_missing_gene_ids,
        num_t_gids          => $t_gids ? scalar @$t_gids : 0,
        num_t_gene_ids      => $t_gene_ids ? scalar @$t_gene_ids : 0,
        num_t_included_gene_ids => $t_included_gene_ids
                                        ? scalar @$t_included_gene_ids : 0,
        num_t_missing_gene_ids  => $t_missing_gene_ids
                                        ? scalar @$t_missing_gene_ids : 0,

        tf_db               => $tf_db,
		cl_db				=> $tfcl_db,
		cl_select_criteria	=> $tfcl_select_criteria,
        tf_cluster_set      => $tf_cluster_set,
		#collection          => \@collections,
		#t_cr_gc_content     => $t_cr_gc_content,
		#bg_cr_gc_content    => $bg_cr_gc_content,
        anchor_tf_name      => $anchor_tf->name,
		anchor_cluster		=> $anchor_cluster,
        max_site_distance   => $max_site_dist,
		#tax_groups          => \@tax_groups,
		#min_ic              => $min_ic,
		tf_families			=> \@tf_families,
        conservation_level  => $conservation_level,
        min_conservation    => $min_conservation,
        threshold           => $threshold,
        upstream_bp         => $upstream_bp,
        downstream_bp       => $downstream_bp,
        results             => $cresults,
        rel_results_dir     => $rel_results_dir,
        result_type         => $result_type,
        num_display_results => $num_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $sort_by,
        warn_zero_bg_hits   => $warn_zero_bg_hits,
        results_file        => RESULTS_TEXT_FILENAME,
        zscore_plot_file    => ZSCORE_PLOT_FILENAME,
        fisher_plot_file    => FISHER_PLOT_FILENAME,
        message             => $message,
		user_tf_family_file	=> $user_tf_family_file,
		user_t_gene_file	=> $user_t_gene_file,
		user_bg_gene_file	=> $user_bg_gene_file,
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

        var_template        => "results_gene_actca.html"
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
# Write the details of the putative TFBS clusters for each TFCluster/gene. 
# Create a html and text file for each TFCluster.
#
sub write_tfbs_cluster_details
{
    my $anchor_cluster_id     = $anchor_cluster->id();
	my $anchor_cluster_tf_ids = $anchor_cluster->tf_ids();

	# if results are truncated, show the relevant tf cluster ids only
    my $tfcl_ids;
	#my $t_gids = $t_counts->gene_ids;
	foreach my $result (@$cresults) {
		push @$tfcl_ids, $result->id;
	}

    foreach my $tfcl_id (@$tfcl_ids) 
	{
        my $tfcl = $tf_cluster_set->get_tf_cluster($tfcl_id);
		#my $tfcl_id = $tfcl->id();
		my $tfcl_tf_ids = $tfcl->tf_ids();

        my $text_filename = "$abs_results_dir/c$tfcl_id.txt";
        my $html_filename = "$abs_results_dir/c$tfcl_id.html";

        # XXX
        # Fetch TFBSs for this TF and all genes and pass to routines
        # below.
        #
        # NOTE: for operon genes only fetch for first gene in operon.
        #
        my @tfcl_genes;
        my %gid_sitepairs;
        my %gid_inc;
        my %operon_genes;
        foreach my $t_gid (@$t_gids) {
            next unless $t_counts->gene_cluster_count($t_gid, $tfcl_id);

            my $gid;
            if ($has_operon) {
                my $fgid = $t_operon_first_gids->{$t_gid};
                if ($fgid) {
                    $gid = $fgid;
                } else {
                    $gid = $t_gid;
                }

                my $gene = $ga->fetch_by_id($t_gid);
                push @{$operon_genes{$gid}}, $gene;
            } else {
                $gid = $t_gid;
            }

            #
            # For operon genes, this gene may have already been included via
            # another gene in the same operon
            #
            next if $gid_inc{$gid};

            $gid_inc{$gid} = 1;

			my $anchor_cluster_tfbss = $ctfbsa->fetch_by_gene(
				-gene_id			=> $gid,
				-tf_ids				=> $anchor_cluster_tf_ids,
				-conservation_level	=> $conservation_level,
				-threshold			=> $threshold,
				-upstream_bp		=> $upstream_bp,
				-downstream_bp		=> $species eq 'yeast' ? 10000 : $downstream_bp
			);

			next if !$anchor_cluster_tfbss;

            my $cluster_tfbss = $ctfbsa->fetch_by_gene(
				-gene_id			=> $gid,
				-tf_ids				=> $tfcl_tf_ids,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $species eq 'yeast' ? 10000 : $downstream_bp
            );

			next if !$cluster_tfbss;

			my $sitepairs = find_cluster_sitepairs(
				$anchor_cluster_tfbss, $anchor_cluster_id, $cluster_tfbss,
                $tfcl_id, $max_site_dist
			);

			if ($sitepairs) {
				my $gene = $ga->fetch_by_id($gid);
				$gene->fetch_promoters();

				push @tfcl_genes, $gene;

				$gid_sitepairs{$gid} = $sitepairs;
			}
        }

        write_tfbs_cluster_details_text(
            $text_filename, 
			$tfcl, $anchor_cluster,
			\@tfcl_genes, \%gid_sitepairs, $t_gid_gene_ids, $t_gene_id_type,
			$has_operon, \%operon_genes,
			%job_args
        );

        write_tfbs_cluster_details_html(
            $html_filename, $rel_results_dir, $species,
			$tfcl, $anchor_cluster, $anchor_tf,
			\@tfcl_genes, \%gid_sitepairs, $t_gid_gene_ids, $t_gene_id_type,
			$has_operon, \%operon_genes,
			%job_args
        );
    }
}

