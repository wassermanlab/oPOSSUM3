#!/usr/bin/env perl

=head1 NAME

opossum_gene_ssa.pl

=head1 SYNOPSIS

  opossum_gene_ssa.pl
      -s species
      -d results_dir
      -g t_gene_file [-bt bg_gene_id_type]
      [-b bg_gene_file] [-gt t_gene_id_type]
          | [-bnr bg_num_rand_genes]
      [-has_operon]
      [-biotype biotype]
      [-db tf_database]
      [-ids tf_ids]
          | [-tff tf_file]
          | ([-co collections] [-tax tax_groups] [-ic min_ic])
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

    -tff tf_file
            File containing list of JASPAR TFBS profile matrix IDs with one
            ID per line. If specified, it takes presedence over any of the
            -tfids, -c, -tax and -ic options below.

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

    -utgf, -ugf FILE
            Original name of the user supplied target gene file for
            informational display purposes only.

    -ubgf, -ubf FILE
            Original name of the user supplied background gene file for
            informational display purposes only.

=head1 DESCRIPTION

Take a list of gene IDs, optional background gene IDs and optional subset
of transcription factors (TFs) either specified in an input file, or limited
by external (JASPAR) database name and information content or taxonomic
supergroup or all TFs in the oPOSSUM database. Also optionally specify PWM
score threshold.

Count the number of TFBSs for each TF which was found at the given
PWM score threshold for both the test and background set of genes. Perform
Fisher exact test and z-score analysis and output these results to the
output file. Optionally write details of TFBSs found in test set to detailed
TFBS hits file.

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

use OPOSSUM::Include::GeneSSAInclude;
use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::CombinedResultSet;
use OPOSSUM::Plot::ScoreVsGC;

use Statistics::Distributions;

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
my @tf_ids;
my $tf_file;
my @collections;
my @tax_groups;
my $min_ic;
my $conservation_level;
my $threshold;
my $threshold_level;
my $search_region_level;
my $upstream_bp;
my $downstream_bp;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $nh;
my $email;
my $user_t_gene_file;
my $user_bg_gene_file;
GetOptions(
    'species|s=s'   => \$species,
    'dir|d=s'       => \$results_dir,
    'tgf|g=s'       => \$t_gene_file,
    'bgf|b=s'       => \$bg_gene_file,
    'tgidt|gt=s'    => \$t_gene_id_type,
    'bgidt|bt=s'    => \$bg_gene_id_type,
    'bnr=i'         => \$bg_num_rand_genes,
    'has_operon'    => \$has_operon,
    'biotype=s'     => \$biotype,
    'tfdb|db=s'     => \$tf_db,
    'tfids|ids=s'   => \@tf_ids,
    'tff=s'         => \$tf_file,
    'co=s'          => \@collections,
    'tax=s'         => \@tax_groups,
    'ic=s'          => \$min_ic,
    'cl=i'          => \$conservation_level,
    'th=s'          => \$threshold,
    'thl=i'         => \$threshold_level,
    'srl=i'         => \$search_region_level,
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
    'help|h|?'      => \$help
);


if ($help) {
    pod2usage(
        -verbose    => 1
    );
}

unless ($species) {
    pod2usage(
        -msg        => "No species specified."
                     . " Use -help switch for detailed option list.\n",
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

unless ($t_gene_file) {
    pod2usage(
        -msg        => "No target gene IDs file specified."
                     . " Use -help switch for detailed option list.\n",
        -verbose    => 1
    );
}

#
# set optional parameters to default values if not provided by the user
#
$tf_db = JASPAR_DB_NAME unless $tf_db;

unless (defined $t_gene_id_type) {
    $t_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $bg_gene_id_type) {
    $bg_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $conservation_level) {
    $conservation_level = DFLT_CONSERVATION_LEVEL;
}

if (!$threshold && !$threshold_level) {
    $threshold_level = DFLT_THRESHOLD_LEVEL;
}

#
# XXX
# NO, in this case use ALL genes in the oPOSSUM DB!
#
#if (!$bg_gene_file and !$bg_num_rand_genes) {
#    $bg_num_rand_genes = DFLT_BG_NUM_RAND_GENES;
#}

if (!$search_region_level and !$upstream_bp and !$downstream_bp) {
	$search_region_level = DFLT_SEARCH_REGION_LEVEL;
}

unless ($job_id) {
    $job_id = 'opossum_gene_ssa';
}

my $heading = sprintf(
    "%s Single Site Analysis", ucfirst $species
);

my %job_args = (
    -job_id             => $job_id,
    -species            => $species,
    -email              => $email,
    -web                => $web,
    -user_t_file        => $user_t_gene_file,
    -user_bg_file       => $user_bg_gene_file,
    -tf_db              => $tf_db,
    -min_ic             => $min_ic,
    -conservation_level => $conservation_level,
    -threshold_level    => $threshold_level,
    -threshold          => $threshold,
    -search_region_level=> $search_region_level,
    -upstream_bp        => $upstream_bp,
    -downstream_bp      => $downstream_bp,
    -num_results        => $num_results,
    -z_cutoff           => $zscore_cutoff,
    -f_cutoff           => $fisher_cutoff,
    -sort_by            => $sort_by,
    -heading            => $heading
);

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

$job_args{-abs_results_dir} = $abs_results_dir;
$job_args{-rel_results_dir} = $rel_results_dir;

#
# Initialize logging
#
my $log_file = get_log_filename("opossum_gene_ssa", $results_dir);

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

$job_args{-logger} = $logger;

#open(ERR, ">>$log_file") || die "Could not open log file $log_file\n";

#carpout(\*ERR);

$logger->info("Starting analysis");

if ($bg_gene_file && $bg_num_rand_genes) {
    fatal("Either provide the background gene file OR provide the number " .
          "of background genes to use", %job_args);
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
        fatal("No background gene IDs read from file $bg_gene_file", %job_args);
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

my $cra = $opdba->get_ConservedRegionAdaptor
    || fatal("Could not get ConservedRegionAdaptor", %job_args);
    
my $crla = $opdba->get_ConservedRegionLengthAdaptor
    || fatal("Could not get ConservedRegionLengthAdaptor", %job_args);

my $aca = $opdba->get_AnalysisCountsAdaptor
    || fatal("Could not get AnalysisCountsAdaptor", %job_args);

my $ctfbsa = $opdba->get_ConservedTFBSAdaptor
    || fatal("Could not get ConservedTFBSAdaptor", %job_args);

my $thla = $opdba->get_ThresholdLevelAdaptor
    || fatal("Could not get ThresholdLevelAdaptor", %job_args);
    
my $srla = $opdba->get_SearchRegionLevelAdaptor
    || fatal("Could not get SearchRegionLevelAdpator", %job_args);



#
# Fetch the target and background gene information for oPOSSUM analysis
#

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

$job_args{-t_num} = scalar @$t_gids;
$job_args{-bg_num} = scalar @$bg_gids;

#
# JASPAR / TF parameter settings
# For gene-based, do not allow custom profiles to be uploaded?
#

# if no collection specified, use the default CORE
my $collections_str;
if (@collections) {
    $collections_str = join(',', @collections);
    @collections = split(/\s*,\s*/, $collections_str);
} else {
    $collections_str = 'CORE';
    push @collections, 'CORE';
}

$job_args{-collections} = $collections_str;

my $tax_groups_str;
if (@tax_groups) {
    $tax_groups_str = join(',', @tax_groups);
    @tax_groups = split(/\s*,\s*/, $tax_groups_str);

    unless (@tax_groups) {
        fatal("Error parsing tax groups", %job_args);
    }
}

$job_args{-tax_groups} = $tax_groups_str;

my $tf_ids_str;
if (@tf_ids) {
    $tf_ids_str = join(',', @tf_ids);
    @tf_ids = split(/\s*,\s*/, $tf_ids_str);
} elsif ($tf_file) {
    $logger->info("Reading TF IDs file");

    @tf_ids = @{read_tf_ids_from_file($tf_file, %job_args)};

    unless (@tf_ids) {
        fatal("No TF IDs read from file $tf_file");
    }

    $tf_ids_str = join(',', @tf_ids);
}

$job_args{-tf_ids} = $tf_ids_str;

#
# Connect to JASPAR database and retrieve the matrices
#

$logger->info("Fetching matrices from JASPAR");

my $jdb = jaspar_db_connect($tf_db)
    || fatal("Could not connect to JASPAR database $tf_db", %job_args);

my %get_matrix_args = (
    -matrixtype => 'PFM'
);

my $tf_select_criteria;
if (@tf_ids) {
    # This takes precendence over tax groups and min IC
    $tf_select_criteria = 'specific';
    $get_matrix_args{-ID} = \@tf_ids;
} else {
    $tf_select_criteria = 'min_ic';
    $get_matrix_args{-collection} = \@collections if @collections;
    $get_matrix_args{-tax_group} = \@tax_groups if @tax_groups;
    $get_matrix_args{-min_ic} = $min_ic || DFLT_CORE_MIN_IC;
}

$logger->info("Fetching TF matrices");
my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

unless ($matrix_set && $matrix_set->size > 0) {
    fatal("Error fetching TF profiles from JASPAR DB", %job_args);
}

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

my $tf_ids = $tf_set->ids();

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
# 	fatal("For custom analysis, both the upstream and downstream " .
#          "search range must be specified", %job_args);
#}

#$logger->info("Analysis Parameters - before");
#$logger->info("threshold level = $threshold_level");
#$logger->info("threshold = $threshold");
#$logger->info("search region level = $search_region_level");
#$logger->info("upstream bp = $upstream_bp");
#$logger->info("downstream bp = $downstream_bp");

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

$job_args{-threshold_level}     = $threshold_level;
$job_args{-threshold}           = $threshold;
$job_args{-search_region_level} = $search_region_level;
$job_args{-upstream_bp}         = $upstream_bp;
$job_args{-downstream_bp}       = $downstream_bp;

$logger->info("Fetching conservation levels");
my $cl_hash = $cla->fetch_conservation_level_hash();
my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();

#
# Analysis Counts
#

my $t_counts;
my $bg_counts;
my $t_cr_length;
my $bg_cr_length;

if ($analysis_type eq 'default') {    
    $logger->info("Fetching target gene TFBS counts");
    $t_counts = $aca->fetch_counts(
        -gene_ids => $t_gids,
        -conservation_level => $conservation_level,
        -threshold_level => $threshold_level,
        -search_region_level => $search_region_level,
        -tf_ids => $tf_ids,
        -operon_gene_ids => $t_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching background gene TFBS counts");
    $bg_counts = $aca->fetch_counts(
        -gene_ids => $bg_gids,
        -conservation_level => $conservation_level,
        -threshold_level => $threshold_level,
        -search_region_level => $search_region_level,
        -tf_ids => $tf_ids,
        -operon_gene_ids => $bg_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching total target conserved region length");
    $t_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -search_region_level    => $search_region_level,
        -gene_ids               => $t_gids,
        -operon_gene_ids        => $t_operon_first_gids,
        -has_operon             => $has_operon
    );
    
    $logger->info("Fetching total background conserved region length");
    $bg_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -search_region_level    => $search_region_level,
        -gene_ids               => $bg_gids,
        -operon_gene_ids        => $bg_operon_first_gids,
        -has_operon             => $has_operon
    );

} else {
    $logger->info("Fetching target gene TFBS counts");
    $t_counts = $aca->fetch_custom_counts(
        -gene_ids => $t_gids,
        -tf_ids => $tf_ids,
        -conservation_level => $conservation_level,
        -threshold => $threshold,
        -upstream_bp => $upstream_bp,
        -downstream_bp => $species eq 'yeast' ? 10000 : $downstream_bp,
        -operon_gene_ids => $t_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching background gene TFBS counts");
    $bg_counts = $aca->fetch_custom_counts(
        -gene_ids               => $bg_gids,
        -tf_ids                 => $tf_ids,
        -conservation_level     => $conservation_level,
        -threshold              => $threshold,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
        -has_operon             => $has_operon,
        -operon_gene_ids        => $bg_operon_first_gids
    );

    $logger->info("Fetching total target conserved length");
    $t_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
        -gene_ids               => $t_gids,
        -operon_gene_ids        => $t_operon_first_gids,
        -has_operon             => $has_operon
    );

    $logger->info("Fetching total background conserved length");
    $bg_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $species eq 'yeast' ? 10000 : $downstream_bp,
        -gene_ids               => $bg_gids,
        -operon_gene_ids        => $bg_operon_first_gids,
        -has_operon             => $has_operon
    );

}

unless ($t_counts) {
    fatal("Error fetching target gene TFBS counts", %job_args);
}

unless ($bg_counts) {
    fatal("Error fetching background gene TFBS counts", %job_args);
}

unless ($t_cr_length) {
    fatal(
        "Error fetching target gene total conserved region length",
        %job_args
    );
}

unless ($bg_cr_length) {
    fatal(
        "Error fetching background gene total conserved region length",
        %job_args
    );
}

# Skipping GC content calculation
=head3
my $t_cr_gc_content = fetch_cr_gc_content(
        $ga, $cra, $t_gids, $conservation_level,
        $upstream_bp, $downstream_bp,
        $biotype
);

unless (defined $t_cr_gc_content) {
    $logger->warn("Could not fetch target gene conserved region GC content");
}

my $bg_cr_gc_content = fetch_cr_gc_content(
        $ga, $cra, $bg_gids, $conservation_level,
        $upstream_bp, $downstream_bp,
        $biotype
);

unless (defined $bg_cr_gc_content) {
    $logger->warn(
        "Could not fetch background gene conserved region GC content"
    );
}
=cut

#
# score calculations
#

my $fisher = OPOSSUM::Analysis::Fisher->new();
fatal("Error initializing Fisher analysis", %job_args)
    unless $fisher;

$logger->info("Computing Fisher scores");
my $fresult_set = $fisher->calculate_Fisher_probability(
    $bg_counts,
    $t_counts
);
fatal("Error performing Fisher analysis", %job_args)
    unless $fresult_set;

my $zscore = OPOSSUM::Analysis::Zscore->new();
fatal("Error initializing z-score analysis", %job_args)
    unless $zscore;

$logger->info("Computing Z-scores");
my $zresult_set = $zscore->calculate_Zscore(
    $bg_counts,
    $t_counts,
    $bg_cr_length,
    $t_cr_length,
    $tf_set
);
fatal("Error computing z-score", %job_args)
    unless $zresult_set;

#
# Use new OPOSSUM::Analysis::CombinedResultSet to combine Fisher and
# Z-score result sets.
#
$logger->info("Combining Fisher and Z-scores");
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

    # Sort z-score from highest to lowest
    $result_params{-reverse} = 1;
}

$logger->info("Getting filtered/sorted result list");
my $cresults = $cresult_set->get_list(%result_params);

my $message = "";
my $ok = 1;
unless ($cresults) {
    $message = "No TFBSs scored above the selected Z-score/Fisher"
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

$job_args{-num_results} = scalar @$cresults;

if ($web) {
    $logger->info("Writing HTML results");
    write_results_html(); 
}

if ($ok) {
    $logger->info("Writing text results");
    my $out_file = "$abs_results_dir/" . RESULTS_TEXT_FILENAME;
    write_results_text($out_file, $cresults, $tf_set, %job_args);
    
    if (!$nh) {
        $logger->info("Writing TFBS details");
        write_tfbs_details(
#                           $t_counts, $cresults, $results_dir,
#                           $tf_set, $tf_db,
#                           $t_gid_gene_ids, $t_gene_id_type,
#                           $has_operon, $t_operon_first_gids,
#                           $species, $ga, $ctfbsa,
#                           $conservation_level, $threshold,
#                           $upstream_bp, $downstream_bp,
#                           $web, $logger
        );
    }
}

$logger->info("Plotting scores vs. profile \%GC content");

my $plotter = OPOSSUM::Plot::ScoreVsGC->new();
unless ($plotter) {
    $logger->error("Could not initialize plotting");
} else {
    my $plot_err;

    my $z_plot_file = "$abs_results_dir/" . ZSCORE_PLOT_FILENAME;
    unless(
        $plotter->plot($cresults, $tf_set, 'Z', ZSCORE_PLOT_SD_FOLD,
        $z_plot_file, \$plot_err
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

#$logger->info("Sending notification email to $email");

send_email(%job_args) if $email;

$logger->info("Finished analysis");

exit;


##########################################


#
# Ouput combined z-score/Fisher results as HTML
#
sub write_results_html
{    
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
        tf_set              => $tf_set,
        tf_select_criteria  => $tf_select_criteria,
        #t_cr_gc_content     => $t_cr_gc_content,
        #bg_cr_gc_content    => $bg_cr_gc_content,
        collections         => \@collections,
        tax_groups          => \@tax_groups,
        tf_ids              => \@tf_ids,
        min_ic              => $min_ic,
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
        warn_zero_bg_gene_hits  => $warn_zero_bg_gene_hits,
        results_file        => RESULTS_TEXT_FILENAME,
        zscore_plot_file    => ZSCORE_PLOT_FILENAME,
        fisher_plot_file    => FISHER_PLOT_FILENAME,
        message             => $message,
        user_t_gene_file    => $user_t_gene_file,
        user_bg_gene_file  => $user_bg_gene_file,
        email               => $email,

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

    my $output = process_template('master.html', $vars);

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
# For each TF/gene, write the details of the putative TFBSs out to text and
# html files.
#
sub write_tfbs_details
{
    # if results are truncated, show the relevant tf ids only
    my $tf_ids;
    my $t_gids = $t_counts->gene_ids;
    foreach my $result (@$cresults) {
        push @$tf_ids, $result->id;
    }
    
    #my $results_dir .= "/tfbs_details";
    #mkdir $results_dir;
    
    foreach my $tf_id (@$tf_ids) {
        my $tf = $tf_set->get_tf($tf_id);
        
        my $tf_name = $tf->name();
        
        my $text_filename = "$abs_results_dir/$tf_id.txt";
        my $html_filename = "$abs_results_dir/$tf_id.html";
    
        # 
        # Fetch TFBSs for this TF and all genes and pass to routines
        # below.
        #
        # NOTE: for operon genes only fetch for first gene in operon.
        #
        my @tf_genes; # also only contain first genes
        my %t_gid_tfbss; # only contain first genes in operons
        my %t_gid_inc; # check which genes have been included
        my %t_operon_genes; # key = first gene id, val = operon genes
    
        foreach my $t_gid (@$t_gids) {
            next unless $t_counts->gene_tfbs_count($t_gid, $tf_id);

            my $gid;
            if ($has_operon) {
                $gid = _fetch_operon_gene_id($t_gid, $t_operon_first_gids);
                my $gene = $ga->fetch_by_id($t_gid);
                push @{$t_operon_genes{$gid}}, $gene;
            } else {
                $gid = $t_gid;
            }
            
            #
            # For operon genes, this gene may have already been included via
            # another gene in the same operon
            #
            next if $t_gid_inc{$gid};
        
            $t_gid_inc{$gid} = 1;
            
            my $tfbss = $ctfbsa->fetch_by_gene(
                -tf_id              => $tf_id,
                -gene_id            => $gid,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $species eq 'yeast' ? 10000 : $downstream_bp
            );

            #$logger->info("Fetching TFBSS for gene = $gid; tf = $tf_id; conservation_level = $conservation_level; threshold = $threshold; upstream_bp = $upstream_bp; downstream_bp = $downstream_bp");
    
            if ($tfbss) {
                $t_gid_tfbss{$gid} = $tfbss;
                my $gene = $ga->fetch_by_id($gid);
                push @tf_genes, $gene;
            }
        }
        
        write_tfbs_details_text(
            $text_filename, $tf, \@tf_genes, \%t_gid_tfbss,
            $t_gid_gene_ids, $t_gene_id_type,
            $has_operon, \%t_operon_genes, %job_args
        );
        
        write_tfbs_details_html(
            $html_filename, $rel_results_dir, $species,
            $tf, \@tf_genes, \%t_gid_tfbss,
            $tf_db, $t_gene_id_type, $t_gid_gene_ids,
            $has_operon, \%t_operon_genes, %job_args
        ) if $web;
        
    }
}
