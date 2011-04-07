#!/usr/local/bin/perl -w

=head1 NAME

opossum_seq_tca.pl

=head1 SYNOPSIS

  opossum_seq_tca.pl
        -s species
        -g t_gene_file
        [-cl conservation_level]
        [-b bg_gene_file]
        [-gt t_gene_id_type]
        [-bt bg_gene_id_type]
        [-bnr bg_num_rand_genes]
        [-has_operon]
        [-biotype biotype]
        [-tdb tf_database]
        [-cdb tf_cluster_database]
        [-co collections]
        [-tax tax_groups]
        [-ic min_ic]
        [-ids tf_ids]
        [-th threshold]
        [-up upstream_bp]
        [-dn downstream_bp]
        [-n num_results | -zcutoff cutoff -fcutoff cutoff]
        [-sr sort_by]
        [-d results_dir]
        [-o out_file]
        [-h hits_file]
        [-l log_file]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -s species           = Species
   -cl conservation_level = Conservation level at which to perform analysis
   -g t_gene_file       = File containing list of target gene IDs
   -b bg_gene_file      = File containing background gene IDs
   -gt t_gene_id_type   = Target gene ID type
   -bt bg ID type       = Background gene ID type
   -bnr num rand genes  = Number of random background genes
   -has_operon          = Boolean indicating whether this species has
                          operons
   -biotype             = String containing biotype name(s) of genes to use
                          as background
   -tdb tf_database = Specify which TFBS database to use
                      (default = JASPAR_2010)
   -cdb tf_cluster_database = Specify which TFBSCluster database to use
                      (default = oPOSSUM_cluster)
   -co collections   = Specify JASPAR collections
   -tax tax_groups  = Specify a comma separated string of tax groups
   -ic min_ic       = Specify the minimum IC
   -ids tfcl_ids    = Specify a comma separated string of TFCluster IDs
   -th threshold    = Minimum relative TFBS position weight matrix
                      (PWM) score to report in the analysis. The
                      thresold may be spefifies as a percentage, i.e.
                      '80%' or a decimal, i.e. 0.8.
                      Default = 80%
                      Min. = 75%
   -n num_results   = Number of results to display
   -zcutoff cutoff  = Z-score cutoff of results to display
   -fcutoff cutoff  = Fisher p-value cutoff of results to display
   -sr sort_by      = Sort results by this value ('zscore', 'fisher')
   -d results_dir       = Name of directory used for input gene ID and TF ID
                          files and output results files
   -o out_file          = opossum results output file
   -h hits_file         = TFBS hits file
   -l log_file          = log file

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
TFBS hits file

=head1 AUTHOR

  Andrew Kwon, modifying script by David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

#use lib '/home/tjkwon/OPOSSUM/oPOSSUM3/scripts/standalone';
use lib '/space/devel/oPOSSUM3/scripts/standalone';

use oPossumGeneInclude;
use oPossumTCAInclude;

#use lib OPOSSUM_LIB_PATH;
use lib TFBS_CLUSTER_LIB_PATH;

use Getopt::Long;
use Pod::Usage;
use File::Temp;
use Carp;
#use CGI::Carp qw(carpout);
#use Template;
use File::Temp qw/ tempfile tempdir /;
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

use Bio::SeqIO;

use TFBS::DB::JASPAR5;
use TFBSCluster::DBSQL::DBAdaptor;

use OPOSSUM::TFSet;
use OPOSSUM::ConservedTFBS;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;
use OPOSSUM::Analysis::Cluster::Zscore;
use OPOSSUM::Analysis::Cluster::Fisher;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::CombinedResultSet;

use Statistics::Distributions;

use constant DEBUG          => 0;

my $results_dir;
my $species;
my $t_gene_file;
my $conservation_level;
my $search_region_level;
my $t_gene_id_type;
my $bg_gene_file;
my $bg_gene_id_type;
my $bg_num_rand_genes;
my $has_operon;
my $biotype;
my $tf_db;
my $cl_db;
my $collections_str;
my $tax_groups_str;
my $min_ic;
my $tf_families_str;
my $threshold;
my $upstream_bp;
my $downstream_bp;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $out_file;
my $hits_file;
my $log_file;
GetOptions(
    's=s'           => \$species,
    'g=s'           => \$t_gene_file,
    'b=s'           => \$bg_gene_file,
    'cl=i'          => \$conservation_level,
    'srl=i'         => \$search_region_level,
    'gt=s'          => \$t_gene_id_type,
    'bt=s'          => \$bg_gene_id_type,
    'bnr=i'         => \$bg_num_rand_genes,
    'has_operon'    => \$has_operon,
    'biotype=s'     => \$biotype,
    'tdb=s'         => \$tf_db,
    'cdb=s'         => \$cl_db,
    'co=s'          => \$collections_str,
    'tax=s'         => \$tax_groups_str,
    'ic=s'          => \$min_ic,
    'fam=s'         => \$tf_families_str,
    'th=s'          => \$threshold,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'sr=s'          => \$sort_by,
    'd=s'           => \$results_dir,
    'o=s'           => \$out_file,
    'h=s'           => \$hits_file,
    'l=s'           => \$log_file
);


if (!$species) {
    pod2usage(
        -msg     => "\nPlease specify the species name",
        -verbose => 1
    );
}

if (!$t_gene_file) {
    pod2usage(
        -msg     => "\nPlease specify an input gene ID file",
        -verbose => 1
    );
}

if ($bg_gene_file && $bg_num_rand_genes) {
    pod2usage(
        -msg => "\nPlease specify EITHER a background gene file"
            . " OR a number of random background genes OR"
            . " nothing (default = all genes in oPOSSUM DB)",
        -verbose => 1
    );
}

#
# set optional parameters to default values if not provided by the user
#

# Has to be JASPAR, unless you are specifying custom cluster db
# To make this really flexible, host should be user-specified as well
# TBD later
$tf_db       = JASPAR_DB_NAME if !$tf_db;

# if this is default, JASPAR has to be default as well
$cl_db = TFBS_CLUSTER_DB_NAME if !$cl_db;

if ($cl_db eq TFBS_CLUSTER_DB_NAME and $tf_db ne JASPAR_DB_NAME) {
    pod2usage(
        -msg => "\nDefault TFBS cluster DB requires the default JASPAR DB",
        -verbose => 1
    );
}

unless (defined $t_gene_id_type) {
    $t_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $bg_gene_id_type) {
    $bg_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $conservation_level) {
    $conservation_level = DFLT_CONSERVATION_LEVEL;
}

if (!$bg_gene_file and !$bg_num_rand_genes) {
    $bg_num_rand_genes = DFLT_BG_NUM_RAND_GENES;
}

$results_dir = "." if !$results_dir;

$out_file = "$results_dir/opossum_gene_tca.out" if !$out_file;

$log_file = "$results_dir/opossum_gene_tca.log" if !$log_file;

#
# Initialize logging
#

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


#$threshold      = DFLT_THRESHOLD . "%" if !$threshold;

my @collections = split /\s*,\s*/, $collections_str if $collections_str;
my @tf_families = split /\s*,\s*/, $tf_families_str if $tf_families_str;
my @tax_groups  = split /\s*,\s*/, $tax_groups_str if $tax_groups_str;


$logger->info("Reading target gene IDs file");
my $t_gene_ids = read_ids_file($t_gene_file, $logger);
unless ($t_gene_ids) {
    $logger->logdie("No target gene IDs read from file $t_gene_file");
}

my $bg_gene_ids;
if ($bg_gene_file) {
    $logger->info("Reading background gene IDs file");
    $bg_gene_ids = read_ids_file($bg_gene_file, $logger);
    unless ($bg_gene_ids) {
        $logger->logdie(
            "No background gene IDs read from file $t_gene_file"
        );
    }
}

my $opdba = opossum_db_connect($species)
    || $logger->logdie("Could not connect to oPOSSUM DB");

my $ga = $opdba->get_GeneAdaptor
    || $logger->logdie("Could not get GeneAdaptor");
    
my $oa = $opdba->get_OperonAdaptor
    || $logger->logdie("Could not get OperonAdaptor");

my $cla = $opdba->get_ConservationLevelAdaptor
    || ("Could not get ConservationLevelAdaptor");

my $cra = $opdba->get_ConservedRegionAdaptor
    || ("Could not get ConservedRegionAdaptor");
    
my $crla = $opdba->get_ConservedRegionLengthAdaptor
    || $logger->logdie("Could not get ConservedRegionLengthAdaptor");

my $aca = $opdba->get_AnalysisClusterCountsAdaptor();
if (!$aca) {
    $logger->logdie("Could not get AnalysisClusterCountsAdaptor");
}

my $ctfbsa = $opdba->get_ConservedTFBSAdaptor();
if (!$ctfbsa) {
    $logger->logdie("Could not get ConservedTFBSAdaptor");
}

my $thla = $opdba->get_ThresholdLevelAdaptor
    || $logger->logdie("Could not get ThresholdLevelAdaptor");
    
my $srla = $opdba->get_SearchRegionLevelAdaptor
    || $logger->logdie("Could not get SearchRegionLevelAdpator");


$logger->info("Fetching conservation levels");
my $cl_hash = $cla->fetch_conservation_level_hash();
my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();

$logger->info("Fetching target oPOSSUM gene IDs");
my (
    $t_gids,
    $t_included_gene_ids,
    $t_missing_gene_ids,
    #$t_gid_ens_ids,
    $t_gid_gene_ids,
    #t_gene_id_gids,
    $t_operon_first_gids,
    $t_operon_unique_gids
) = fetch_opossum_gene_ids(
    $ga, $oa, $t_gene_id_type, $t_gene_ids, $has_operon
);

if (!$t_gids || !$t_gids->[0]) {
    $logger->logdie("Error fetching target gene oPOSSUM gene IDs");
}

my $bg_gids;
my $bg_included_gene_ids;
my $bg_missing_gene_ids;
my $bg_gid_gene_ids;
my $bg_operon_first_gids;
my $bg_operon_unique_gids;
if ($bg_gene_ids) {
    $logger->info("Fetching background oPOSSUM gene IDs");
    ($bg_gids,
     $bg_included_gene_ids,
     $bg_missing_gene_ids,
     $bg_gid_gene_ids,
     $bg_operon_first_gids,
     $bg_operon_unique_gids
    ) = fetch_opossum_gene_ids(
        $ga, $oa, $bg_gene_id_type, $bg_gene_ids, $has_operon
    );

    if (!$bg_gids || !$bg_gids->[0]) {
        $logger->logdie("Error fetching background oPOSSUM gene IDs");
    }
} elsif ($bg_num_rand_genes) {
    $logger->info("Fetching random background oPOSSUM gene IDs");
    ($bg_gids,
     $bg_operon_first_gids,
     $bg_operon_unique_gids
    ) = fetch_random_opossum_gene_ids(
        $ga, $oa, $bg_num_rand_genes, $has_operon, $biotype
    );

    if (!$bg_gids || !$bg_gids->[0]) {
        $logger->logdie(
            "Error fetching random background oPOSSUM gene IDs"
        );
    }
} else {
    $logger->info("Fetching all background oPOSSUM gene IDs");
    $bg_gids = $ga->fetch_gene_ids();
}


my %get_matrix_args = (
    -matrixtype => 'PFM'
);


# AK: No longer either/or: you can have specific families AND min_ic
# The question is: what do you call tf_select_criteria now?
# Do I actually need this parameter?
# Hmm
# even if you specify the families, there's no reason it shouldn't also
# filter on min_ic and tax groups
# something that might be worthwhie porting over to gene based tca
# This should work in JASPAR5

$get_matrix_args{-family} = \@tf_families if @tf_families;
$get_matrix_args{-min_ic} = $min_ic if defined $min_ic;
$get_matrix_args{-tax_group} = \@tax_groups if @tax_groups;
$get_matrix_args{-collection} = \@collections if (@collections);

#
# Connect to JASPAR database and retrieve the matrices
#

my $jdb = jaspar_db_connect($tf_db)
    || $logger->logdie("Could not connect to JASPAR DB");

$logger->info("Fetching TF matrices");
my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

$logger->debug("Fetch matrix arguments:\n"
    . Data::Dumper::Dumper(%get_matrix_args));

unless ($matrix_set && $matrix_set->size > 0) {
    $logger->logdie("Error fetching TF profiles from JASPAR DB");
}

$logger->debug("Matrix set:\n" . Data::Dumper::Dumper($matrix_set));

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

my $cdb = tfbs_cluster_db_connect($cl_db);
if (!$cdb) {
    $logger->logdie("Could not connect to $cl_db");
}

#
# retrieve the TFBS cluster set
#
my $tf_cluster_set = fetch_tf_cluster_set($cdb, %get_matrix_args);

unless ($tf_cluster_set and $tf_cluster_set->size > 0) {
    $logger->logdie("Error fetching TFBSClusters from oPOSSUM_cluster DB");
}
$logger->debug("TFBSClusterSet:\n" . Data::Dumper::Dumper($tf_cluster_set));

#
# set default or custom parameter settings
#

my $analysis_type = 'default';

if (!$conservation_level) {
    $conservation_level = DFLT_CONSERVATION_LEVEL;
}

my $thresh_level;
my $thresh_level_hash = $thla->fetch_threshold_level_hash();
if ($threshold) {
    if ($thresh_level_hash->{$threshold}) {
        $thresh_level = $thresh_level_hash->{$threshold};
    } else {
        $analysis_type = 'custom';
    }
} else {
    $thresh_level = DFLT_THRESHOLD_LEVEL;
}

my $sr_level;
$sr_level = $search_region_level if $search_region_level;
if (!$sr_level)
{
    if (!$upstream_bp and !$downstream_bp) {
        $sr_level = DFLT_SEARCH_REGION_LEVEL;
    } else {
        my $found = 0;
        my $sr_level_hash = $srla->fetch_search_region_level_hash();
        foreach my $srl (sort keys %$sr_level_hash) {
            my $sr_up = $sr_level_hash->{$srl}->upstream_bp();
            my $sr_down = $sr_level_hash->{$srl}->downstream_bp();
            
            if ($sr_up == $upstream_bp and $sr_down == $downstream_bp) {
                $found = 1;
                $sr_level = $srl;
            }
        }
        $analysis_type = 'custom' if !$found;
    }
}

my $tf_ids = $tf_set->ids();
my $cl_ids = $tf_cluster_set->ids();

my $t_counts;
my $bg_counts;
my $t_cr_length;
my $bg_cr_length;

if ($analysis_type eq 'default')
{    
    $logger->info("Fetching target gene TFBS cluster counts");
    $t_counts = $aca->fetch_counts(
        -gene_ids => $t_gids,
        -conservation_level => $conservation_level,
        -threshold_level => $thresh_level,
        -search_region_level => $sr_level,
        -cluster_ids => $cl_ids,
        -operon_gene_ids => $t_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching background gene TFBS cluster counts");
    $bg_counts = $aca->fetch_counts(
        -gene_ids => $bg_gids,
        -conservation_level => $conservation_level,
        -threshold_level => $thresh_level,
        -search_region_level => $sr_level,
        -cluster_ids => $cl_ids,
        -operon_gene_ids => $bg_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching total target conserved region length");
    $t_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -search_region_level    => $sr_level,
        -gene_ids               => $t_gids,
        -operon_gene_ids        => $t_operon_first_gids,
        -has_operon             => $has_operon
    );
    
    $logger->info("Fetching total background conserved region length");
    $bg_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -search_region_level    => $sr_level,
        -gene_ids               => $bg_gids,
        -operon_gene_ids        => $bg_operon_first_gids,
        -has_operon             => $has_operon
    );

} else {
    
    $logger->info("Fetching target gene TFBS cluster counts");
    $t_counts = $aca->fetch_custom_counts(
        -gene_ids => $t_gids,
        -clusters => $tf_cluster_set->get_tf_cluster_list(),
        -conservation_level => $conservation_level,
        -threhsold => $threshold,
        -upstream_bp => $upstream_bp,
        -downstream_bp => $downstream_bp,
        -has_operon => $has_operon,
        -operon_gene_ids => $t_operon_first_gids
    );
    
    $logger->info("Fetching background gene TFBS cluster counts");
    my $bg_counts = $aca->fetch_custom_counts(
        -gene_ids               => $bg_gids,
        -clusters               => $tf_cluster_set->get_tf_cluster_list(),
        -conservation_level     => $conservation_level,
        -threshold              => $threshold,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp,
        -has_operon             => $has_operon,
        -operon_gene_ids        => $bg_operon_first_gids
    );

    $logger->info("Fetching total target conserved length");
    $t_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp,
        -gene_ids               => $t_gids,
        -operon_gene_ids        => $t_operon_first_gids,
        -has_operon             => $has_operon
    );

    $logger->info("Fetching total background conserved length");
    $bg_cr_length = $crla->fetch_total_length(
        -conservation_level     => $conservation_level,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp,
        -gene_ids               => $bg_gids,
        -operon_gene_ids        => $bg_operon_first_gids,
        -has_operon             => $has_operon
    );

}

unless ($t_counts) {
    $logger->logdie("Error fetching target gene TFBS cluster counts");
}

unless ($bg_counts) {
    $logger->logdie("Error fetching background gene anchored TFBS cluster counts");
}

$logger->logdie(
    "Error fetching target gene total conserved region length"
) unless $t_cr_length;

$logger->logdie(
    "Error fetching background gene total conserved region length"
) unless $bg_cr_length;

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

my $fisher = OPOSSUM::Analysis::Cluster::Fisher->new();
$logger->info("Computing Fisher scores");
my $fresults = $fisher->calculate_Fisher_probability($bg_counts, $t_counts);

$logger->debug("Fisher results:\n" . Data::Dumper::Dumper($fresults));

my $zscore = OPOSSUM::Analysis::Cluster::Zscore->new();
$logger->info("Computing z-scores");
my $zresults = $zscore->calculate_Zscore(
    $bg_counts, $t_counts, $bg_cr_length, $t_cr_length
);

$logger->debug("Z-score results:\n" . Data::Dumper::Dumper($zresults));

#
# Use new OPOSSUM::Analysis::Cluster::CombinedResultSet to combine Fisher and
# Z-score result sets.
#
my $cresults = OPOSSUM::Analysis::Cluster::CombinedResultSet->new(
    -fisher_result_set  => $fresults,
    -zscore_result_set  => $zresults
);

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

    if ($sort_by eq 'zscore') {
        # Sort z-score from highest to lowest
        $result_params{-reverse} = 1;
    }
}

my $results = $cresults->get_list(%result_params);

if ($results && $results->[0]) {
    $logger->info("Writing text results");
    write_results_text($out_file, $tf_cluster_set, $results);
}

if ($hits_file) {
    $logger->info("Writing TFBS cluster hits");
    write_hits_text($hits_file, $t_counts, $threshold);
}

$logger->info("Finished analysis");

exit;



#
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
    my ($filename, $tf_cluster_set, $results) = @_;

    return unless $results && $results->[0];

    open(FH, ">$filename")
        || $logger->logdie("Could not create analysis results file $filename");

    $logger->info("Writing analysis results to $filename\n");

    printf FH "TFBS Cluster ID\tClass\tFamily\t";
    printf FH "Target gene hits\tTarget gene non-hits\t";
    printf FH "Background gene hits\tBackground gene non-hits\t";
    printf FH "Target cluster hits\tBackground cluster hits\t";
    printf FH "Target cluster nucleotide rate\tBackground cluster nucleotide rate\t";
    printf FH "Z-score\tFisher score\n";

    foreach my $result (@$results) {
        my $cl = $tf_cluster_set->get_tf_cluster($result->id());

        printf FH 
            "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $cl->id(),
            $cl->class() || 'N/A',
            $cl->family() || 'N/A',
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_cluster_hits() || 0,
            $result->bg_cluster_hits() || 0,
            defined $result->t_cluster_rate()
                ? sprintf("%.3f", $result->t_cluster_rate()) : 'N/A',
            defined $result->bg_cluster_rate()
                ? sprintf("%.3f", $result->bg_cluster_rate()) : 'N/A',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'N/A',
            defined $result->fisher_p_value()
                ? sprintf("%.3g", $result->fisher_p_value()) : 'N/A';
    }

    close(FH);
}



#
# For each TFBS cluster/gene, write the details of the putative TFBS clusters out
# to a text file.
#
sub write_hits_text
{
    my ($filename, $ac, $threshold) = @_;

    open(FH, ">$filename") || $logger->logdie(
        "Could not create TFBS cluster details text file $filename"
    );

    $logger->info("Writing TFBS cluster details to $filename");

    my $text = "TFBS Cluster ID";
    
    unless (defined $t_gene_id_type
            && $t_gene_id_type == DFLT_GENE_ID_TYPE)
    {
        $text .= "\tGene ID(s)";
    }

    $text .= "\tEnsembl ID";

    if ($has_operon) {
        $text .= "\tOperon ID";
    }

    $text .= qq{\tChr\tProm.Start\tProm.End\tProm.Strand\tNearest TSS\t};
    $text .= qq{Site Start\tSite End\tRel. Start\tRel. End\tStrand\t%%Score\tSequence\n};

    # if results are truncated, show the relevant cluster ids only
    my $cluster_ids;
    my $gids = $ac->gene_ids;
    foreach my $result (@$results) {
        push @$cluster_ids, $result->id;
    }
    
    foreach my $cl_id (@$cl_ids)
    {
        my $cl = $tf_cluster_set->get_tf_cluster($cl_id);
        my $tf_ids = $cl->tf_ids;
        
        $text .= $cl->id . "\t";
        
        my @cl_genes;
        my %gid_cluster_sites;
        my %operon_genes; # key = first gene id, val = operon genes
        if ($has_operon) {
            foreach my $gid (@$t_operon_unique_gids) {
                my $tfbss = $ctfbsa->fetch_by_gene(
                    -tf_ids             => $tf_ids,
                    -gene_id            => $gid,
                    -conservation_level => $conservation_level,
                    -threshold          => $threshold,
                    -upstream_bp        => $upstream_bp,
                    -downstream_bp      => $downstream_bp
                );
    
                if ($tfbss) {
                    my $merged_sites = merge_cluster_sites($tfbss, $cl_id);
                    $gid_cluster_sites{$gid} = $merged_sites;

                    my $gene = $ga->fetch_by_id($gid);
                    push @cl_genes, $gene;
                }
            }

            foreach my $gid (@$t_gids) {
                my $fgid = $t_operon_unique_gids->{$gid};
    
                if ($fgid && $gid_cluster_sites{$fgid}) {
                    my $gene = $ga->fetch_by_id($gid);
                    push @{$operon_genes{$fgid}}, $gene;
                }
            }
        } else {
            #
            # Species does not have operons.
            #
            foreach my $gid (@$t_gids) {
                my $tfbss = $ctfbsa->fetch_by_gene(
                    -tf_ids             => $tf_ids,
                    -gene_id            => $gid,
                    -conservation_level => $conservation_level,
                    -threshold          => $threshold,
                    -upstream_bp        => $upstream_bp,
                    -downstream_bp      => $downstream_bp
                );
    
                if ($tfbss) {
                    my $merged_sites = merge_cluster_sites($tfbss, $cl_id);
                    $gid_cluster_sites{$gid} = $merged_sites;
    
                    my $gene = $ga->fetch_by_id($gid);
                    push @cl_genes, $gene;
                    
                }
            }
        }
        
        foreach my $gene (@cl_genes)
        {
            my $sites = $gid_cluster_sites{$gene->id};
            my $opgenes = $operon_genes{$gene->id};
            $text .= write_tfbs_cluster_gene_details($gene, $sites, $t_gid_gene_ids,
                                            $t_gene_id_type, $opgenes, $has_operon);
        }
    }
    
    print FH $text;
    
    close(FH);
}

#
# given a particular gene and the set of TFBSs for a particular TF, write out
# the locations and return the text string
#
sub write_tfbs_cluster_gene_details
{
    my ($gene, $sites, $gid_gene_ids, $gene_id_type, $operon_genes, $has_operon) = @_;
    
    my $gid = $gene->id();

    my $text;
    
    my $ensembl_id  = $gene->ensembl_id();
    my $gene_start  = $gene->start();
    my $gene_end    = $gene->end();
    my $strand      = $gene->strand();
    my $chr         = $gene->chr();
    my $tss         = $gene->tss();
    my $promoters   = $gene->promoters();

    my $prom_start;
    my $prom_end;
    if ($strand == 1) {
        $prom_start = $tss;
        $prom_end   = $gene_end;
    } else {
        $prom_start = $gene_start;
        $prom_end   = $tss;
    }

    if ($has_operon and defined $operon_genes)
    {        
        unless (defined $gene_id_type
                && $gene_id_type == DFLT_GENE_ID_TYPE
        ) {
            my $count = 0;
            foreach my $in_gene (@$operon_genes) {
                $text .= ',' if $count > 0;
                $text .= join(',', @{$gid_gene_ids->{$in_gene->id()}});
                $count++;
            }
            $text .= "\t";
        }
        
        my @ensembl_ids;
        foreach my $in_gene (@{$operon_genes}) {
            #push @ensembl_ids, $t_gid_ens_ids->{$in_gene->id()};
            push @ensembl_ids, $in_gene->ensembl_id();
        }
        $text .= join(',', @ensembl_ids);
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
        $text .= "$ensembl_id";
    }
    
    if (defined $has_operon) {
        my $symbol = "N/A";
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

    # gene_start/gene_end should refer to the first gene here
    # what about tss? which one should I use? show both?
    # actually, no. closest tss refers to the actual promoter tss
    # since the first gene tss is used, should be kept.
    
    my $first = 1;
    foreach my $site (@$sites) {

        my $site_start          = $gene_start + $site->start() - 1;
        my $site_end            = $gene_start + $site->end() - 1;
        my $site_seq            = $site->seq();
        my $site_strand         = $site->strand();
        my $site_score          = $site->score();
        my $site_rel_score      = $site->rel_score();

        my $site_closest_tss;
        my $min_site_tss_dist = 999999;
        foreach my $promoter (@$promoters) {
            my $tss = $promoter->tss();

            my $site_start_tss_dist = abs($site_start - $tss);
            my $site_end_tss_dist   = abs($site_end - $tss);
            
            if ($site_start_tss_dist < $min_site_tss_dist) {
                $min_site_tss_dist = $site_start_tss_dist;
                $site_closest_tss = $tss;
            }
            
            if ($site_end_tss_dist < $min_site_tss_dist) {
                $min_site_tss_dist = $site_end_tss_dist;
                $site_closest_tss = $tss;
            }
        }
        
        my ($site_rel_start, $site_rel_end);
        if ($strand == 1) {
            $site_rel_start = $site_start - $site_closest_tss;
            if ($site_start >= $site_closest_tss) {
                $site_rel_start++;
            }
            
            $site_rel_end = $site_end - $site_closest_tss;
            if ($site_end >= $site_closest_tss) {
                $site_rel_end++;
            }
        } else {
            $site_rel_start = $site_closest_tss - $site_start;
            if ($site_start <= $site_closest_tss) {
                $site_rel_start++;
            }

            $site_rel_end = $site_closest_tss - $site_end;
            if ($site_end <= $site_closest_tss) {
                $site_rel_end++;
            }

            # swap coords so start is more upstream than end
            ($site_rel_start, $site_rel_end) = ($site_rel_end, $site_rel_start);
        }

        unless ($first) {
            $text .= "\t\t\t\t\t";
            if ($has_operon) {
                $text .= "\t";
            }

            unless (   defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE
            ) {
                $text .= "\t";
            }
        }

        $text .= sprintf("\t%d\t%d\t%d\t%d\t%d\t%s\t%.3f\t%.1f\t%s\n",
            $site_closest_tss,
            $site_start,
            $site_end,
            $site_rel_start,
            $site_rel_end,
            $site_strand == 1 ? '+' : '-',
            $site_score,
            $site_rel_score * 100,
            $site_seq,
        );

        $first = 0;
    } # end foreach tfbs

    return $text;
}
