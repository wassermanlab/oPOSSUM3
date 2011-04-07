#!/usr/local/bin/perl -w

=head1 NAME

opossum_gene_ssa.pl

=head1 SYNOPSIS

  opossum_gene_ssa.pl
        -s species
        -g t_gene_file
        [-cl conservation_level]
        [-b bg_gene_file]
        [-gt t_gene_id_type]
        [-bt bg_gene_id_type]
        [-bnr bg_num_rand_genes]
        [-has_operon]
        [-biotype biotype]
        [-db tf_database]
        [-co collections]
        [-tf tf_file]
        [-tax tax_groups]
        [-ic min_ic]
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
   -db tf_database      = Specify which TF database to use
                         (default = JASPAR_2010)
   -co collections      = Specify which JASPAR collections to use
                         (default = CORE)
   -tf tf_file          = File containing list of TF IDs
   -tax tax_groups      = Specify a comma separated string of tax groups
   -ic min_ic           = Specify minimum IC of TFs
   -th threshold        = Minimum relative TFBS position weight matrix
                          (PWM) score to report in the analysis. The
                          thresold may be spefifies as a percentage, i.e.
                          '80%' or a decimal, i.e. 0.8.
                          Default = 80%
                          Min. = 75%
   -up upstream         = Amount of sequence upstream of TSS to include in
                          analysis
   -dn downstream       = Amount of sequence downstream of TSS to include in
                          analysis
   -n num_results       = Number of results to display
   -zcutoff cutoff      = Z-score cutoff of results to display
   -fcutoff cutoff      = Fisher p-value cutoff of results to display
   -sr sort_by          = Sort results by this value ('zscore', 'fisher')
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
TFBS hits file.

=head1 AUTHOR

  Andrew Kwon, based on script by David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

#use lib '/home/tjkwon/OPOSSUM/oPOSSUM3/scripts/standalone';
use lib '/space/devel/oPOSSUM3/scripts/standalone';

use oPossumGeneInclude;

#use lib OPOSSUM_LIB_PATH;

use Getopt::Long;
use Pod::Usage;
use File::Temp;
use Carp;
#use CGI::Carp qw(carpout);
#use Template;
use File::Temp qw/ tempfile tempdir /;
use Log::Log4perl qw(get_logger :levels);

use Bio::SeqIO;

use TFBS::DB::JASPAR5;

use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::CombinedResultSet;

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
my $collections_str;
my $tax_groups_str;
my $tf_file;
my $min_ic;
my $threshold;
my $upstream_bp;
my $downstream_bp;
#my $conservation;
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
    'db=s'          => \$tf_db,
    'co=s'          => \$collections_str,
    'tf=s'          => \$tf_file,
    'tax=s'         => \$tax_groups_str,
    'ic=s'          => \$min_ic,
    'th=s'          => \$threshold,
    'up=s'          => \$upstream_bp,
    'dn=s'          => \$downstream_bp,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'sr=s'          => \$sort_by,
    'd=s'           => \$results_dir,
    'o=s'           => \$out_file,
    'h=s'           => \$hits_file,
    'l=s'           => \$log_file
);

#my $exclude_single_hits = $no_exclude_single_hits ? 0 : 1;


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

# cannot have only 1 of up/down specified
if ((!$upstream_bp and $downstream_bp) or ($upstream_bp and !$downstream_bp)) {
    pod2usage(
        -msg => "For custom analysis, both the upstream and downstream "
            . "search range must be specified",
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

$out_file = "$results_dir/opossum_gene_ssa.out" if !$out_file;

$log_file = "$results_dir/opossum_gene_ssa.log" if !$log_file;

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

my @collections = split /\s*,\s*/, $collections_str if $collections_str;
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

my $aca = $opdba->get_AnalysisCountsAdaptor();
if (!$aca) {
    $logger->logdie("Could not get AnalysisCountsAdaptor");
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

if ($tf_file) {
    # Specific IDs take precendence over tax group and min IC
    $logger->info("Reading TF IDs file");
    my $tf_ids = read_ids_file($tf_file, $logger);

    unless ($tf_ids) {
        $logger->logdie("No TF IDs read from file $tf_file");
    }

    $get_matrix_args{-ID} = $tf_ids;
} else {
    $get_matrix_args{-collection} = \@collections if (@collections);
    
    $get_matrix_args{-tax_group} = \@tax_groups if (@tax_groups);

    $get_matrix_args{-min_ic} = $min_ic if defined $min_ic;
}

#
# Connect to JASPAR database and retrieve the matrices
#

my $jdb = jaspar_db_connect($tf_db)
    || $logger->logdie("Could not connect to JASPAR DB");

$logger->info("Fetching TF matrices");
my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

unless ($matrix_set && $matrix_set->size > 0) {
    $logger->logdie("Error fetching TF profiles from JASPAR DB");
}

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);


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

my $t_counts;
my $bg_counts;
my $t_cr_length;
my $bg_cr_length;
if ($analysis_type eq 'default')
{    
    $logger->info("Fetching target gene TFBS counts");
    $t_counts = $aca->fetch_counts(
        -gene_ids => $t_gids,
        -conservation_level => $conservation_level,
        -threshold_level => $thresh_level,
        -search_region_level => $sr_level,
        -tf_ids => $tf_ids,
        -operon_gene_ids => $t_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching background gene TFBS counts");
    $bg_counts = $aca->fetch_counts(
        -gene_ids => $bg_gids,
        -conservation_level => $conservation_level,
        -threshold_level => $thresh_level,
        -search_region_level => $sr_level,
        -tf_ids => $tf_ids,
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
    
    $logger->info("Fetching target gene TFBS counts");
    $t_counts = $aca->fetch_custom_counts(
        -gene_ids => $t_gids,
        -tf_ids => $tf_ids,
        -conservation_level => $conservation_level,
        -threhsold => $threshold,
        -upstream_bp => $upstream_bp,
        -downstream_bp => $downstream_bp,
        -operon_gene_ids => $t_operon_first_gids,
        -has_operon => $has_operon
    );
    
    $logger->info("Fetching background gene TFBS counts");
    my $bg_counts = $aca->fetch_custom_counts(
        -gene_ids               => $bg_gids,
        -tf_ids                 => $tf_ids,
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
    $logger->logdie("Error fetching target gene TFBS counts");
}

unless ($bg_counts) {
    $logger->logdie("Error fetching background gene anchored TFBS counts");
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

my $fisher = OPOSSUM::Analysis::Fisher->new();
$logger->logdie("Error initializing Fisher analysis") unless $fisher;

$logger->info("Computing Fisher scores");
my $fresult_set = $fisher->calculate_Fisher_probability(
    $bg_counts,
    $t_counts
);
$logger->logdie("Error performing Fisher analysis") unless $fresult_set;

my $zscore = OPOSSUM::Analysis::Zscore->new();
$logger->logdie("Error initializing z-score analysis") unless $zscore;

$logger->info("Computing Z-scores");
my $zresult_set = $zscore->calculate_Zscore(
    $bg_counts,
    $t_counts,
    $bg_cr_length,
    $t_cr_length,
    $tf_set
);
$logger->logdie("Error computing z-score") unless $zresult_set;

$logger->info("Combining Fisher and Z-scores");
my $cresult_set = OPOSSUM::Analysis::CombinedResultSet->new(
    -fisher_result_set  => $fresult_set,
    -zscore_result_set  => $zresult_set
);
$logger->logdie("Error combining Fisher and z-score result_set")
    unless $cresult_set;

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

$logger->info("Getting filtered/sorted result list");
my $results = $cresult_set->get_list(%result_params);

if ($results && $results->[0]) {
    $logger->info("Writing text results");
    write_results_text($out_file, $tf_set, $results);
}

if ($hits_file) {
    write_hits_text($hits_file, $t_counts, $results);
}

$logger->info("Finished analysis");

exit;

#
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
    my ($filename, $tf_set, $results) = @_;

    return unless $results && $results->[0];

    open(FH, ">$filename")
        || $logger->logdie("Could not create analysis results file $filename");

    $logger->info("Writing analysis results to $filename\n");

    printf FH "TF Name\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tf = $tf_set->get_tf($result->id());

        printf FH 
            "%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $tf->name(),
            $tf->ID(),
            $tf->class() || 'N/A',
            $tf->tag('family') || 'N/A',
            $tf->tag('tax_group') || 'N/A',
            $tf->to_ICM->total_ic(),
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_tfbs_hits() || 0,
            $result->bg_tfbs_hits() || 0,
            defined $result->t_tfbs_rate()
                ? sprintf("%.3f", $result->t_tfbs_rate()) : 'N/A',
            defined $result->bg_tfbs_rate()
                ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'N/A',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'N/A',
            defined $result->fisher_p_value()
                ? sprintf("%.3g", $result->fisher_p_value()) : 'N/A';
    }
    close(FH);
}


#
# Write the details of the putative TFBSs for each TF/gene. Create a
# text file for each TF.
#
sub write_hits_text
{
    my ($filename, $ac, $results) = @_;
    
    open(FH, ">$filename") || $logger->logdie(
        "Could not create TFBS details text file $filename"
    );

    my $text = "TF ID";
    
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
    $text .= qq{Site Start\tSite End\tRel. Start\tRel. End\tStrand\tScore\t%%Score\tSequence\n};

    # if results are truncated, show the relevant tf ids only
    my $tf_ids;
    my $gids = $ac->gene_ids;
    foreach my $result (@$results) {
        push @$tf_ids, $result->id;
    }
    #my $tf_ids = $ac->tf_ids(); # should restrict to high scoring ones TBD

    foreach my $tf_id (@$tf_ids)
    {
        my $tf = $tf_set->get_tf($tf_id);

        my $tf_name = $tf->name();
        
        $text .= "$tf_name";

        # 
        # Fetch TFBSs for this TF and all genes and pass to routines
        # below.
        #
        # NOTE: for operon genes only fetch for first gene in operon.
        #
        my @tf_genes; # also only contain first genes
        my %gid_tfbss; # only contain first genes in operons
        my %operon_genes; # key = first gene id, val = operon genes
        if ($has_operon) {
            foreach my $gid (@$t_operon_unique_gids) {
                my $tfbss = $ctfbsa->fetch_by_gene(
                    -tf_id              => $tf_id,
                    -gene_id            => $gid,
                    -conservation_level => $conservation_level,
                    -threshold          => $threshold,
                    -upstream_bp        => $upstream_bp,
                    -downstream_bp      => $downstream_bp
                );
    
                if ($tfbss) {
                    $gid_tfbss{$gid} = $tfbss;

                    my $gene = $ga->fetch_by_id($gid);
                    push @tf_genes, $gene;
                }
            }

            foreach my $gid (@$t_gids) {
                my $fgid = $t_operon_unique_gids->{$gid};
    
                if ($fgid && $gid_tfbss{$fgid}) {
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
                    -tf_id              => $tf_id,
                    -gene_id            => $gid,
                    -conservation_level => $conservation_level,
                    -threshold          => $threshold,
                    -upstream_bp        => $upstream_bp,
                    -downstream_bp      => $downstream_bp
                );
    
                if ($tfbss) {
                    $gid_tfbss{$gid} = $tfbss;
    
                    my $gene = $ga->fetch_by_id($gid);
                    push @tf_genes, $gene;
                    
                }
            }
        }
        
        foreach my $gene (@tf_genes)
        {
            my $tfbss = $gid_tfbss{$gene->id};
            my $opgenes = $operon_genes{$gene->id};
            $text .= "\t";
            $text .= write_tfbs_gene_details($gene, $tfbss, $t_gid_gene_ids,
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
sub write_tfbs_gene_details
{
    my ($gene, $tfbss, $gid_gene_ids, $gene_id_type, $operon_genes, $has_operon) = @_;
    
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
    foreach my $tfbs (@$tfbss) {

        my $site_start          = $gene_start + $tfbs->start() - 1;
        my $site_end            = $gene_start + $tfbs->end() - 1;
        my $site_seq            = $tfbs->seq();
        my $site_strand         = $tfbs->strand();
        my $site_score          = $tfbs->score();
        my $site_rel_score      = $tfbs->rel_score();

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

