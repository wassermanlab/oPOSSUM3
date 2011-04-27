#!/usr/local/bin/perl -w

=head1 NAME

opossum_gene_acsa.pl

=head1 SYNOPSIS

  opossum_gene_acsa.pl
        -j job_id
        -d results_dir
        -s species
        -g t_gene_file
        -id tf_id
        -cl conservation_level
        -dist site_distance
        [-b bg_gene_file]
        [-gt t_gene_id_type]
        [-bt bg_gene_id_type]
        [-bnr num_rand_bg_genes]
        [-has_operon]
        [-biotype biotype]
        [-db tf_database]
        [-co collection]
        [-tf tf_file]
        [-tax tax_groups]
        [-ic min_ic]
        [-th threshold]
        [-up upstream_bp]
        [-dn downstream_bp]
        [-n num_results | -zcutoff cutoff -fcutoff cutoff]
        [-sr sort_by]
        [-m email]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -j job ID            = The oPOSSUM job ID. This is also the same as the
                          temporary subdirectory created in the results
                          directory.
   -d results_dir       = Name of directory used for input gene ID and TF ID
                          files and output results files
   -s species           = Species
   -id tf_id            = Anchor TF ID
   -cl conservation_level = Conservation level at which to perform analysis
   -dist site_distance  = Maximum inter-binding site distance
   -g t_gene_file       = File containing list of target gene IDs
   -b bg_gene_file      = File containing background gene IDs
   -gt t_gene_id_type   = Target gene ID type
   -bt bg ID type       = Background gene ID type
   -bnr num rand bg genes
                        = Number of random background genes
   -has_operon          = Boolean indicating whether this species has
                          operons
   -biotype             = String containing biotype name(s) of genes to use
                          as background
   -db tf_database      = Specify which TF database to use
                         (default = JASPAR_2010)
   -co collection       = Specify which JASPAR collection to use
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
   -m email             = E-mail address of user

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

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use oPossumWebOpt;
use oPossumGeneWebOpt;
use oPossumACSAWebOpt;

use lib OPOSSUM_LIB_PATH;

use Getopt::Long;
use Pod::Usage;
use File::Temp;
use Carp;
#use CGI::Carp qw(carpout);
use Template;
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
use constant BG_COLOR_CLASS => 'bgc_gene_acsa';

my $job_id;
my $results_dir;
my $species;
my $t_gene_file;
my $anchor_tf_id;
my $conservation_level;
my $max_site_dist;
my $t_gene_id_type;
my $bg_gene_file;
my $bg_gene_id_type;
my $bg_num_rand_genes;
my $has_operon;
my $biotype;
my $tf_db;
my $collection;
my $tax_groups_str;
my $tf_file;
my $min_ic;
my $threshold;
my $upstream_bp;
my $downstream_bp;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $email;
GetOptions(
    'j=s'           => \$job_id,
    'd=s'           => \$results_dir,
    's=s'           => \$species,
    'g=s'           => \$t_gene_file,
    'id=s'          => \$anchor_tf_id,
    'b=s'           => \$bg_gene_file,
    'cl=i'          => \$conservation_level,
    'dist=i'        => \$max_site_dist,
    'gt=s'          => \$t_gene_id_type,
    'bt=s'          => \$bg_gene_id_type,
    'bnr=i'         => \$bg_num_rand_genes,
    'has_operon'    => \$has_operon,
    'biotype=s'     => \$biotype,
    'db=s'          => \$tf_db,
    'co=s'          => \$collection,
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
    'm=s'           => \$email
);

#my $exclude_single_hits = $no_exclude_single_hits ? 0 : 1;

die "No results directory specified\n" if !$results_dir;

# Create relative results dir name from abs results dir
my $rel_results_dir = $results_dir;

# Remove absolute path
$rel_results_dir =~ s/.*\///; 

# Add relative path
$rel_results_dir = REL_HTDOCS_RESULTS_PATH . "/$rel_results_dir";

#
# Initialize logging
#
my $log_file = "$results_dir/opossum_gene_acsa.log";

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

$tf_db = JASPAR_DB_NAME if !$tf_db;

$collection = 'CORE' if !$collection;

$threshold = DFLT_THRESHOLD . "%" if !$threshold;

unless ($species) {
    fatal("No species specified");
}

unless ($anchor_tf_id) {
    fatal("No anchor TF ID specified");
}

unless ($t_gene_file) {
    fatal("No target gene IDs file specified");
}

unless (defined $t_gene_id_type) {
    $t_gene_id_type = DFLT_GENE_ID_TYPE;
}

unless (defined $bg_gene_id_type) {
    $bg_gene_id_type = DFLT_GENE_ID_TYPE;
}

$logger->info("Reading target gene IDs file");
my $t_gene_ids = read_ids_file($t_gene_file);
unless ($t_gene_ids) {
    fatal("No target gene IDs read from file $t_gene_file");
}

my $bg_gene_ids;
if ($bg_gene_file) {
    $logger->info("Reading background gene IDs file");
    $bg_gene_ids = read_ids_file($bg_gene_file);
    unless ($bg_gene_ids) {
        fatal(
            "No background gene IDs read from file $t_gene_file"
        );
    }
}

my $db_name = sprintf("%s_%s", OPOSSUM_DB_NAME, $species);

my $opdba = OPOSSUM::DBSQL::DBAdaptor->new(
    -host     => OPOSSUM_DB_HOST,
    -dbname   => $db_name,
    -user     => OPOSSUM_DB_USER,
    -password => OPOSSUM_DB_PASS
);

unless ($opdba) {
    fatal("Could not connect to oPOSSUM database $db_name");
}

my $ga = $opdba->get_GeneAdaptor
    || fatal("Could not get GeneAdaptor");

my $cla = $opdba->get_ConservationLevelAdaptor
    || fatal("Could not get ConservationLevelAdaptor");

my $crla = $opdba->get_ConservedRegionLengthAdaptor
    || fatal("Could not get ConservedRegionLengthAdaptor");

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
    $ga, $t_gene_id_type, $t_gene_ids, $has_operon
);

if (!$t_gids || !$t_gids->[0]) {
    fatal("Error fetching target gene oPOSSUM gene IDs");
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
        $ga, $bg_gene_id_type, $bg_gene_ids, $has_operon
    );

    if (!$bg_gids || !$bg_gids->[0]) {
        fatal("Error fetching background oPOSSUM gene IDs");
    }
} elsif ($bg_num_rand_genes) {
    $logger->info("Fetching random background oPOSSUM gene IDs");
    ($bg_gids,
     $bg_operon_first_gids,
     $bg_operon_unique_gids
    ) = fetch_random_opossum_gene_ids(
        $ga, $bg_num_rand_genes, $has_operon, $biotype
    );

    if (!$bg_gids || !$bg_gids->[0]) {
        fatal(
            "Error fetching random background oPOSSUM gene IDs"
        );
    }
} else {
    $logger->info("Fetching all background oPOSSUM gene IDs");
    $bg_gids = $ga->fetch_gene_ids();
}


my %get_matrix_args = (
    -collection => $collection,
    -matrixtype => 'PFM'
);

my @tax_groups;
if ($tf_file) {
    # Specific IDs take precendence over tax group and min IC
    $logger->info("Reading TF IDs file");
    my $tf_ids = read_ids_file($tf_file);

    unless ($tf_ids) {
        fatal("No TF IDs read from file $tf_file");
    }

    $get_matrix_args{-ID} = $tf_ids;
} else {
    if ($tax_groups_str) {
        @tax_groups  = split /\s*,\s*/, $tax_groups_str;

        unless (@tax_groups) {
            fatal("Error parsing tax groups");
        }

        $get_matrix_args{-tax_group} = \@tax_groups;
    }

    if (defined $min_ic) {
        $get_matrix_args{-min_ic} = $min_ic;
    }
}

#
# Connect to JASPAR database
#
my $jdb = TFBS::DB::JASPAR5->connect(
    "dbi:mysql:" . $tf_db . ":" . JASPAR_DB_HOST,
    JASPAR_DB_USER,
    JASPAR_DB_PASS
);

fatal("Could not connect to JASPAR database $tf_db") if !$jdb;

$logger->info("Fetching anchoring TF matrix");
my $anchor_matrix = $jdb->get_Matrix_by_ID($anchor_tf_id, 'PFM');
unless ($anchor_matrix) {
    fatal("Error fetching anchoring TF profile from JASPAR DB");
}

$logger->info("Fetching TF matrices");
my $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

unless ($matrix_set && $matrix_set->size > 0) {
    fatal("Error fetching TF profiles from JASPAR DB");
}

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

my $aca = $opdba->get_AnalysisCountsAdaptor();
if (!$aca) {
    return fatal("Could not get AnalysisCountsAdaptor");
}

my $ctfbsa = $opdba->get_ConservedTFBSAdaptor();
if (!$ctfbsa) {
    return fatal("Could not get ConservedTFBSAdaptor");
}

my $tf_ids = $tf_set->ids();

$logger->info("Fetching target gene anchored TFBS counts");
my $t_counts = $aca->fetch_anchored_counts(
        -gene_ids               => $t_gids,
        -anchor_tf_id           => $anchor_tf_id,
        -tf_ids                 => $tf_ids,
        -distance               => $max_site_dist,
        -conservation_level     => $conservation_level,
        -threshold              => $threshold,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp,
        -has_operon             => $has_operon,
        -operon_gene_ids        => $t_operon_first_gids
);

unless ($t_counts) {
    fatal("Error fetching target gene anchored TFBS counts");
}

$logger->info("Fetching background gene anchored TFBS counts");
my $bg_counts = $aca->fetch_anchored_counts(
        -gene_ids               => $bg_gids,
        -anchor_tf_id           => $anchor_tf_id,
        -tf_ids                 => $tf_ids,
        -distance               => $max_site_dist,
        -conservation_level     => $conservation_level,
        -threshold              => $threshold,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp,
        -has_operon             => $has_operon,
        -operon_gene_ids        => $bg_operon_first_gids
);

unless ($bg_counts) {
    fatal("Error fetching background gene anchored TFBS counts");
}

$logger->info("Fetching target total conserved region length");
my $t_cr_length = $crla->fetch_total_length(
    -conservation_level     => $conservation_level,
    -upstream_bp            => $upstream_bp,
    -downstream_bp          => $downstream_bp,
    -gene_ids               => $t_gids
);

unless ($t_cr_length) {
    fatal("Error fetching target gene total conserved region length");
}

$logger->info("Fetching background total conserved region length");
my $bg_cr_length = $crla->fetch_total_length(
    -conservation_level     => $conservation_level,
    -upstream_bp            => $upstream_bp,
    -downstream_bp          => $downstream_bp,
    -gene_ids               => $bg_gids
);

unless ($bg_cr_length) {
    fatal(
        "Error fetching background gene total conserved region length"
    );
}

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

my $fisher = OPOSSUM::Analysis::Fisher->new();
fatal("Error initializing Fisher analysis") unless $fisher;

$logger->info("Computing Fisher scores");
my $fresult_set = $fisher->calculate_Fisher_probability(
    $bg_counts,
    $t_counts
);
fatal("Error performing Fisher analysis") unless $fresult_set;

my $zscore = OPOSSUM::Analysis::Zscore->new();
fatal("Error initializing z-score analysis") unless $zscore;

$logger->info("Computing Z-scores");
my $zresult_set = $zscore->calculate_Zscore(
    $bg_counts,
    $t_counts,
    $bg_cr_length,
    $t_cr_length,
    $tf_set
);
fatal("Error computing z-score") unless $zresult_set;

$logger->info("Combining Fisher and Z-scores");
my $cresult_set = OPOSSUM::Analysis::CombinedResultSet->new(
    -fisher_result_set  => $fresult_set,
    -zscore_result_set  => $zresult_set
);
fatal("Error combining Fisher and z-score result_set")
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
my $results = $cresult_set->get_list(%result_params);

$logger->info("Writing HTML results");
write_results_html();

if ($results && $results->[0]) {
    $logger->info("Writing text results");
    write_results_text();

    $logger->info("Writing TFBS details");
    write_tfbs_details();
}

$logger->info("Finished analysis");

send_email($email) if $email;

exit;

#
# For the provided list of Ensembl / external gene IDs, fetch the associated
# oPOSSUM gene IDs from the DB. Also keep track of which provided gene IDs
# mapped (included) or did not map (missing) to oPOSSUM gene IDs.
#
# NOTE: There can be a 1-to-many mapping of oPOSSUM gene IDs to external gene
# IDs. Make sure all of these are retained in the included external gene IDs
# list.
#
# Andrew's NOTE: This function now takes an additional argument $has_operon
# This determines whether the operon information should be used when fetching
# the genes.
#
# If operons are used, then for each gene which is in an operon, store the
# mapping of the genes ID to the ID of the first gene in the operon.
#
# Also return an additional array which contains the actual gene IDs to use
# in fetching TFBSs, TFBS counts, conserved regions lengths etc. This list
# will contain:
# - for a gene which is NOT in an operon, it's own gene ID
# - for a gene which is in an operon, the gene ID of the first gene in the
#   operon
#
sub fetch_opossum_gene_ids
{
    my ($ga, $id_type, $ids, $has_operon) = @_;

    my @mapped_ids;
    my @unmapped_ids;
    my %gene_id_ext_ids;

    #
    # This will be set to a mapping of operon gene IDs to the first gene ID
    # in the operon.
    #
    my %operon_first_gene_id;

    #
    # This will be set to a list of gene IDs for both genes which do not
    # belong to an operon and genes which are the first gene in an operon 
    # that one or more of the input genes belongs to.
    #
    my @operon_unique_gene_ids;

    #
    # These are never used
    #
    #my %ext_id_gene_ids;
    #my %gene_id_ensembl_ids; 

    my $gene_ids;
    if (defined $id_type && $id_type == DFLT_GENE_ID_TYPE) {
        #
        # Ensembl IDs provided.
        # 
        $gene_ids = $ga->fetch_gene_ids_by_ensembl_ids(
            $ids,
            -mapped_ensembl_ids     => \@mapped_ids,
            -unmapped_ensembl_ids   => \@unmapped_ids,
            #-gene_id_ensembl_ids    => \%gene_id_ensembl_ids
            #-ensembl_id_gene_ids    => \%ext_id_gene_ids
        );
        #%gene_id_ext_ids = %gene_id_ensembl_ids;
    } else { 
        #
        # Fetch by external gene IDs
        # 
        $gene_ids = $ga->fetch_gene_ids_by_external_ids(
            $ids,
            $id_type,
            -mapped_ext_ids         => \@mapped_ids,
            -unmapped_ext_ids       => \@unmapped_ids,
            -gene_id_ext_ids        => \%gene_id_ext_ids
            #-ext_id_gene_ids        => \%ext_id_gene_ids
            #-gene_id_ensembl_ids    => \%gene_id_ensembl_ids
        );
    } 

    if ($has_operon) {
        my $oa = $opdba->get_OperonAdaptor();
        unless ($oa) {
            fatal("Error getting OperonAdaptor");
        }

        my %fgene_id_inc;
        foreach my $gene_id (@$gene_ids) {
            # first check to see if this gene belongs to an operon
            # it is possible to have multiple genes belonging to the same
            # operon to be included in the genes to be fetched.
            # map each gene id in the operon to the first gene id
            my $operon = $oa->fetch_by_gene_id($gene_id);

            if ($operon) {
                #print STDERR "Operon " . $operon->id . " for gene $gene_id\n";
                my $fgene = $operon->fetch_first_gene();

                $operon_first_gene_id{$gene_id} = $fgene->id;

                unless ($fgene_id_inc{$fgene->id}) {
                    push @operon_unique_gene_ids, $fgene->id;
                }
            } else {
                push @operon_unique_gene_ids, $gene_id;
            }
        }
    } 

    return (
        $gene_ids,
        @mapped_ids                 ? \@mapped_ids : undef,
        @unmapped_ids               ? \@unmapped_ids : undef,
        #%gene_id_ensembl_ids        ? \%gene_id_ensembl_ids : undef,
        %gene_id_ext_ids            ? \%gene_id_ext_ids : undef,
        #%ext_id_gene_ids            ? \%ext_id_gene_ids : undef,
        #%operon_first_gene_to_gene_ids
        #                            ? \%operon_first_gene_to_gene_ids : undef
        %operon_first_gene_id       ? \%operon_first_gene_id : undef,
        @operon_unique_gene_ids     ? \@operon_unique_gene_ids : undef
    );
}

sub fetch_random_opossum_gene_ids
{
    my ($ga, $num_genes, $has_operon, $biotype) = @_;

    return if !$num_genes;

    $biotype = _parse_biotype($biotype) if $biotype;

    my $where = _biotype_where_clause($biotype);

    my $total_num_genes = $ga->fetch_gene_count($where);

    $num_genes = $total_num_genes if $num_genes > $total_num_genes;

    my %fetch_args = (-num_genes => $num_genes);

    if ($biotype) {
        if (ref $biotype eq 'ARRAY') {
            $fetch_args{-biotypes} = $biotype;
        } else {
            $fetch_args{-biotype} = $biotype;
        }
    }

    my $rand_gids = $ga->fetch_random_gene_ids(%fetch_args);

    #if ($has_operon) {
    #    my %incl_operon;
    #    my %incl_gid;
    #    my @passed_gids;

    #    my $oa = $opdba->get_OperonAdaptor();
    #    unless ($oa) {
    #        $self->_error("Error getting OperonAdaptor");
    #        return;
    #    }

    #    # so that we can skip them if they turn up again
    #    foreach my $rand_gid (@$rand_gids) {
    #        $incl_gid{$rand_gid} = 1;
    #    }

    #    my $gid_opid_map = $oa->fetch_gene_ids_to_operon_ids_hash($rand_gids);

    #    foreach my $gid (@$rand_gids) {
    #        my $operon = $$gid_opid_map{$gid};
    #        if (!$operon) {
    #            push @passed_gids, $gid;
    #            next;
    #        }
    #        if (!$incl_operon{$operon->id}) {
    #            $incl_operon{$operon->id} = $gid;
    #            push @passed_gids, $gid;
    #        }
    #    }

    #    # just in case num_genes > scalar(@$rand_gids), don't use num_genes
    #    while (scalar(@passed_gids) < scalar(@$rand_gids)) {
    #        my $new_gid;
    #        if (@biotypes) {
    #            $new_gid = $ga->fetch_random_gene_ids(
    #                -num_genes => 1,
    #                -biotypes   => \@biotypes
    #            );
    #        } else {
    #            $new_gid = $ga->fetch_random_gene_ids(
    #                -num_genes => 1,
    #                -biotype   => $biotype
    #            );
    #        }
    #        next if $incl_gid{$new_gid};
    #        my $op = $oa->fetch_by_gene_id($$new_gid[0]);
    #        next if $incl_operon{$op->id};
    #        push @passed_gids, $$new_gid[0];
    #    }

    #    $rand_gids = \@passed_gids;
    #}

    my %operon_first_gene_id;
    my @operon_unique_gene_ids;

    if ($has_operon) {
        my $oa = $opdba->get_OperonAdaptor();
        unless ($oa) {
            fatal("Error getting OperonAdaptor");
        }

        # if dflt_gene_id_type, # genes = 1
        my %fgene_id_inc;
        foreach my $gid (@$rand_gids) {
            # first check to see if this gene belongs to an operon
            # it is possible to have multiple genes belonging to the same operon
            # to be included in the genes to be fetched.
            # thus need to keep an array of gene ids
            # map each non-first gene id in the operon to the first gene id
            my $operon = $oa->fetch_by_gene_id($gid);
            if ($operon) {
                my $fgene = $operon->fetch_first_gene();
                #if ($$rand_gids[$i] != $fgene->id) {
                    $operon_first_gene_id{$gid} = $fgene->id;
                #}

                unless ($fgene_id_inc{$fgene->id}) {
                    push @operon_unique_gene_ids, $fgene->id;
                }
            } else {
                push @operon_unique_gene_ids, $gid;
            }
        }
    }

    return (
        $rand_gids,
        %operon_first_gene_id   ? \%operon_first_gene_id : undef,
        @operon_unique_gene_ids ? \@operon_unique_gene_ids : undef
    );
} 

#
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
    my $filename = "$results_dir/" . RESULTS_TEXT_FILENAME;

    return unless $results && $results->[0];

    open(FH, ">$filename")
        || fatal("Could not create analysis results file $filename");

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
# Ouput combined z-score/Fisher results as HTML
#
sub write_results_html
{
    my $warn_zero_bg_gene_hits = 0;
    foreach my $result (@$results) {
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

    my $heading = sprintf(
        "%s Anchored Combination Site Analysis", ucfirst $species
    );

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
        species             => $species,
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
        collection          => $collection,
        tf_set              => $tf_set,
        t_cr_gc_content     => $t_cr_gc_content,
        bg_cr_gc_content    => $bg_cr_gc_content,

        anchor_matrix       => $anchor_matrix,
        tax_groups          => \@tax_groups,
        min_ic              => $min_ic,
        conservation_level  => $conservation_level,
        min_conservation    => $min_conservation,
        threshold           => $threshold,
        upstream_bp         => $upstream_bp,
        downstream_bp       => $downstream_bp,
        results             => $results,
        rel_results_dir     => $rel_results_dir,
        result_type         => $result_type,
        num_display_results => $num_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $sort_by,
        warn_zero_bg_gene_hits  => $warn_zero_bg_gene_hits,
        results_file        => RESULTS_TEXT_FILENAME,

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

        var_template        => "results_gene_acsa.html"
    };

    my $output = process_template('master.html', $vars);

    my $html_filename = "$results_dir/" . RESULTS_HTDOCS_FILENAME;

    open(OUT, ">$html_filename")
        || fatal("Could not create HTML results file $html_filename");

    print OUT $output;

    close(OUT);

    $logger->info("Wrote HTML formatted results to $html_filename");

    return $html_filename;
}

#
# Write the details of the putative TFBSs for each TF/gene. Create a
# text file for each TF.
#
sub write_tfbs_details
{
    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $tf = $tf_set->get_tf($tf_id);

        my $tf_name = $tf->name();

        my $text_filename = "$results_dir/$tf_id.txt";
        my $html_filename = "$results_dir/$tf_id.html";

        # XXX
        # Fetch TFBSs for this TF and all genes and pass to routines
        # below.
        #
        # NOTE: for operon genes only fetch for first gene in operon.
        #
        my @tf_genes;
        my %gid_sitepairs;
        my %gid_inc;
        my %operon_genes;
        foreach my $t_gid (@$t_gids) {
            next unless $t_counts->gene_tfbs_count($t_gid, $tf_id);

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

            my $sitepairs = $ctfbsa->fetch_anchored_by_tf_id(
                -tf_id              => $tf_id,
                -anchor_tf_id       => $anchor_tf_id,
                -gene_id            => $gid,
                -distance           => $max_site_dist,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );

            if ($sitepairs) {
                #
                # This should always be true, since we previously checked
                # the count was not zero.
                #
                my $gene = $ga->fetch_by_id($gid);
                $gene->fetch_promoters();

                push @tf_genes, $gene;

                $gid_sitepairs{$gid} = $sitepairs;
            }
        }

        write_tfbs_details_text(
            $text_filename, $tf, \@tf_genes, \%gid_sitepairs, \%operon_genes
        );

        write_tfbs_details_html(
            $html_filename, $tf, \@tf_genes, \%gid_sitepairs, \%operon_genes
        );
    }
}

#
# Write the details of the putative TFBSs for each TF/gene. Create a
# text file for each TF.
#
sub write_tfbs_details_text
{
    my ($filename, $tf, $genes, $gid_sitepairs, $operon_genes) = @_;

    my $tf_name = $tf->name();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details text file $filename"
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    #
    # It's not really necessary to check this for gene based oPOSSUM as
    # all profiles come from JASPAR and are retrieved as PFMs.
    #
    my $tf_total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $tf_total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $tf_total_ic = 'N/A';
    }

    my $anchor_total_ic;
    if ($anchor_matrix->isa("TFBS::Matrix::PFM")) {
        $anchor_total_ic = sprintf("%.3f", $anchor_matrix->to_ICM->total_ic());
    } else {
        $anchor_total_ic = 'N/A';
    }

    my $text = sprintf("%s\n\n", $anchor_matrix->name());

    $text .= sprintf("JASPAR ID:\t%s\n", $anchor_matrix->ID());
    $text .= sprintf("Class:    \t%s\n", $anchor_matrix->class() || 'N/A');
    $text .= sprintf("Family:   \t%s\n", $anchor_matrix->tag('family') || 'N/A');
    $text .= sprintf("Tax group:\t%s\n",
        $anchor_matrix->tag('tax_group') || 'N/A');
    $text .= sprintf("IC:       \t%s\n", $tf_total_ic);

    $text .= "\n\n";

    $text .= sprintf("%s\n\n", $tf->name());

    $text .= sprintf("JASPAR ID:\t%s\n", $tf->ID());
    $text .= sprintf("Class:    \t%s\n", $tf->class() || 'N/A');
    $text .= sprintf("Family:   \t%s\n", $tf->tag('family') || 'N/A');
    $text .= sprintf("Tax group:\t%s\n", $tf->tag('tax_group') || 'N/A');
    $text .= sprintf("IC:       \t%s\n", $tf_total_ic);

    $text .= sprintf("\n\n%s - %s Conserved Binding Site Combinations\n\n",
        $anchor_matrix->name(), $tf->name());

    unless (defined $t_gene_id_type
            && $t_gene_id_type == DFLT_GENE_ID_TYPE)
    {
        $text .= "Gene ID(s)";
    }

    $text .= "\tEnsembl ID";

    if ($has_operon) {
        $text .= "\tOperon ID";
    }

    $text .= qq{\tChr\tStart\tEnd\tStrand\tNearest TSS\tAnchoring TF\tStart\tEnd\tRel. Start\tRel. End\tStrand\tScore\t%%Score\tSequence\tAnchored TF\tStart\tEnd\tRel. Start\tRel. End\tStrand\tScore\t%%Score\tSequence\tDistance\n};

    foreach my $gene (@$genes) {
        my $gid = $gene->id();

        my $sitepairs = $gid_sitepairs->{$gid};

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

        if ($has_operon && defined $operon_genes->{$gid}) {
            unless (defined $t_gene_id_type
                    && $t_gene_id_type == DFLT_GENE_ID_TYPE
            ) {
                my $count = 0;
                foreach my $in_gene (@$operon_genes->{$gid}) {
                    $text .= ',' if $count > 0;
                    $text .= join(',', @{$t_gid_gene_ids->{$in_gene->id()}});
                    $count++;
                }
                $text .= "\t";
            }
            
            my @ensembl_ids;
            foreach my $in_gene (@{$operon_genes->{$gid}}) {
                #push @ensembl_ids, $t_gid_ens_ids->{$in_gene->id()};
                push @ensembl_ids, $in_gene->ensembl_id();
            }
            $text .= join(',', @ensembl_ids) . "\t";
        } else {
            unless (defined $t_gene_id_type
                    && $t_gene_id_type == DFLT_GENE_ID_TYPE)
            {
                $text .= join(',', @{$t_gid_gene_ids->{$gid}}) . "\t";
            }

            $ensembl_id  = $gene->ensembl_id();
            $strand      = $gene->strand();
            $chr         = $gene->chr();
            $tss         = $gene->tss(); # this tss is not used.
                                         # overwritten later.
            $text .= "$ensembl_id\t";
        }
        
        if ($has_operon) {
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
        foreach my $sitepair (@$sitepairs) {
            my $tfbs        = $sitepair->{tf_site};
            my $anchor      = $sitepair->{anchor_site};
            my $distance    = $sitepair->{distance};

            my $site_start          = $gene_start + $tfbs->start() - 1;
            my $site_end            = $gene_start + $tfbs->end() - 1;
            my $site_seq            = $tfbs->seq();
            my $site_strand         = $tfbs->strand();
            my $site_score          = $tfbs->score();
            my $site_rel_score      = $tfbs->rel_score();

            my $anchor_start        = $gene_start + $anchor->start() - 1;
            my $anchor_end          = $gene_start + $anchor->end() - 1;
            my $anchor_seq          = $anchor->seq();
            my $anchor_strand       = $anchor->strand();
            my $anchor_score        = $anchor->score();
            my $anchor_rel_score    = $anchor->rel_score();

            my $anchor_closest_tss;
            my $site_closest_tss;
            my $min_anchor_tss_dist = 999999;
            my $min_site_tss_dist = 999999;
            foreach my $promoter (@$promoters) {
                my $tss = $promoter->tss();

                my $anchor_start_tss_dist = abs($anchor_start - $tss);
                my $anchor_end_tss_dist   = abs($anchor_end - $tss);

                if ($anchor_start_tss_dist < $min_anchor_tss_dist) {
                    $min_anchor_tss_dist = $anchor_start_tss_dist;
                    $anchor_closest_tss = $tss;
                }

                if ($anchor_end_tss_dist < $min_anchor_tss_dist) {
                    $min_anchor_tss_dist = $anchor_end_tss_dist;
                    $anchor_closest_tss = $tss;
                }

                #my $site_start_tss_dist = abs($site_start - $tss);
                #my $site_end_tss_dist   = abs($site_end - $tss);
                #
                #if ($site_start_tss_dist < $min_site_tss_dist) {
                #    $min_site_tss_dist = $site_start_tss_dist;
                #    $site_closest_tss = $tss;
                #}
                #
                #if ($site_end_tss_dist < $min_site_tss_dist) {
                #    $min_site_tss_dist = $site_end_tss_dist;
                #    $site_closest_tss = $tss;
                #}
            }
            
            my ($site_rel_start, $site_rel_end);
            my ($anchor_rel_start, $anchor_rel_end);
            if ($strand == 1) {
                $anchor_rel_start = $anchor_start - $anchor_closest_tss;
                if ($anchor_start >= $anchor_closest_tss) {
                    $anchor_rel_start++;
                }

                $anchor_rel_end = $anchor_end - $anchor_closest_tss;
                if ($anchor_end >= $anchor_closest_tss) {
                    $anchor_rel_end++;
                }

                #$site_rel_start = $site_start - $site_closest_tss;
                #if ($site_start >= $site_closest_tss) {
                #    $site_rel_start++;
                #}
                #
                #$site_rel_end = $site_end - $site_closest_tss;
                #if ($site_end >= $site_closest_tss) {
                #    $site_rel_end++;
                #}

                $site_rel_start = $site_start - $anchor_closest_tss;
                if ($site_start >= $anchor_closest_tss) {
                    $site_rel_start++;
                }

                $site_rel_end = $site_end - $anchor_closest_tss;
                if ($site_end >= $anchor_closest_tss) {
                    $site_rel_end++;
                }
            } else {
                $anchor_rel_start = $anchor_closest_tss - $anchor_start;
                if ($anchor_start <= $anchor_closest_tss) {
                    $anchor_rel_start++;
                }

                $anchor_rel_end = $anchor_closest_tss - $anchor_end;
                if ($anchor_end <= $anchor_closest_tss) {
                    $anchor_rel_end++;
                }

                # swap coords so start is more upstream than end
                ($anchor_rel_start, $anchor_rel_end)
                    = ($anchor_rel_end, $anchor_rel_start);

                #$site_rel_start = $site_closest_tss - $site_start;
                #if ($site_start <= $site_closest_tss) {
                #    $site_rel_start++;
                #}
                #
                #$site_rel_end = $site_closest_tss - $site_end;
                #if ($site_end <= $site_closest_tss) {
                #    $site_rel_end++;
                #}

                $site_rel_start = $anchor_closest_tss - $site_start;
                if ($site_start <= $anchor_closest_tss) {
                    $site_rel_start++;
                }

                $site_rel_end = $anchor_closest_tss - $site_end;
                if ($site_end <= $anchor_closest_tss) {
                    $site_rel_end++;
                }

                # swap coords so start is more upstream than end
                ($site_rel_start, $site_rel_end)
                    = ($site_rel_end, $site_rel_start);
            }

            unless ($first) {
                $text .= "\t\t\t\t\t";
                if ($has_operon) {
                    $text .= "\t";
                }

                unless (   defined $t_gene_id_type
                        && $t_gene_id_type == DFLT_GENE_ID_TYPE
                ) {
                    $text .= "\t";
                }
            }

            $text .= sprintf("\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\t%d\n",
                $anchor_closest_tss,
                $anchor_matrix->name(),
                $anchor_start,
                $anchor_end,
                $anchor_rel_start,
                $anchor_rel_end,
                $anchor_strand == 1 ? '+' : '-',
                $anchor_score,
                $anchor_rel_score * 100,
                $anchor_seq,
                $tf->name(),
                $site_start,
                $site_end,
                $site_rel_start,
                $site_rel_end,
                $site_strand == 1 ? '+' : '-',
                $site_score,
                $site_rel_score * 100,
                $site_seq,
                $distance
            );

            $first = 0;
        } # end foreach tfbs
    } # end foreach gid

    unless (open(FH, ">$filename")) {
        $logger->error("Unable to create TFBS details file $filename - $!");
        return;
    }

    print FH $text;

    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/gene. Create an
# HTML file for each TF.
#
sub write_tfbs_details_html
{
    my ($filename, $tf, $genes, $gid_sitepairs, $operon_genes) = @_;

    my $tf_name = $tf->name();
    my $tf_id   = $tf->ID();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details html file $filename"
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    my $heading = sprintf "%s Anchored Combination Site Analysis",
        ucfirst $species;

    my $title = "oPOSSUM $heading";

    my $section = sprintf(
        "%s - %s Conserved Binding Site Combinations",
        $anchor_matrix->name(), $tf_name
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
        low_matrix_ic       => LOW_MATRIX_IC,
        high_matrix_ic      => HIGH_MATRIX_IC,
        low_matrix_gc       => LOW_MATRIX_GC,
        high_matrix_gc      => HIGH_MATRIX_GC,
        low_seq_gc          => LOW_SEQ_GC,
        high_seq_gc         => HIGH_SEQ_GC,

        jaspar_url          => JASPAR_URL,

        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'N/A'
                               },

        tf_db               => $tf_db,
        anchor_matrix       => $anchor_matrix,
        tf                  => $tf,
        gene_id_type        => $t_gene_id_type,
        dflt_gene_id_type   => DFLT_GENE_ID_TYPE,
        genes               => $genes,
        gid_sitepairs       => $gid_sitepairs,
        has_operon          => $has_operon,
        operon_genes        => $operon_genes,
        rel_results_dir     => $rel_results_dir,
        tfbs_details_file   => "$tf_id.txt",
        var_template        => "tfbs_details_gene_acsa.html"
    };

    my $output = process_template('master.html', $vars);

    print FH $output;

    close(FH);
}

sub read_ids_file
{
    my ($file) = @_;

    open(FH, $file) || fatal("Could not open IDs file $file - $!");

    my @ids;
    while (my $line = <FH>) {
        chomp $line;

        push @ids, $line;
    }

    close(FH);

    return @ids ? \@ids : undef;
}

sub fetch_cr_gc_content
{
    my ($gids, $clevel, $upstream_bp, $downstream_bp, $biotype) = @_;
    
    if (!$gids or scalar(@$gids) == 0) {
        # all gene ids
        my $ga = $opdba->get_GeneAdaptor();

        my $where = _biotype_where_clause($biotype);

        $gids = $ga->fetch_gene_ids($where);
    }
    
    my $cra = $opdba->get_ConservedRegionAdaptor();
    if (!$cra) {
        fatal("Could not get ConservedRegionAdaptor");
    }
    
    my $sum_gc = 0;
    my $sum_length = 0;
    foreach my $gid (@$gids) {
        my $gc_content = $cra->fetch_gc_content_by_upstream_downstream(
            $gid, $clevel, $upstream_bp, $downstream_bp
        );

        my $cr_length = $cra->fetch_length_by_upstream_downstream(
            $gid, $clevel, $upstream_bp, $downstream_bp
        );
        
        $sum_gc += $gc_content * $cr_length;
        $sum_length += $cr_length;
    }
    
    #print STDERR "fetch_cr_gc_content: # gids = " . scalar(@$gids) . "\n";
    #print STDERR "fetch_cr_gc_content: sum(gc_content) = $sum_gc_content\n";
    my $avg_gc_content = $sum_gc / $sum_length;
    $avg_gc_content = sprintf("%.2f", $avg_gc_content);
    
    return $avg_gc_content;
}

sub process_template
{
    my ($template_name, $vars) = @_;

    my $config = {
        ABSOLUTE        => 1,
        INCLUDE_PATH    => ABS_HTDOCS_TEMPLATE_PATH . "/", # or list ref
        INTERPOLATE     => 1,           # expand "$var" in plain text
        POST_CHOMP      => 1,           # cleanup whitespace
        #PRE_PROCESS     => 'header',    # prefix each template
        EVAL_PERL => 1                  # evaluate Perl code blocks
    };

    my $string   = '';
    my $template = Template->new($config);
    my $input    = ABS_HTDOCS_TEMPLATE_PATH . "/$template_name";
    $template->process($input, $vars, \$string)
        || fatal($template->error());

    return $string;
}

sub send_email
{
    my ($email) = @_;

    return if !$email;

    my $results_url = sprintf "%s%s/%s",
        WEB_SERVER_URL,
        "$rel_results_dir",
        RESULTS_HTDOCS_FILENAME;

    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg = "Your oPOSSUM $species Anchored Combination Site Analysis results are available.\n\n";
    $msg .= "They can be viewed at $results_url\n\n";
    $msg .= "Thank-you,\n";
    $msg .= "the oPOSSUM development team\n";
    $msg .= ADMIN_EMAIL . "\n";

    if (!open(SM, "|" . $cmd)) {
        $logger->error("Could not open sendmail - $!");
        return;
    }

    printf SM "To: %s\n", $email;
    printf SM "From: %s\n", ADMIN_EMAIL;
    print SM "Subject: oPOSSUM $species aCSA results\n\n";
    print SM "$msg" ;

    close(SM);
}

sub fatal
{
    my ($error) = @_;

    $error = 'Unknown error' if !$error;

    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg = "oPOSSUM $species aCSA analysis failed\n";
    $msg .= "\nError: $error\n";
    $msg .= "\nJob ID: $job_id\n";
    $msg .= "\nUser e-mail: $email\n" if $email;

    if (open(SM, "|" . $cmd)) {
        printf SM "To: %s\n", ADMIN_EMAIL;
        print SM "Subject: oPOSSUM $species aCSA fatal error\n\n";
        print SM "$msg" ;

        close(SM);
    } else {
        $logger->error("Could not open sendmail - $!") if $logger;
    }

    $logger->logdie("$error") if $logger;
}

#
# Biotype can be a single value, a comma separated string or an array ref
#
sub _biotype_where_clause
{
    my ($biotype) = @_;

    return unless $biotype;

    $biotype = _parse_biotype($biotype);

    my $where;
    if ($biotype) {
        if (ref $biotype eq 'ARRAY') {
            $where = "where biotype in ('"
                . join("','", @$biotype)
                . "')";
        } elsif ($biotype !~ /^all$/i) {
            # don't set if biotype is 'all'
            $where = "where biotype = '$biotype'";
        }
    }

    return $where;
}

#
# Parse biotype argument. Biotype may be a single value, a multi-value string
# separated by spaces or commas or an array ref. If it is a multi-value string, # convert it to an array and return an array ref. If it is a single value and
# it is equal to 'all' return as undef. otherwise return passed argument
# as is.
#
sub _parse_biotype
{
    my $biotype = shift;

    return unless $biotype;

    return $biotype if ref $biotype eq 'ARRAY';

    $biotype =~ s/^\s+//;   # remove leading spaces
    $biotype =~ s/\s+$//;   # remove trailing spaces

    if ($biotype =~ /,/) {
        my @biotypes = split /\s*,\s*/, $biotype;

        return \@biotypes;
    } elsif ($biotype =~ /\s+/) {
        my @biotypes = split /\s+/, $biotype;

        return \@biotypes;
    } elsif ($biotype =~ /^all$/i) {
        return undef;
    }

    return $biotype;
}
