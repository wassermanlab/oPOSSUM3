#!/usr/local/bin/perl -w

=head1 NAME

default_analysis.pl

=head1 SYNOPSIS

  default_analysis.pl -g gene_file
        [-d opossum_db]
        [-s species] [-i gene_ID_type]
        [-bg bg_gene_file | -n num_rand_bg_genes]
        [-t tf_file | ([-ic IC] [-tax tax_groups] [-c collection])
        [-cl conservation_level]
        [-thl threshold_level]
        [-srl srch_reg_level]
        [-h tfbs_hits_file] [-nx]
        [-o out_file]


=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -d opossum_db    = Specify the name of the oPOSSUM DB to use. If not
                      specified and the species is, the database name
                      is generated from the species name, e.g.
                      oPOSSUM_2010_<species>
   -g gene_file     = File containing gene IDs to analyze. Format
                      should be one ID per line. The gene ID type
                      may be specified by the -i switch, i.e.
                      'ensembl', 'hugo' etc.
   -s species       = Species on which to base the analysis, i.e.
                      'human', 'mouse'. The gene IDs/types must
                      be valid for this species.
   -i gene_ID_type  = Type of gene ID contained in the input gene_file.
                      Should be a numeric corresponding to a value in
                      the id_type column of the external_gene_id_types
                      table. e.g. 1 equals MGI for mouse HNGC for human.
   -bg bg_gene_file = File containing gene IDs of background set of
                      genes. Format is one ID per line. Assumes the
                      same gene ID type as the target genes. If
                      neither -bg nor -n is provided all genes in the
                      oPOSSUM DB are used as background.
   -n num_rand_bg_genes = Generate a random background set of genes from
                      the oPOSSUM database of size num_rand_bg_genes.
                      If neither -bg nor -n is provided all genes in the
                      oPOSSUM DB are used as background.
   -t tf_file       = File containing TF IDs/names. The IDs/names
                      should appear one per line. Either internal
                      oPOSSUM DB IDs (numeric), external (JASPAR)
                      IDs of the form MA#### or names (anything else)
                      may be used. The -t and -c, -ic, -tax options are
                      exclusive. If none are specified, use all TFBS
                      profiles in the oPOSSUM database.
   -c collection    = Name of JASPAR collection to use.
                      default = 'CORE'
   -ic IC           = Minimum information content (specificity) of
                      TFBS profile matrices to use for analysis.
   -tax tax_groups  = List of taxonomic supergroups of TFBS profiles
                      to use for analysis, e.g.:
                        'plants', 'vertebrates', 'insects' etc...
   -cl cons_level   = Level of conservation for which to perform the
                      analysis.
                      Levels are:
                        1 = Top 30% of conserved regions with an
                        absolute minimum percent identity of 60%
                        2 = Top 20% of conserved regions with an
                        absolute minimum percent identity of 65%
                        3 = Top 10% of conserved regions with an
                        absolute minimum percent identity of 70%
                      Default = 3
   -thl threshold_level = Minimum relative TFBS position weight matrix
                      (PWM) score to include in the analysis.
                      Levels are:
                        1 = 75% of maximum possible PWM score for the
                            given TFBS
                        2 = 80% of maximum possible PWM score for the
                            given TFBS
                        3 = 85% of maximum possible PWM score for the
                            given TFBS
                      Default = 2
   -srl srch_reg_level  = Search region level at which to perform the
                      analysis. Specifies the amount of sequence
                      upstream (-) and downstream (+) of the
                      transcription start site to include in the
                      analysis.
                      Levels are:
                        1 = -10000 to +10000 bp
                        2 = -10000 to +5000 bp
                        3 = -5000 to +5000 bp
                        4 = -5000 to +2000 bp
                        5 = -2000 to +2000 bp
                        6 = -2000 to +0 bp
                      Default = 3
   -h tfbs_hits_file    = Name of file to which details of the putative
                      TFBSs are written.
   -nx              = If switch is provided, do NOT exclude results
                      in which a given TFBS was only found for a
                      single gene.
                      Default = exclude single gene hits
   -o out_file      = Name of output file to which analysis results
                      are written.
                      Default = standard output


=head1 DESCRIPTION

Perform the same analysis as through the oPOSSUM website using default
parameter selection (conservation level, transcription factor binding site
(TFBS) position weight matrix (PWM) score threshold level and search region
level.

Take a list of genes as provided in the input genes file. Optionally also
take a list of background genes or a number of randomly generated genes or
all the genes in the oPOSSUM database. Optionally take a subset of
transcription factors (TFs) either specified in an input file, or limited
by external (JASPAR) database name and information content or taxonomic
supergroup or all TFs in the oPOSSUM database. Also optionally take arguments
specifying conservation level, PWM score threshold level and search region
level (specifying amount of upstream/downstream sequence to include in
searching for TFBSs). 

Count the number of TFBSs for each TF which was found at the given
conservation level, for the given PWM score threshold level and search
region level for both the foreground and background set of genes. Perform
Fisher exact test and z-score analysis and output these results to the
output file.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

#
# If you install the oPOSSUM modules in the standard perl lib tree,
# comment out this line to use the installed lib.
#
use lib '/space/devel/oPOSSUM_2010/lib';

use Getopt::Long;
use Pod::Usage;

use TFBS::DB::JASPAR5;

use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::CombinedResultSet;
use OPOSSUM::TFSet;

use constant OPOSSUM_DB_HOST => "vm5.cmmt.ubc.ca";
use constant OPOSSUM_DB_USER => "opossum_r";
use constant OPOSSUM_DB_PASS => undef;

use constant JASPAR_DB_HOST  => "vm5.cmmt.ubc.ca";
use constant JASPAR_DB_NAME  => "JASPAR_2010";
use constant JASPAR_DB_USER  => "jaspar_r";
use constant JASPAR_DB_PASS  => undef;

use constant DFLT_CONSERVATION_LEVEL     => 3;
use constant DFLT_THRESHOLD_LEVEL        => 2;
use constant DFLT_SEARCH_REGION_LEVEL    => 3;
use constant DFLT_ID_TYPE                => 0;
use constant DFLT_NO_EXCLUDE_SINGLE_HITS => 0;
use constant DFLT_COLLECTION             => 'CORE';

my $opossum_db_name;
my $t_gene_file;
my $species_name;
my $id_type;
my $bg_gene_file;
my $num_rand_bg_genes;
my $collection;
my $tf_file;
my $min_ic;
my @tax_groups;
my $cons_level;
my $thresh_level;
my $sr_level;
my $out_file;
my $hits_file;
my $no_exclude_single_hits = DFLT_NO_EXCLUDE_SINGLE_HITS;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'g=s'   => \$t_gene_file,
    's=s'   => \$species_name,
    'i=s'   => \$id_type,
    'bg=s'  => \$bg_gene_file,
    'n=i'   => \$num_rand_bg_genes,
    'c=s'   => \$collection,
    't=s'   => \$tf_file,
    'ic=f'  => \$min_ic,
    'tax=s' => \@tax_groups,
    'cl=i'  => \$cons_level,
    'thl=i' => \$thresh_level,
    'srl=i' => \$sr_level,
    'o=s'   => \$out_file,
    'h=s'   => \$hits_file,
    'nx'    => \$no_exclude_single_hits
);

my $exclude_single_hits = $no_exclude_single_hits ? 0 : 1;

if (!$opossum_db_name && !$species_name) {
    pod2usage(
        -msg     => "\nPlease specify an oPOSSUM DB name or species name",
        -verbose => 1
    );
}

if (!$opossum_db_name && $species_name) {
    $opossum_db_name = "oPOSSUM_2010_$species_name";
}

if (!$t_gene_file) {
    pod2usage(
        -msg     => "\nPlease specify an input gene ID file",
        -verbose => 1
    );
}

if ($bg_gene_file && $num_rand_bg_genes) {
    pod2usage(
        -msg => "\nPlease specify EITHER a background gene file"
            . " OR a number of random background genes OR"
            . " nothing (default = all genes in oPOSSUM DB)",
        -verbose => 1
    );
}

if ($tf_file) {
    if ($collection || $min_ic || @tax_groups) {
        die "\nWhen specifying TFBSs, please provide only one of either a TFBS"
            . " file, a minimum information content or one or more taxonomic"
            . " supergroups\n";
    }
} else {
    if (!$collection) {
        $collection = DFLT_COLLECTION;
    }
}

#
# Multiple taxonomic supergroups allowed. Can be specified as
# -tax vertebrate -tax plant... OR -tax "vertebrate, plant"
#
if (@tax_groups) {
    @tax_groups = split(/[,;:\s]+/, join('[,;:\s]+', @tax_groups));
}

#
# Connect to oPOSSUM database
#
my $opdb = OPOSSUM::DBSQL::DBAdaptor->new(
    -host     => OPOSSUM_DB_HOST,
    -user     => OPOSSUM_DB_USER,
    -password => OPOSSUM_DB_PASS,
    -dbname   => $opossum_db_name
);

die "Could not connect to oPOSSUM database\n" if !$opdb;

#
# Connect to JASPAR database
#
my $jdb = TFBS::DB::JASPAR5->connect(
    'dbi:mysql:' . JASPAR_DB_NAME . ':' . JASPAR_DB_HOST,
    JASPAR_DB_USER,
    JASPAR_DB_PASS
);

die "Could not connect to JASPAR database\n" if !$jdb;

#
# Get various database adaptors
#
my $dbia = $opdb->get_DBInfoAdaptor;
die "Error getting DBInfoAdaptor\n" if !$dbia;

my $cla = $opdb->get_ConservationLevelAdaptor;
die "Error getting ConservationLevelAdaptor\n" if !$cla;

my $thla = $opdb->get_ThresholdLevelAdaptor;
die "Error getting ThresholdLevelAdaptor\n" if !$thla;

my $srla = $opdb->get_SearchRegionLevelAdaptor;
die "Error getting SearchRegionLevelAdaptor\n" if !$srla;

my $ga = $opdb->get_GeneAdaptor;
die "Error getting GeneAdaptor\n" if !$ga;

my $pa = $opdb->get_PromoterAdaptor;
die "Error getting PromoterAdaptor\n" if !$pa;

my $xgidta = $opdb->get_ExternalGeneIDTypeAdaptor();
die "Error getting ExternalGeneIDTypeAdaptor\n" if !$xgidta;

my $crla = $opdb->get_ConservedRegionLengthAdaptor();
die "Error getting ConservedTFBSAdaptor\n" if !$crla;

my $ctfsa = $opdb->get_ConservedTFBSAdaptor();
die "Error getting ConservedTFBSAdaptor\n" if !$ctfsa;

my $aca = $opdb->get_AnalysisCountsAdaptor();
die "Error getting AnalysisCountsAdaptor\n" if !$aca;

my $db_info = $dbia->fetch_db_info();
die "Error fetching DB info\n" if !$db_info;

#
# Check arguments or set to defaults
#
my $gidts = $xgidta->fetch_external_gene_id_type_hash();
if ($id_type) {
    unless ($gidts->{$id_type}) {
        my $err = "Invalid gene ID type $id_type specified;";
        $err .= " valid types are:\n";
        foreach my $gidt (sort keys %$gidts) {
            $err .= sprintf "%d - %s\n", $gidt, $gidts->{$gidt}->name();
        }

        die $err;
    }
} else {
    $id_type = DFLT_ID_TYPE;
}

my $cons_level_hash = $cla->fetch_conservation_level_hash();
if ($cons_level) {
    unless ($cons_level_hash->{$cons_level}) {
        my $err = "Invalid conservation level $cons_level; valid levels are:\n";
        foreach my $cl (sort keys %$cons_level_hash) {
            $err .= sprintf "%d - %s\n",
                $cl, $cons_level_hash->{$cl}->min_conservation();
        }

        die $err;
    }
} else {
    $cons_level = DFLT_CONSERVATION_LEVEL;
}

my $thresh_level_hash = $thla->fetch_threshold_level_hash();
if ($thresh_level) {
    unless ($thresh_level_hash->{$thresh_level}) {
        my $err = "Invalid threshold level $thresh_level; valid levels are:\n";
        foreach my $thl (sort keys %$thresh_level_hash) {
            $err .= sprintf "%d - %s\n",
                $thl, $thresh_level_hash->{$thl}->threshold();
        }

        die $err;
    }
} else {
    $thresh_level = DFLT_THRESHOLD_LEVEL;
}

my $sr_level_hash = $srla->fetch_search_region_level_hash();
if ($sr_level) {
    unless ($sr_level_hash->{$sr_level}) {
        my $err = "Invalid search region level $sr_level; valid levels are:\n";
        foreach my $srl (sort keys %$sr_level_hash) {
            $err .= sprintf "%d - %d-%d\n",
                $srl,
                $sr_level_hash->{$srl}->upstream_bp(),
                $sr_level_hash->{$srl}->downstream_bp();
        }

        die $err;
    }
} else {
    $sr_level = DFLT_SEARCH_REGION_LEVEL;
}

#
# Fetch TF set.
#
my $tf_set = fetch_tf_set($collection, \@tax_groups, $min_ic);

#
# Get the target set of genes
#
my $t_in_gene_ids = read_gene_ids($t_gene_file);
if (!$t_in_gene_ids) {
    die "Error reading gene IDs from $t_gene_file\n";
}

my ($t_gids,
    $t_included_gene_ids,
    $t_missing_gene_ids,
    $t_gid_ensids,
    $t_gid_gene_ids,
    $t_gene_id_gids)    = fetch_opossum_gene_ids($t_in_gene_ids, $id_type);

#
# Get the background set of genes
#
my $bg_gids;
my $bg_included_gene_ids;
my $bg_missing_gene_ids;
if ($num_rand_bg_genes) {
    ($bg_gids,
     $bg_included_gene_ids) = fetch_random_gene_ids($num_rand_bg_genes);

    if (!$bg_gids) {
        die "Error fetching random set of background genes\n";
    }
} elsif ($bg_gene_file) {
    my $bg_in_gene_ids = read_gene_ids($bg_gene_file);
    ($bg_gids,
     $bg_included_gene_ids,
     $bg_missing_gene_ids)  = fetch_opossum_gene_ids($bg_in_gene_ids, $id_type);
} else {
    #$bg_gids = $ga->fetch_gene_ids();
}

#
# Get the target and background TFBS counts (number of times for which
# each TFBS appears for each gene at the conservation, PWM threshold and
# search region levels specifies by the corresponding levels).
#
my $tf_ids = $tf_set->ids;

my $t_ac = $aca->fetch_counts(
    -tf_ids              => $tf_ids,
    -gene_ids            => $t_gids,
    -conservation_level  => $cons_level,
    -threshold_level     => $thresh_level,
    -search_region_level => $sr_level
);

if (!$t_ac) {
    die "Error fetching default TFBS counts for target genes\n";
}

my $bg_ac = $aca->fetch_counts(
    -tf_ids              => $tf_ids,
    -gene_ids            => $bg_gids,
    -conservation_level  => $cons_level,
    -threshold_level     => $thresh_level,
    -search_region_level => $sr_level
);

if (!$bg_ac) {
    die "Error fetching default TFBS counts for background genes\n";
}

my $t_cons_reg_len = $crla->fetch_total_length(
    -conservation_level     => $cons_level,
    -search_region_level    => $sr_level,
    -gene_ids               => $t_gids
);

my $bg_cons_reg_len = $crla->fetch_total_length(
    -conservation_level     => $cons_level,
    -search_region_level    => $sr_level,
    -gene_ids               => $bg_gids
);

my $fisher = OPOSSUM::Analysis::Fisher->new();
die "Error initializing Fisher analysis\n" if !$fisher;

my $fresults = $fisher->calculate_Fisher_probability($bg_ac, $t_ac);
die "Error computing Fisher probability\n" if !$fresults;

my $zscore = OPOSSUM::Analysis::Zscore->new();
die "Error creating Zscore analysis object\n" if !$zscore;

my $zresults = $zscore->calculate_Zscore(
    $bg_ac,
    $t_ac,
    $bg_cons_reg_len,
    $t_cons_reg_len,
    $tf_set
);
die "Error computing Zscores\n" if !$zresults;

my $combined_results = OPOSSUM::Analysis::CombinedResultSet->new(
    -fisher_result_set  => $fresults,
    -zscore_result_set  => $zresults
);

write_results($out_file, $combined_results, $tf_set);

if ($hits_file) {
    write_tfbs_details($hits_file, $t_ac, $tf_set, $cons_level, $thresh_level,
        $sr_level, $t_included_gene_ids, $t_gid_ensids, $t_gid_gene_ids,
        $id_type);
}

exit;

sub fetch_tf_set
{
    my ($collection, $tax_groups, $min_ic) = @_;

    my $matrix_set;

    if ($tf_file) {
        open(FH, $tf_file) || die "Error opening TF file $tf_file - $!\n";

        $matrix_set = TFBS::MatrixSet->new();

        while (my $line = <FH>) {
            my $matrix;
            my $matrix_id;

            #
            # Try to automatically figure out whether TFBSs are specified as
            # internal oPOSSUM database IDs (numeric), external (JASPAR) IDs
            # of the form MA####, MF#### or PF#### or names (anything which
            # doesn't fit above two patterns)
            #
            chomp $line;
            if (   $line =~ /^\s*(MA\d{4})\s*$/
                || $line =~ /^\s*(MF\d{4})\s*$/
                || $line =~ /^\s*(PF\d{4})\s*$/)
            {
                # TFBSs specified by external (JASPAR) ID
                $matrix_id = $1;
                $matrix = $jdb->get_Matrix_by_ID($matrix_id);
            } elsif ($line =~ /^\s*(\S+)\s*$/) {
                # TFBSs specified by name
                $matrix_id = $1;
                $matrix = $jdb->get_Matrix_by_name($matrix_id);
            }

            if ($matrix) {
                $matrix_set->add_matrix($matrix);
            } else {
                print "Could not fetch TFBS profile $matrix_id\n";
            }
        }
        close(FH);
    } else {
        my %matrix_args;
        $matrix_args{-collection} = $collection if $collection;
        $matrix_args{-tax_group} = $tax_groups if $tax_groups;
        $matrix_args{-min_ic} = $min_ic if $min_ic;
        $matrix_args{-matrixtype} = 'PFM';

        $matrix_set = $jdb->get_MatrixSet(%matrix_args);
    }

    die "Could not fetch TFBS profiles\n"
        if !$matrix_set || $matrix_set->size == 0;

    my $tf_set = OPOSSUM::TFSet->new(
        -matrix_set => $matrix_set
    );

    return $tf_set;
}

#
# For the provided list of Ensembl / external gene IDs, fetch the associated
# oPOSSUM gene IDs from the DB. Also keep track of which provided gene IDs
# mapped (included) or did not map (missing) to oPOSSUM gene IDs.
#
# NOTE: There can be a 1-to-many mapping of oPOSSUM gene IDs to external gene
# IDs. Make sure all of these are retained in the included external gene IDs
# list.
#
sub fetch_opossum_gene_ids
{
    my ($in_gene_ids, $id_type) = @_;

    my $sql;
    if ($id_type == DFLT_ID_TYPE) {
        #
        # Ensembl IDs provided.
        #
        $sql = qq{select gene_id, ensembl_id from genes where ensembl_id = ?};
    } else {
        #
        # Fetch by external gene IDs
        #
        $sql = qq{
            select g.gene_id, g.ensembl_id
            from genes g, external_gene_ids xgi
            where xgi.gene_id = g.gene_id
            and xgi.id_type = $id_type
            and xgi.external_id = ?
        };
    }

    my $sth = $ga->prepare($sql);
    if (!$sth) {
        die "Error preparing SQL statement:\n\t$sql\n";
    }

    my @gids;
    my @included_gene_ids;
    my @missing_gene_ids;
    my %gid_included;
    my %gene_id_included;
    my %gid_ensids;     # mapping of oPOSSUM gene IDs to Ensembl IDs
    my %gid_gene_ids;   # mapping of oPOSSUM gene IDs to external IDs
    my %gene_id_gids;   # mapping of external / Ensembl IDs to oPOSSUM gene IDs
    foreach my $in_gene_id (@$in_gene_ids) {
        if (!$sth->execute($in_gene_id)) {
            die "Error executing fetch gene IDs for $in_gene_id:\n"
                . $sth->errstr . "\n";
        }

        my $gid;
        my $ensid;
        while (my @row = $sth->fetchrow_array()) {
            $gid   = $row[0];
            $ensid = $row[1];

            unless ($gid_included{$gid}) {
                push @gids, $gid;

                #
                # There is a 1-to-1 mapping of oPOSSUM gene IDs
                # and Ensembl IDs.
                #
                $gid_ensids{$gid} = $ensid;

                $gid_included{$gid} = 1;
            }

            unless ($gene_id_included{$in_gene_id}) {
                push @included_gene_ids, $in_gene_id;

                #
                # Each external gene ID can have only one oPOSSUM gene ID
                # mapped to it.
                #
                $gene_id_gids{$in_gene_id} = $gid;

                $gene_id_included{$in_gene_id} = 1;
            }

            #
            # A single oPOSSUM gene ID can have multiple external gene IDs
            # mapped to it.
            #
            unless ($id_type == DFLT_ID_TYPE) {
                push @{$gid_gene_ids{$gid}}, $in_gene_id;
            }
        }
    }
    $sth->finish();

    #
    # Determine which of the input Ensembl / external gene IDs where not mapped
    # to an oPOSSUM gene ID (missing).
    #
    foreach my $in_gene_id (@$in_gene_ids) {
        unless ($gene_id_included{$in_gene_id}) {
            push @missing_gene_ids, $in_gene_id;
        }
    }

    return (
        @gids              ? \@gids              : undef,
        @included_gene_ids ? \@included_gene_ids : undef,
        @missing_gene_ids  ? \@missing_gene_ids  : undef,
        %gid_ensids        ? \%gid_ensids        : undef,
        %gid_gene_ids      ? \%gid_gene_ids      : undef,
        %gene_id_gids      ? \%gene_id_gids      : undef
    );
}

#
# Generate a list of random gene IDs from the oPOSSUM DB
#
sub fetch_random_gene_ids
{
    my ($num_genes) = @_;

    return if !$num_genes;

    my $sql = qq{
        select gene_id, ensembl_id from genes
        order by rand() limit $num_genes
    };

    my $sth = $ga->prepare($sql)
        || die "Error preparing SQL statement:\n\t$sql\n";

    $sth->execute() || die "Error executing SQL statement:\n\t$sql\n";

    my @gids;
    my @ensids;
    while (my @row = $sth->fetchrow_array()) {
        push @gids,   $row[0];
        push @ensids, $row[1];
    }

    return (
        @gids   ? \@gids : undef,
        @ensids ? \@ensids : undef
    );
}

#
# Read the gene IDs from the specified file
#
sub read_gene_ids
{
    my ($gene_file) = @_;

    open(FH, $gene_file) || die "Error opening gene file $gene_file - $!\n";

    my @gene_ids;
    my %gene_ids;
    while (my $line = <FH>) {
        chomp $line;
        $line =~ tr/\r//d;  ## get rid of the windows line-ender, if necessary
        if ($id_type == DFLT_ID_TYPE) {
            if ($line =~ /^(ENS\S*G\d+)/) {
                my $ens_id = $1;
                #
                # If a duplicate Ensembl ID is in the input file
                # it is only included once in the analysis.
                #
                if ($gene_ids{$ens_id}) {
                    print "Duplicate Ensembl ID $ens_id in input file"
                        . " - ignoring\n";
                    next;
                }
                push @gene_ids, $ens_id;
                $gene_ids{$ens_id} = 1;
            } else {
                die "$gene_file contains the following line which does not"
                    . " appear to contain a valid Ensembl ID:\n$line\n";
            }
        } else {
            if ($line =~ /^(\S+)$/) {
                my $ext_id = $1;
                #
                # If a duplicate external ID is in the input file
                # it is only included once in the analysis.
                #
                if ($gene_ids{$ext_id}) {
                    print "duplicate ID $ext_id in input file"
                        . " - ignoring\n";
                    next;
                }
                push @gene_ids, $ext_id;
                $gene_ids{$ext_id} = 1;
            } else {
                die "$gene_file contains the following line which does not"
                    . " appear to be a valid external ID:\n$line\n";
            }
        }
    }

    close(FH);

    return @gene_ids ? \@gene_ids : undef;
}

#
# Ouput combined z-score/Fisher results
#
sub write_results
{
    my ($filename, $result_set, $tf_set) = @_;

    if ($filename) {
        open(FH, ">$filename")
            || die "Error opening output results file $filename - $!\n";
    } else {
        open(FH, ">-");
    }

    my $results = $result_set->get_list(-sort_by => 'zscore', -reverse => 1);

    printf FH "TF Name\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";
    foreach my $result (@$results) {
        my $tf = $tf_set->get_tf($result->id());

        printf FH 
            "%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",
            $tf->{name},
            $tf->ID,
            $tf->class() || 'N/A',
            $tf->tag('family') || 'N/A',
            $tf->tag('tax_group') || 'N/A',
            $tf->to_ICM->total_ic(),
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_tfbs_hits() || 0,
            defined $result->t_tfbs_rate()
                ? sprintf("%.3f", $result->t_tfbs_rate()) : 'N/A',
            $result->bg_tfbs_hits() || 0,
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
# Write the details of the putative TFBSs for each TF/gene.
#
sub write_tfbs_details
{
    my ($filename, $ac, $tf_set, $cons_level, $thresh_level, $sr_level,
        $gene_ids, $gid_ensids, $gid_gene_ids, $gene_id_type) = @_;

    open(HFH, ">$filename") || die "Error creating TFBS details file - $!\n";

    my $gids          = $ac->gene_ids;
    my $tf_ids        = $tf_set->ids;

    my $threshold     = $thresh_level_hash->{$thresh_level}->threshold();
    my $srl           = $sr_level_hash->{$sr_level};
    my $upstream_bp   = $srl->upstream_bp;
    my $downstream_bp = $srl->downstream_bp;

    my %gid_sites;
    my %gid_proms;
    my %gid_srs;
    my %gid_chrs;
    my %gid_starts;
    my %gid_ends;
    my %gid_strands;
    foreach my $gid (@$t_gids) {
        my $gene = $ga->fetch_by_gene_id($gid);

        $gid_chrs{$gid}    = $gene->chr;
        $gid_starts{$gid}  = $gene->start;
        $gid_ends{$gid}    = $gene->end;
        $gid_strands{$gid} = $gene->strand;

        my $promoters = $gene->fetch_promoters();

        $gid_proms{$gid} = $promoters;

        #my $search_regions = $gene->promoter_search_regions(
        #    $upstream_bp, $downstream_bp
        #);
        #
        #$gid_srs{$gid} = $search_regions;
    }

    foreach my $tf_id (@$tf_ids) {
        my $tf = $tf_set->get_tf($tf_id);

        my $tf_name = $tf->name();

        print HFH "$tf_name associated genes:\n";

        if ($id_type == DFLT_ID_TYPE) {
            print HFH "Ensembl ID";
        } else {
            printf HFH "Gene ID\t%s", $gidts->{$id_type}->name();
        }

        print HFH "\tChr\tStrand\tPromoter Start\tPromoter End\tTSS\tTFBS Start\tTFBS End\tTFBS Rel. Start\tTFBS Rel. End\tTFBS Orientation\tScore\tRel. Score\tSeq\n";

        foreach my $gid (@$gids) {
            next if !$ac->gene_tfbs_count($gid, $tf_id);

            my $gene_start = $gid_starts{$gid};
            my $gene_end   = $gid_ends{$gid};
            my $chr        = $gid_chrs{$gid};
            my $strand     = $gid_strands{$gid};
            my $ens_id     = $gid_ensids->{$gid};

            my $sites = $ctfsa->fetch_by_gene(
                -gene_id            => $gid,
                -tf_id              => $tf_id,
                -conservation_level => $cons_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );

            my $first = 1;
            foreach my $site (@$sites) {
                my $site_start      = $site->start;
                my $site_end        = $site->end;
                my $site_seq        = $site->seq;
                my $site_strand     = $site->strand;
                my $site_score      = $site->score;
                my $site_rel_score  = $site->rel_score * 100;

                my $prev_start      = 0;
                my $prev_end        = 0;
                my $prev_site_start = 0;
                my $printed         = 0;

                foreach my $prom (@{$gid_proms{$gid}}) {
                    my $tss = $prom->tss;
                    my $prom_start;
                    my $prom_end;
                    if ($strand == 1) {
                        $prom_start = $tss - $upstream_bp;
                        $prom_end   = $tss + $downstream_bp - 1;
                    } else {
                        $prom_start = $tss - $downstream_bp + 1;
                        $prom_end   = $tss + $upstream_bp;
                    }

                    $prom_start = $gene_start if $prom_start < $gene_start;
                    $prom_end   = $gene_end if $prom_end > $gene_end;

                    # 1-based sequence coords
                    my $reg_start = $prom_start - $gene_start + 1;
                    my $reg_end   = $prom_end - $gene_start + 1;
                    my $rel_tss   = $tss - $gene_start + 1;

                    next if (   $reg_start  == $prev_start
                             && $reg_end    == $prev_end
                             && $site_start == $prev_site_start
                    );

                    if (   $site_start >= $reg_start
                        && $site_end <= $reg_end)
                    {
                        if ($id_type == DFLT_ID_TYPE) {
                            if ($first) {
                                #printf HFH "%-18s  %2s  %2d",
                                printf HFH "%s\t%s\t%d\t",
                                    $ens_id,
                                    $chr,
                                    $strand;
                                $first = 0;
                            } else {
                                #printf HFH "                          ";
                                printf HFH "\t\t\t";
                            }
                        } else {
                            if ($first) {
                                #printf HFH "%-15s  %-18s  %2s  %2d",
                                printf HFH "%s\t%s\t%s\t%d\t",
                                    join ',', @{$gid_gene_ids->{$gid}},
                                    $ens_id,
                                    $chr,
                                    $strand;
                                $first = 0;
                            } else {
                #printf HFH "                                           ";
                                printf HFH "\t\t\t\t";
                            }
                        }

                        my ($rel_start, $rel_end);
                        if ($strand == 1) {
                            $rel_start = $site_start - $rel_tss;
                            if ($site_start >= $rel_tss) {
                                $rel_start++;
                            }

                            $rel_end = $site_end - $rel_tss;
                            if ($site_end >= $rel_tss) {
                                $rel_end++;
                            }
                        } else {
                            $rel_start = $rel_tss - $site_start;
                            if ($site_start <= $rel_tss) {
                                $rel_start++;
                            }

                            $rel_end = $rel_tss - $site_end;
                            if ($site_end <= $rel_tss) {
                                $rel_end++;
                            }

                            ($rel_start, $rel_end) = ($rel_end, $rel_start);
                        }

                        if ($printed) {
                            #print HFH '  *';
                            print HFH '*';
                        } else {
                            #print HFH '   ';
                        }

                        #printf HFH "%9d  %9d  %9d  %25s  %9d  %5d  %9d  %5d  %2d  %.3g\n",
                        printf HFH
                            "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3g\t%.1f%%\t%s\n",
                            $prom_start,
                            $prom_end,
                            $tss,
                            $site_start + $gene_start - 1,
                            $site_end + $gene_start - 1,
                            $rel_start,
                            $rel_end,
                            $site_strand,
                            $site_score,
                            $site_rel_score,
                            $site_seq;

                        $prev_start      = $reg_start;
                        $prev_end        = $reg_end;
                        $prev_site_start = $site_start;

                        $printed = 1;
                    }
                }
            }
        }
        print HFH "\n";
    }
    close(HFH);
}
