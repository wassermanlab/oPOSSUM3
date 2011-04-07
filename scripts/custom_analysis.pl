#!/usr/local/bin/perl -w

=head1 NAME

custom_analysis.pl

=head1 SYNOPSIS

  custom_analysis.pl -g gene_file [-s species] [-i gene_ID_type]
        [-bg bg_gene_file | -n num_rand_bg_genes]
        [-db tf_database]
        [-t tf_file | -ic IC | -tax tax_groups]
        [-cl cons_level] [-th threshold]
        [-u upstream_bp -d downstream_bp]
        [-h tfbs_hits_file] [-nx]
        [-o out_file]


=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -g gene_file		= File containing gene IDs to analyze. Format
   			  should be one ID per line. The gene ID type
			  may be specified by the -i switch, i.e.
			  'ensembl', 'hugo' etc.
   -s species		= Species on which to base the analysis, i.e.
   			  'human', 'mouse'. The gene IDs/types must
			  be valid for this species.
			  Default = 'human'
   -i gene_ID_type	= Type of gene ID contained in the gene_file,
   			  i.e. 'ensembl', 'hugo' etc.
			  Default = 'ensembl'
   -bg bg_gene_file	= File containing gene IDs of background set of
   			  genes. Format is one ID per line. Assumes the
			  same gene ID type as the target genes. If
			  neither -bg nor -n is provided all genes in the
			  oPOSSUM DB are used as background.
   -n num_rand_bg_genes	= Generate a random background set of genes from
   			  the oPOSSUM database of size num_rand_bg_genes.
			  If neither -bg nor -n is provided all genes in the
			  oPOSSUM DB are used as background.
   -db tf_database	= Specify which TF database to use
   			  (default = JASPAR_CORE)
   -t tf_file		= File containing TF IDs/names. The IDs/names
   			  should appear one per line. Either internal
			  oPOSSUM DB IDs (numeric), external (JASPAR)
			  IDs of the form MA#### or names (anything else)
			  may be used. The -t, -ic and -tax options are
			  exclusive. If none are specified, use all TFBS
			  profiles in the oPOSSUM database.
   -ic IC		= Minimum information content (specificity) of
   			  TFBS profile matrices to use for analysis. The
			  -t, -ic and -tax options are exclusive. If none
			  are specified, use all TFBS profiles in the
			  oPOSSUM database.
			  Minimum valid IC = 8 bits
   -tax tax_groups	= List of taxonomic supergroups of TFBS profiles
   			  to use for analysis. The -t, -ic and -tax options
			  are exclusive. If none are specified, use all
			  TFBS profiles in the oPOSSUM database.
			  Valid options are:
			  	'plant', 'vertebrate' or 'insect'
   -cl cons_level	= Level of conservation for which to perform the
   			  analysis.
			  Levels are:
			    1 = Top 30% of conserved regions with an
				absolute minimum percent identity of 60%
			    2 = Top 20% of conserved regions with an
				absolute minimum percent identity of 65%
			    3 = Top 10% of conserved regions with an
				absolute minimum percent identity of 70%
			  Default = 3
   -th threshold	= Minimum relative TFBS position weight matrix
   			  (PWM) score to report in the analysis. The
			  thresold may be spefifies as a percentage, i.e.
			  '80%' or a decimal, i.e. 0.8.
			  Default = 80%
			  Min. = 75%
   -u upstream_bp	= Amount of upstream sequence to include in the
   			  search region for TFBSs.
			  Default = 5000
			  Max. = 10000
   -d downstream_bp	= Amount of downstream sequence to include in the
   			  search region for TFBSs.
			  Default = 5000
			  Max. = 10000
   -h tfbs_hits_file    = Name of file to which details of the putative
   			  TFBSs are written.
   -nx			= If switch is provided, do NOT exclude results
   			  in which a given TFBS was only found for a
			  single gene.
			  Default = exclude single gene hits
   -o out_file		= Name of output file to which analysis results
			  are written.
			  Default = standard output
   			  

=head1 DESCRIPTION

Perform the same analysis as through the oPOSSUM website using custom
parameter selection (conservation level, transcription factor binding site
(TFBS) position weight matrix (PWM) score threshold and search region).

Take a list of genes as provided in the input genes file. Optionally also
take a list of background genes or a number of randomly generated genes or
all the genes in the oPOSSUM database. Optionally take a subset of
transcription factors (TFs) either specified in an input file, or limited
by external (JASPAR) database name and information content or taxonomic
supergroup or all TFs in the oPOSSUM database. Also optionally take arguments
specifying conservation level, PWM score threshold and amount of
upstream/downstream sequence to include in searching for TFBSs. 

Count the number of TFBSs for each TF which was found at the given
conservation level and PWM score threshold within the given search region
for both the foreground and background set of genes. Perform Fisher exact
test and z-score analysis and output these results to the output file.

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
use lib '/space/devel/oPOSSUM/lib';

use Getopt::Long;
use Pod::Usage;

use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFInfoSet;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Counts;

use constant OPOSSUM_DB_HOST => "vm5.cmmt.ubc.ca";
use constant OPOSSUM_DB_USER => "opossum_r";
use constant OPOSSUM_DB_PASS => undef;
use constant OPOSSUM_DB_NAME => "oPOSSUM_2008";

use constant DFLT_CONSERVATION_LEVEL     => 3;
use constant DFLT_THRESHOLD              => 0.80;
use constant MIN_THRESHOLD               => 0.75;
use constant DFLT_UPSTREAM_BP            => 5000;
use constant DFLT_DOWNSTREAM_BP          => 5000;
use constant MAX_UPSTREAM_BP             => 10000;
use constant MAX_DOWNSTREAM_BP           => 10000;
use constant DFLT_SPECIES_NUM            => 1;
use constant DFLT_ID_TYPE_NAME           => 'ensembl';
use constant DFLT_NO_EXCLUDE_SINGLE_HITS => 0;
use constant TF_DBS     => ['jaspar_core', 'jaspar_phylofacts'];
use constant DFLT_TF_DB => 'jaspar_core';

my $gene_file;
my $species_name;
my $id_type_name;
my $bg_gene_file;
my $num_rand_bg_genes;
my $tf_db;
my $tf_file;
my $min_ic;
my @tax_groups;
my $cons_level;
my $threshold;
my $upstream_bp;
my $downstream_bp;
my $out_file;
my $hits_file;
my $no_exclude_single_hits = DFLT_NO_EXCLUDE_SINGLE_HITS;
GetOptions(
    'g=s'   => \$gene_file,
    's=s'   => \$species_name,
    'i=s'   => \$id_type_name,
    'bg=s'  => \$bg_gene_file,
    'n=i'   => \$num_rand_bg_genes,
    'db=s'  => \$tf_db,
    't=s'   => \$tf_file,
    'ic=f'  => \$min_ic,
    'tax=s' => \@tax_groups,
    'cl=i'  => \$cons_level,
    'th=s'  => \$threshold,
    'u=i'   => \$upstream_bp,
    'd=i'   => \$downstream_bp,
    'o=s'   => \$out_file,
    'h=s'   => \$hits_file,
    'nx'    => \$no_exclude_single_hits
);

my $exclude_single_hits = $no_exclude_single_hits ? 0 : 1;

if (!$gene_file) {
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

if (!$tf_db) {
    $tf_db = DFLT_TF_DB;
}

#$tf_db = uc $tf_db;

# In case TF DB was just specified as 'CORE' or 'PHYLOFACTS'
if ($tf_db =~ /^core/i || $tf_db =~ /^phylofacts/i) {
    $tf_db = 'jaspar_' . $tf_db;
    $tf_db = lc $tf_db;
}

#my $tf_dbs = TF_DBS;
#if (!grep /^$tf_db$/, @$tf_dbs) {
#    pod2usage(
#        -msg		=> sprintf("\nUnknown TF database $tf_db"
#				    . "\nMust be one of: %s\n\t",
#				    join ', ', @$tf_dbs),
#	-verbose	=> 1);
#}

if ($threshold) {
    if ($threshold =~ /(\S+)%/) {
        $threshold = $1 / 100;
    }
} else {
    $threshold = DFLT_THRESHOLD;
}

if ($upstream_bp) {
    if ($upstream_bp > MAX_UPSTREAM_BP) {
        die "Upstream search region is larger than maximum allowed "
            . MAX_UPSTREAM_BP . "\n";
    } elsif ($upstream_bp < 0) {
        die "Upstream search region cannot be negative\n";
    }
} else {
    $upstream_bp = DFLT_UPSTREAM_BP;
}

if ($downstream_bp) {
    if ($downstream_bp > MAX_DOWNSTREAM_BP) {
        die "Downstream search region is larger than maximum allowed "
            . MAX_DOWNSTREAM_BP . "\n";
    } elsif ($downstream_bp < 0) {
        die "Downstream search region cannot be negative\n";
    }
} else {
    $downstream_bp = DFLT_DOWNSTREAM_BP;
}

if ($upstream_bp + $downstream_bp <= 0) {
    die "Search region specified must be greater than 0 length\n";
}

if ($tf_file) {
    if ($min_ic || @tax_groups) {
        die "\nIf specifying TFBSs, please provide only one of either a TFBS"
            . " file, a minimum information content or one or more taxonomic"
            . " supergroups\n";
    }
} elsif ($min_ic && @tax_groups) {
    die "\nIf specifying TFBSs, please provide only one of either a TFBS"
        . " file, a minimum information content or one or more taxonomic"
        . " supergroups\n";
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
    -dbname   => OPOSSUM_DB_NAME
);

die "Could not connect to oPOSSUM database\n" if !$opdb;

#
# Get various database adaptors
#
my $dbbia = $opdb->get_DBBuildInfoAdaptor;
die "Error getting DBBuildInfoAdaptor\n" if !$dbbia;

my $cla = $opdb->get_ConservationLevelAdaptor;
die "Error getting ConservationLevelAdaptor\n" if !$cla;

my $thla = $opdb->get_ThresholdLevelAdaptor;
die "Error getting ThresholdLevelAdaptor\n" if !$thla;

my $srla = $opdb->get_SearchRegionLevelAdaptor;
die "Error getting SearchRegionLevelAdaptor\n" if !$srla;

my $gpa = $opdb->get_GenePairAdaptor;
die "Error getting GenePairAdaptor\n" if !$gpa;

my $ppa = $opdb->get_PromoterPairAdaptor;
die "Error getting PromoterPairAdaptor\n" if !$ppa;

my $tia = $opdb->get_TFInfoAdaptor();
die "Error getting TFInfoAdaptor\n" if !$tia;

my $xgita = $opdb->get_ExternalGeneIDTypeAdaptor;
die "Error getting ExternalGeneIDTypeAdaptor\n" if !$xgita;

#my $ctfsa = $opdb->get_ConservedTFBSAdaptor(lc $tf_db);
my $ctfsa = $opdb->get_ConservedTFBSAdaptor($tf_db);
die "Error getting ConservedTFBSAdaptor\n" if !$ctfsa;

#my $aca = $opdb->get_AnalysisCountsAdaptor(lc $tf_db);
my $aca = $opdb->get_AnalysisCountsAdaptor($tf_db);
die "Error getting AnalysisCountsAdaptor\n" if !$aca;

my $db_build_info = $dbbia->fetch_db_build_info();
die "Error fetching DB build info\n" if !$db_build_info;

#
# Check arguments or set to defaults
#
my $species_num;
if ($species_name) {
    $species_name = lc $species_name;
    if ($species_name eq $db_build_info->species1) {
        $species_num = 1;
    } elsif ($species_name eq $db_build_info->species2) {
        $species_num = 2;
    } else {
        die sprintf("Species is invalid; please specify %s or %s\n",
            $db_build_info->species1, $db_build_info->species2);
    }
} else {
    $species_num = DFLT_SPECIES_NUM;
}

my $id_type;
if ($id_type_name) {
    my $gidts = $xgita->fetch_names;
    if (!grep(/$id_type_name/i, @$gidts)) {
        die sprintf(
            "Invalid gene ID type $id_type_name; valid types are %s\n",
            join ", ", @$gidts);
    }

    my $ext_gene_id_type = $xgita->fetch_by_name($id_type_name);
    if (!$ext_gene_id_type) {
        die "Error fetching external gene ID type for $id_type_name\n";
    }
    $id_type = $ext_gene_id_type->id_type;
} else {
    $id_type_name = DFLT_ID_TYPE_NAME;
    $id_type      = 0;
}

if ($cons_level) {
    my $levels = $cla->fetch_levels;
    if (!grep(/^$cons_level$/, @$levels)) {
        die sprintf(
            "Invalid conservation level $cons_level; valid levels are %s\n",
            join ", ", @$levels);
    }
} else {
    $cons_level = DFLT_CONSERVATION_LEVEL;
}

#
# Fetch TFBS info set.
#
my $tf_info_set = fetch_tf_info_set();

#
# Get the target set of genes
#
my ($t_gpids, $t_gene_ids, $t_missing_gene_ids) =
    get_gene_ids_from_file($gene_file, $id_type);

#
# Get the background set of genes
#
my $bg_gpids;
my $bg_gene_ids;
my $bg_missing_gene_ids;
if ($num_rand_bg_genes) {
    $bg_gpids = get_random_gene_pair_ids($num_rand_bg_genes);
    if (!$bg_gpids) {
        die "Error generating background set of genes\n";
    }

    @$bg_gene_ids = @$bg_gpids;
} elsif ($bg_gene_file) {
    ($bg_gpids, $bg_gene_ids, $bg_missing_gene_ids) =
        get_gene_ids_from_file($bg_gene_file, $id_type);
} else {
    $bg_gpids     = $gpa->fetch_gene_pair_ids();
    @$bg_gene_ids = @$bg_gpids;
}

#
# Get the target and background TFBS counts (number of times for which
# each TFBS appears for each PromoterPair at the conservation, PWM
# threshold and search region levels specifies by the corresponding levels).
#
my @tf_ids = $tf_info_set->tf_ids;

my $t_ac = $aca->fetch_custom_counts(
    -tf_ids             => \@tf_ids,
    -gene_pair_ids      => $t_gpids,
    -conservation_level => $cons_level,
    -threshold          => $threshold,
    -upstream_bp        => $upstream_bp,
    -downstream_bp      => $downstream_bp
);
if (!$t_ac) {
    die "Error fetching custom TFBS counts for target genes\n";
}

$t_ac->tf_info_set($tf_info_set);

my $bg_ac = $aca->fetch_custom_counts(
    -tf_ids             => \@tf_ids,
    -gene_pair_ids      => $bg_gpids,
    -conservation_level => $cons_level,
    -threshold          => $threshold,
    -upstream_bp        => $upstream_bp,
    -downstream_bp      => $downstream_bp
);
if (!$bg_ac) {
    die "Error fetching custom TFBS counts for background genes\n";
}

$bg_ac->tf_info_set($tf_info_set);

my $fisher = OPOSSUM::Analysis::Fisher->new();
die "Error initializing Fisher analysis\n" if !$fisher;

my $fresults = $fisher->calculate_Fisher_probability($bg_ac, $t_ac);
die "Error computing Fisher probability\n" if !$fresults;

my $zscore = OPOSSUM::Analysis::Zscore->new();
die "Error creating Zscore analysis object\n" if !$zscore;

my $zresults = $zscore->calculate_Zscore($bg_ac, $t_ac);
die "Error computing Zscores\n" if !$zresults;

# sort results by z-score
$zresults->sort_by('z_score', 1);

my $formatted_results = format_results($fresults, $zresults, $tf_info_set);

if ($out_file) {
    write_results($out_file, $formatted_results);
} else {
    write_results(undef, $formatted_results);
}

if ($hits_file) {
    my $sorted_tf_ids = sort_tf_ids_by_zscore_result($zresults);
    write_tfbs_details(
        $hits_file,     $id_type,   $sorted_tf_ids,
        $t_gene_ids,    $t_ac,      $upstream_bp,
        $downstream_bp, $threshold, $cons_level
    );
}

exit;

sub fetch_tf_info_set
{
    my $tf_info_set;

    if ($tf_file) {
        open(FH, $tf_file) || die "Error opening TFBS file $tf_file - $!\n";

        $tf_info_set = OPOSSUM::TFInfoSet->new();
        die "Error fetching master TFBS set\n" if !$tf_info_set;

        while (my $line = <FH>) {
            chomp $line;
            #
            # Try to automatically figure out whether TFBSs are specified as
            # internal oPOSSUM database IDs (numeric), external (JASPAR) IDs
            # of the form MA####, MF#### or PF#### or names (anything which
            # doesn't fit above two patterns)
            #
            if (   $line =~ /^\s*(MA\d{4})\s*$/
                || $line =~ /^\s*(MF\d{4})\s*$/
                || $line =~ /^\s*(PF\d{4})\s*$/)
            {
                # TFBSs specified by external (JASPAR) ID
                my $tf_info = $tia->fetch_by_external_id($1);
                if ($tf_info) {
                    $tf_info_set->add_tf_info($tf_info);
                } else {
                    print "Could not fetch TFBS with external ID $1 from"
                        . " database - skipping\n";
                }
            } elsif ($line =~ /^\s*(\d+)\s*$/) {
                # TFBSs specified by oPOSSUM internal database (numeric) ID
                my $tf_info = $tia->fetch_by_id($1);
                if ($tf_info) {
                    $tf_info_set->add_tf_info($tf_info);
                } else {
                    print "Could not fetch TFBS with ID $1 from"
                        . " database - skipping\n";
                }
            } elsif ($line =~ /^\s*(\S+)\s*$/) {
                # TFBSs specified by name
                my $tf_info = $tia->fetch_by_name($1);
                if ($tf_info) {
                    $tf_info_set->add_tf_info($tf_info);
                } else {
                    print "Could not fetch TFBS with name $1 from"
                        . " database - skipping\n";
                }
            }
        }
        close(FH);
    } elsif ($min_ic) {
        $tf_info_set = $tia->fetch_tf_info_set_by_min_ic($min_ic, $tf_db);
        die "Could not fetch $tf_db TFs with min. IC $min_ic\n"
            if !$tf_info_set || $tf_info_set->size == 0;
    } elsif (@tax_groups) {
        $tf_info_set =
            $tia->fetch_tf_info_set_by_phylums(\@tax_groups, $tf_db);
        die "Could not fetch $tf_db TFs by taxonomic supergroup(s) '"
            . join(', ', @tax_groups)
            . "' - please make sure the supergroup(s) chosen are valid for"
            . " this TF database\n"
            if !$tf_info_set || $tf_info_set->size == 0;
    } else {
        # Get all TFBSs in the database
        $tf_info_set = $tia->fetch_tf_info_set_by_external_db($tf_db);
        die "Could not fetch $tf_db TFs\n"
            if !$tf_info_set || $tf_info_set->size == 0;
    }

    return $tf_info_set;
}

#
# Read a list of gene identifiers from the specified file and then fetch the
# associated GenePair IDs from the oPOSSUM DB. Also keep track of any gene
# IDs which were not found in the oPOSSUM DB.
#
sub get_gene_ids_from_file
{
    my ($filename, $id_type) = @_;

    my @gene_pair_ids;
    my @gene_ids;
    my @missing_gene_ids;
    my $gene_list;

    $gene_list = read_gene_ids($filename);
    if (!$gene_list) {
        die "Error fetching gene IDs from $filename\n";
    }

    my $sql;
    if ($id_type != 0) {
        #
        # Fetch by external gene IDs
        #
        $sql = qq{select gp.gene_pair_id, xgi.external_id
		    from gene_pairs gp, promoter_pairs pp,
		    	external_gene_ids xgi
		    where gp.gene_pair_id = pp.gene_pair_id
		    	and xgi.id_type = $id_type
			and xgi.species = $species_num
			and xgi.gene_pair_id = gp.gene_pair_id
			and xgi.external_id in ('};
    } else {
        #
        # Fetch by Ensembl IDs
        #
        $sql = qq{select gp.gene_pair_id, gp.ensembl_id$species_num
		    from gene_pairs gp, promoter_pairs pp
		    where gp.gene_pair_id = pp.gene_pair_id
		    	and gp.ensembl_id$species_num in ('};

    }
    $sql .= join '\',\'', @$gene_list;
    $sql .= '\')';

    my $sth = $gpa->prepare($sql);
    if (!$sth) {
        die "Error preparing fetch gene pair IDs:\n" . $sth->errstr . "\n";
    }

    if (!$sth->execute) {
        die "Error executing fetch gene pair IDs:\n" . $sth->errstr . "\n";
    }

    my %gene_pair_id_exists;
    my %gene_id_exists;
    while (my @row = $sth->fetchrow_array) {
        # If multiple (external) gene IDs map to a single gene pair ID
        # only the first one gets included in the gene_ids list and the
        # others will be placed on the missing_gene_ids list. This is
        # misleading.
        #next if grep /^$row[0]$/, @gene_pair_ids;
        next if $gene_pair_id_exists{$row[0]};

        push @gene_pair_ids, $row[0];
        $gene_pair_id_exists{$row[0]} = 1;
        push @gene_ids, $row[1];
        $gene_id_exists{$row[1]} = 1;
    }
    $sth->finish;

    foreach my $gene_id (@$gene_list) {
        #if (!grep /^$gene_id$/i, @gene_ids) {
        unless ($gene_id_exists{$gene_id}) {
            push @missing_gene_ids, $gene_id;
        }
    }

    return (
        @gene_pair_ids    ? \@gene_pair_ids    : undef,
        @gene_ids         ? \@gene_ids         : undef,
        @missing_gene_ids ? \@missing_gene_ids : undef
    );
}

#
# Generate a list of random GenePair IDs from the oPOSSUM DB
#
sub get_random_gene_pair_ids
{
    my ($num_genes) = @_;

    return if !$num_genes;

    my $gpids = $gpa->fetch_gene_pair_ids;

    my $total_num_genes = scalar @$gpids;
    $num_genes = $total_num_genes if $num_genes > $total_num_genes;

    #
    # Note: could use MySQL built-in RAND() function, e.g.:
    # 	select gene_pair_id from gene_pairs order by rand() limit $num_genes;
    #
    my @rand_gpids;
    my %rand_gpid_exists;
    my $i = 0;
    while ($i < $num_genes) {
        my $idx  = int(rand($total_num_genes - 1));
        my $gpid = $gpids->[$idx];

        #if (!grep /^$gpid$/, @rand_gpids) {
        unless ($rand_gpid_exists{$gpid}) {
            push @rand_gpids, $gpid;
            $rand_gpid_exists{$gpid} = 1;
            $i++;
        }
    }

    return @rand_gpids ? \@rand_gpids : undef;
}

#
# Read the gene IDs from the specified file
#
sub read_gene_ids
{
    my ($gene_file) = @_;

    open(FH, $gene_file) || die "Error opening gene file $gene_file - $!\n";

    my @gene_ids;
    my %gene_id_exists;
    while (my $line = <FH>) {
        chomp $line;
        if ($id_type == 0) {
            if ($line =~ /^(ENS\S*G\d+)/) {
                my $ens_id = $1;
                #
                # If a duplicate Ensembl ID is in the input file
                # it is only included once in the analysis.
                #
                #if (grep(/^$ens_id$/, @gene_ids)) {
                if ($gene_id_exists{$ens_id}) {
                    print "Duplicate Ensembl ID $ens_id in input file"
                        . " - ignoring\n";
                    next;
                }
                push @gene_ids, $ens_id;
                $gene_id_exists{$ens_id} = 1;
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
                #if (grep(/^$ext_id$/, @gene_ids)) {
                if ($gene_id_exists{$ext_id}) {
                    print "duplicate ID $ext_id in input file"
                        . " - ignoring\n";
                    next;
                }
                push @gene_ids, $ext_id;
                $gene_id_exists{$ext_id} = 1;
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
# Combine and format Fisher and z-score results
#
sub format_results
{
    my ($fisher_results, $zscore_results, $tf_set) = @_;

    return if !$fisher_results || !$zscore_results;

    my $num_results = $zscore_results->num_results;

    my @fresults;
    my $result_idx = 0;
    while ($result_idx < $num_results) {
        my $zscore_result = $zscore_results->get_result($result_idx++);
        last if !$zscore_result;

        if ($zscore_result->t_gene_hits > 1) {
            if ($zscore_result->bg_gene_hits == 0) {
                printf STDERR "Warning: no background gene hits for TF %s\n",
                    $zscore_result->id;
            }

            my $name        = '';
            my $external_id = '';
            my $class       = '';
            my $phylum      = '';
            my $ic          = '';
            my $tf          = $tf_set->get_tf_info($zscore_result->id);
            $name        = $tf->name;
            $external_id = $tf->external_id || '';
            $class       = $tf->class || '';
            $phylum      = $tf->phylum || '';
            $ic          = $tf->ic;

            my $zscore = $zscore_result->z_score;
            my $zscore_str;
            if (defined $zscore) {
                if ($zscore =~ /^inf/i || $zscore =~ /^nan$/i) {
                    $zscore_str = "N/A";
                } else {
                    $zscore_str = sprintf("%.4g", $zscore);
                }
            } else {
                $zscore_str = "";
            }

            my $fisher_result =
                $fisher_results->get_result_by_id($zscore_result->id);
            my $fisher_p_value = $fisher_result->p_value;

            push @fresults,
                {
                id              => $zscore_result->id,
                t_tfbs_hits     => $zscore_result->t_hits,
                bg_tfbs_hits    => $zscore_result->bg_hits,
                t_gene_hits     => $zscore_result->t_gene_hits,
                bg_gene_hits    => $zscore_result->bg_gene_hits,
                t_gene_no_hits  => $fisher_result->t_no_hits,
                bg_gene_no_hits => $fisher_result->bg_no_hits,
                t_rate          => sprintf("%.4f", $zscore_result->t_rate),
                bg_rate         => sprintf("%.4f", $zscore_result->bg_rate),
                z_score         => $zscore_str,
                p_value         => defined $zscore_result->p_value
                ? sprintf("%.3e", $zscore_result->p_value)
                : "",
                fisher_p_value => sprintf("%.3e", $fisher_p_value),
                name           => $name,
                external_id    => $external_id,
                class  => $class  || '',
                phylum => $phylum || '',
                ic => sprintf("%.3f", $ic)
                };
        }
    }

    return @fresults ? \@fresults : undef;
}

#
# Ouput combined z-score/Fisher results
#
sub write_results
{
    my ($filename, $summary) = @_;

    if ($filename) {
        open(FH, ">$filename")
            || die "Error opening output results file $filename - $!\n";
    } else {
        open(FH, ">-");
    }

#
# Multi-line (tab-delimited) header format
#
#printf FH "TF\tTF\tTF\tIC\tBackgrd\tBackgrd\tTarget\tTarget\tBackgrd\tBackgrd\tTarget\tTarget\tZ-score\tFisher\n";
#printf FH "\tClass\tS.group\t\tgene\tgene\tgene\tgene\tTFBS\tTFBS\tTFBS\tTFBS\t\tscore\n";
#printf FH "\t\t\t\thits\tno-hits\thits\tno-hits\thits\trate\thits\trate\t\t\n";
#
# Single line (tab-delimited) header format
#
    printf FH
        "TF\tTF Class\tTF Supergroup\tIC\tBackground gene hits\tBackground gene non-hits\tTarget gene hits\tTarget gene non-hits\tBackground TFBS hits\tBackground TFBS rate\tTarget TFBS hits\tTarget TFBS rate\tZ-score\tFisher score\n";
    foreach my $result (@$summary) {
        printf FH "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",
            $result->{name},
            $result->{class},
            $result->{phylum},
            $result->{ic},
            $result->{bg_gene_hits},
            $result->{bg_gene_no_hits},
            $result->{t_gene_hits},
            $result->{t_gene_no_hits},
            $result->{bg_tfbs_hits},
            $result->{bg_rate},
            $result->{t_tfbs_hits},
            $result->{t_rate},
            $result->{z_score},
            $result->{fisher_p_value};
    }

    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/gene.
#
sub write_tfbs_details
{
    my ($hits_file,     $id_type,   $tf_ids,
        $t_gene_ids,    $t_ac,      $upstream_bp,
        $downstream_bp, $threshold, $cons_level
    ) = @_;

    open(HFH, ">$hits_file") || die "Error creating TFBS details file - $!\n";

    my $t_gpids     = $t_ac->gene_pair_ids;
    my $tf_info_set = $t_ac->tf_info_set;

    my $sql;
    if ($id_type == 0) {
        $sql = qq{select ensembl_id$species_num
		    from gene_pairs
		    where gene_pair_id = ?};
    } else {
        #
        # The external_id portion of the where clause seems redundent but is
        # there because there are possible one to many relationships between
        # gene_pair_ids and external_ids and we have to make sure
        # the external_ids we get back are the same as entered by the user
        # (not just the first one selected which is associated with the
        # gene_pair_id.
        #
        $sql = qq{select gp.ensembl_id$species_num, xgi.external_id
		    from external_gene_ids xgi, gene_pairs gp
		    where xgi.gene_pair_id = gp.gene_pair_id
			and xgi.id_type = $id_type
			and xgi.species = $species_num
			and gp.gene_pair_id = ?
			and xgi.external_id in ('};
        $sql .= join "','", @$t_gene_ids;
        $sql .= "')";
    }
    my $sth = $opdb->prepare($sql);

    my %gpid_ensids;
    my %gpid_gids;
    foreach my $gpid (@$t_gpids) {
        $sth->execute($gpid);
        if (my ($ensid, $gid) = $sth->fetchrow_array) {
            $gpid_ensids{$gpid} = $ensid;
            if (defined $gid) {
                $gpid_gids{$gpid} = $gid;
            }
        }
    }
    $sth->finish;

    my %gpid_sites;
    my %gpid_pps;
    my %gpid_srs;
    my %gpid_chrs;
    my %gpid_strands;
    foreach my $gpid (@$t_gpids) {
        my $gp = $gpa->fetch_by_gene_pair_id($gpid);
        if ($species_num == 1) {
            $gpid_chrs{$gpid}    = $gp->chr1;
            $gpid_strands{$gpid} = $gp->strand1;
        } elsif ($species_num == 2) {
            $gpid_chrs{$gpid}    = $gp->chr2;
            $gpid_strands{$gpid} = $gp->strand2;
        }

        my $pps = $ppa->fetch_by_gene_pair_id($gpid);
        my @search_regions;
        foreach my $pp (@$pps) {
            my ($tss_start1, $tss_end1) = define_search_region(
                $pp->tss1,    $pp->start1,  $pp->end1,
                $gp->strand1, $upstream_bp, $downstream_bp
            );

            my ($tss_start2, $tss_end2) = define_search_region(
                $pp->tss2,    $pp->start2,  $pp->end2,
                $gp->strand2, $upstream_bp, $downstream_bp
            );

            push @search_regions,
                {
                start1 => $tss_start1,
                end1   => $tss_end1,
                start2 => $tss_start2,
                end2   => $tss_end2,
                tss1   => $pp->tss1,
                tss2   => $pp->tss2
                };
        }
        $gpid_pps{$gpid} = $pps;
        $gpid_srs{$gpid} = \@search_regions;
    }

    foreach my $tf_id (@$tf_ids) {
        my $tf_name = $tf_info_set->get_tf_info($tf_id)->name;
        print HFH "$tf_name associated genes:\n";

        if ($id_type == 0) {
            #printf HFH "%-18s", 'Ensembl ID';
            print HFH "Ensembl ID";
        } else {
            #printf HFH "%-18s  %-18s", 'Gene ID', 'Ensembl ID';
            print HFH "Gene ID\tEnsembl ID";
        }
        #print HFH " Chr      Start        End Strand   Score Seq\n";
        print HFH "\tChr\tStart\tEnd\tStrand\tScore\tSeq\n";

        foreach my $gpid (@$t_gpids) {
            if ($t_ac->gene_tfbs_count($gpid, $tf_id) > 0) {
                my $sites = $ctfsa->fetch_by_gene_pair(
                    -gene_pair_id       => $gpid,
                    -tf_id              => $tf_id,
                    -conservation_level => $cons_level,
                    -threshold          => $threshold,
                    -upstream_bp        => $upstream_bp,
                    -downstream_bp      => $downstream_bp
                );

                my $first = 1;
                foreach my $site (@$sites) {
                    my $prev_start1 = 0;
                    my $prev_start2 = 0;
                    my $prev_end1   = 0;
                    my $prev_end2   = 0;
                    my $prev_site1  = 0;
                    my $prev_site2  = 0;
                    my $printed     = 0;
                    foreach my $region (@{$gpid_srs{$gpid}}) {

                        if ($species_num == 1) {
                            next
                                if ($region->{start1} == $prev_start1
                                && $region->{end1} == $prev_end1
                                && $site->start1 == $prev_site1);
                        } elsif ($species_num == 2) {
                            next
                                if ($region->{start2} == $prev_start2
                                && $region->{end2} == $prev_end2
                                && $site->start2 == $prev_site2);
                        }

                        if (   $site->start1 >= $region->{start1}
                            && $site->end1 <= $region->{end1})
                        {
                            if ($id_type == 0) {
                                if ($first) {
                                    #printf HFH "%-18s  %2s  %2d",
                                    printf HFH "%s\t%s\t%d\t",
                                        $gpid_ensids{$gpid},
                                        $gpid_chrs{$gpid},
                                        $gpid_strands{$gpid};
                                    $first = 0;
                                } else {
                                    #printf HFH "                          ";
                                    printf HFH "\t\t\t";
                                }
                            } else {
                                if ($first) {
                                    #printf HFH "%-15s  %-18s  %2s  %2d",
                                    printf HFH "%s\t%s\t%s\t%d\t",
                                        $gpid_gids{$gpid},
                                        $gpid_ensids{$gpid},
                                        $gpid_chrs{$gpid},
                                        $gpid_strands{$gpid};
                                    $first = 0;
                                } else {
                    #printf HFH "                                           ";
                                    printf HFH "\t\t\t\t";
                                }
                            }

                            if ($species_num == 1) {
                                my ($rel_start1, $rel_end1);
                                if ($gpid_strands{$gpid} == 1) {
                                    $rel_start1 =
                                        $site->start1 - $region->{tss1};
                                    if ($site->start1 >= $region->{tss1}) {
                                        $rel_start1++;
                                    }
                                    $rel_end1 = $site->end1 - $region->{tss1};
                                    if ($site->end1 >= $region->{tss1}) {
                                        $rel_end1++;
                                    }
                                } else {
                                    $rel_start1 =
                                        $region->{tss1} - $site->start1;
                                    if ($site->start1 <= $region->{tss1}) {
                                        $rel_start1++;
                                    }
                                    $rel_end1 = $region->{tss1} - $site->end1;
                                    if ($site->end1 <= $region->{tss1}) {
                                        $rel_end1++;
                                    }
                                }

                                if ($printed) {
                                    #print HFH '  *';
                                    print HFH '*';
                                } else {
                                    #print HFH '   ';
                                }

           #printf HFH "%9d  %9d  %9d  %25s  %9d  %5d  %9d  %5d  %2d  %.3g\n",
                                printf HFH
                                    "%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%.3g\n",
                                    $region->{tss1},
                                    $region->{start1},
                                    $region->{end1},
                                    $site->seq1,
                                    $site->start1,
                                    $rel_start1,
                                    $site->end1,
                                    $rel_end1,
                                    $site->strand1,
                                    $site->score1;

                                $prev_start1 = $region->{start1};
                                $prev_end1   = $region->{end1};
                                $prev_site1  = $site->start1;
                            } else {
                                my ($rel_start2, $rel_end2);
                                if ($gpid_strands{$gpid} == 1) {
                                    $rel_start2 =
                                        $site->start2 - $region->{tss2};
                                    if ($site->start2 >= $region->{tss2}) {
                                        $rel_start2++;
                                    }
                                    $rel_end2 = $site->end2 - $region->{tss2};
                                    if ($site->end2 >= $region->{tss2}) {
                                        $rel_end2++;
                                    }
                                } else {
                                    $rel_start2 =
                                        $region->{tss2} - $site->start2;
                                    if ($site->start2 <= $region->{tss2}) {
                                        $rel_start2++;
                                    }
                                    $rel_end2 = $region->{tss2} - $site->end2;
                                    if ($site->end2 <= $region->{tss2}) {
                                        $rel_end2++;
                                    }
                                }

           #printf HFH "%9d  %9d  %9d  %25s  %9d  %5d  %9d  %5d  %2d  %.3g\n",
                                printf HFH
                                    "%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%.3g\n",
                                    $region->{tss2},
                                    $region->{start2},
                                    $region->{end2},
                                    $site->seq2,
                                    $site->start2,
                                    $rel_start2,
                                    $site->end2,
                                    $rel_end2,
                                    $site->strand2,
                                    $site->score2;

                                $prev_start2 = $region->{start2};
                                $prev_end2   = $region->{end2};
                                $prev_site2  = $site->start2;
                            }
                            $printed = 1;
                        }
                    }
                }
            }
        }
        print HFH "\n";
    }
    close(HFH);
}

#
# Subroutine to determine the search region start/end based on the
# promoter and amount of upstream/downstream sequence
#
sub define_search_region
{
    my ($pp_tss, $pp_start, $pp_end, $strand, $upstream_bp, $downstream_bp) =
        @_;

    my $tss_start;
    my $tss_end;
    if ($strand == 1) {
        if (defined $upstream_bp) {
            $tss_start = $pp_tss - $upstream_bp;
            $tss_start = $pp_start if $pp_start > $tss_start;
        } else {
            $tss_start = $pp_start;
        }
        if (defined $downstream_bp) {
            $tss_end = $pp_tss + $downstream_bp - 1;
            $tss_end = $pp_end if $pp_end < $tss_end;
        } else {
            $tss_end = $pp_end;
        }
    } elsif ($strand == -1) {
        if (defined $upstream_bp) {
            $tss_end = $pp_tss + $upstream_bp;
            $tss_end = $pp_end if $pp_end < $tss_end;
        } else {
            $tss_end = $pp_end;
        }

        if (defined $downstream_bp) {
            $tss_start = $pp_tss - $downstream_bp + 1;
            $tss_start = $pp_start if $pp_start > $tss_start;
        } else {
            $tss_start = $pp_start;
        }
    } else {
        warn "Error determining gene pair strand";
        return;
    }

    return ($tss_start, $tss_end);
}

#
# Return a list of TF IDs in the order of the z-score results
#
sub sort_tf_ids_by_zscore_result
{
    my ($zscore_results) = @_;

    return if !$zscore_results;

    my $num_results = $zscore_results->num_results;

    my @tf_ids;
    my $result_idx = 0;
    while ($result_idx < $num_results) {
        my $zscore_result = $zscore_results->get_result($result_idx++);

        push @tf_ids, $zscore_result->id;
    }

    return @tf_ids ? \@tf_ids : undef;
}
