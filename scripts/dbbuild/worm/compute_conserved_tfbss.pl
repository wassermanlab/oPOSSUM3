#!/usr/bin/perl -w

=head1 NAME

compute_conserved_tfbss.pl

=head1 SYNOPSIS

  compute_conserved_tfbss.pl -d opossum_db_name -h opossum_db_host
    [-s start] [-e end] -o out_tfbs_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -d opossum_db_name   = Name of ORCA db we are working on.
   -h opossum_db_host   = Host name of the ORCA db.
   -s start             = Start computing with this gene id.
   -e end               = Finish computing with this gene id.
   -o out_tfbs_file     = Ouput conserved TFBSs file (for import into
                          oPOSSUM conserved_tfbss table using mysqlimport).
   -l log_file          = Name of log file to which processing and error
                          messages are written.
                          (Default = compute_conserved_tfbss.log)

=head1 DESCRIPTION

This is the oPOSSUM_2010 script for computing conserved TFBSs. The genes,
exons(?) and conserved regions tables must already be populated.

=head1 ALGORITHM

Foreach gene in the oPOSSUM database, use the
OPOSSUM::DBSQL::Analysis::PhastConsAdaptor module to read in an
ORCA::Analysis::PhastCons object. Use this object to compute the conserved
TFBSs and write them out to a file. The TFBSs are computed once for the
lowest level of conservation and then assigned the conservation score of
the highest scoring conserved region which they fall into.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

# use most current (development) libs
# comment out to use installed libs
use lib '/home/tjkwon/OPOSSUM/oPOSSUM3/lib';
use lib '/home/dave/devel/ORCAtk/lib';
use lib '/home/tjkwon/ensembl/ensembl-54/modules';

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use OPOSSUM::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use TFBS::DB::JASPAR5;
use ORCA::Analysis::PhastCons;


use constant DEBUG                      => 0;

use constant LOG_FILE 		            => 'compute_conserved_tfbss.log';

use constant OPOSSUM_DB_USER	        => "opossum_r";

use constant JASPAR_DB_HOST             => 'vm5.cmmt.ubc.ca';
use constant JASPAR_DB_NAME             => 'JASPAR_2010';
use constant JASPAR_DB_USER             => 'jaspar_r';
use constant JASPAR_DB_PASS             => '';

use constant ENSEMBL_DB_HOST            => 'vm2.cmmt.ubc.ca';
use constant ENSEMBL_DB_USER            => 'ensembl_r';
use constant ENSEMBL_DB_PASS            => '';

# Just use CORE/vertebrate matrices for now
use constant JASPAR_COLLECTIONS         => ['CORE'];
use constant JASPAR_TAX_GROUPS          => ['vertebrates'];

use constant MIN_TFBS_CR_OVERLAP        => 1;
use constant FILTER_OVERLAPPING_TFBSS   => 1;


my $log_file = LOG_FILE;
my $opossum_db_name;
my $opossum_db_host;
my $start_gid;
my $end_gid;
my $out_tfbs_file;
my $tax_group_str;
my $collection_str;
my $min_ic;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
    's=i'	=> \$start_gid,
    'e=i'	=> \$end_gid,
    'ic=s'	=> \$min_ic,
    'c=s'	=> \$collection_str,
    't=s'	=> \$tax_group_str,
    'o=s'	=> \$out_tfbs_file,
    'l=s'	=> \$log_file
);

if (!$opossum_db_name) {
    pod2usage(
        -msg        => "No oPOSSUM DB name specified",
        -verbose    => 1
    );
}

if (!$opossum_db_host) {
    pod2usage(
        -msg        => "No oPOSSUM DB host specified",
        -verbose    => 1
    );
}

if (!$out_tfbs_file) {
    pod2usage(
        -msg        => "No output conserved TFBSs file specified",
        -verbose    => 1
    );
}

#
# Initialize logging
#
my $logger = get_logger();
if (DEBUG) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}
my $appender = Log::Log4perl::Appender->new("Log::Dispatch::File",
                                            filename    => $log_file,
                                            mode        => "write");
#my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("%M:%L %p: %m%n");
$appender->layout($layout);
$logger->add_appender($appender);

my $start_time = time;
my $localtime = localtime($start_time);

$logger->info("compute_conserved_tfbss started on $localtime\n");

my $jdbh = TFBS::DB::JASPAR5->connect(
    "dbi:mysql:" . JASPAR_DB_NAME . ":" . JASPAR_DB_HOST,
    JASPAR_DB_USER,
    JASPAR_DB_PASS
);

if (!$jdbh) {
    $logger->logdie("connecting to JASPAR database - $DBI::errstr");
}

my $opdba = OPOSSUM::DBSQL::DBAdaptor->new(
    -host       => $opossum_db_host,
    -dbname     => $opossum_db_name,
    -user       => OPOSSUM_DB_USER,
    -password   => undef
);

if (!$opdba) {
    $logger->logdie("connecting to oPOSSUM database - $DBI::errstr");
}

#
# Get some adaptors up front
#
my $dbia = $opdba->get_DBInfoAdaptor;
if (!$dbia) {
    $logger->logdie("getting DBInfoAdaptor");
}

my $ga = $opdba->get_GeneAdaptor;
if (!$ga) {
    $logger->logdie("getting GeneAdaptor");
}

#my $phcaa = $opdba->get_PhastConsAdaptor;
#if (!$phcaa) {
#    $logger->logdie("getting PhastConsAdaptor");
#}

#my $tfia = $opdba->get_TFInfoAdaptor;
#if (!$tfia) {
#    $logger->logdie("getting TFInfoAdaptor");
#}

my $cra = $opdba->get_ConservedRegionAdaptor;
if (!$cra) {
    $logger->logdie("getting ConservedRegionAdaptor");
}

my $cla = $opdba->get_ConservationLevelAdaptor;
if (!$cla) {
    $logger->logdie("getting ConservationLevelAdaptor");
}

my $thla = $opdba->get_ThresholdLevelAdaptor;
if (!$thla) {
    $logger->logdie("getting ThresholdLevelAdaptor");
}

my $srla = $opdba->get_SearchRegionLevelAdaptor;
if (!$srla) {
    $logger->logdie("getting SearchRegionLevelAdaptor");
}


my $db_info = $dbia->fetch_db_info;
if (!$db_info) {
    $logger->logdie("fetching DB info");
}

my $ens_db_name = $db_info->ensembl_db()
    || $logger->logdie("Ensembl DB name not set in db_info table");

#my $ens_ver;
#if ($ens_db_name =~ /_core_(\d+)_/) {
#    $ens_ver = $1;
#}
#
#$logger->logdie("determining Ensembl version") if !$ens_ver;
#
#my $ens_lib = "/usr/local/src/ensembl-${ens_ver}/ensembl/modules";
#
#use lib $ens_lib;

my $ensdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => ENSEMBL_DB_HOST,
    -user    => ENSEMBL_DB_USER,
    -pass    => ENSEMBL_DB_PASS,
    -dbname  => $ens_db_name,
    -species => $db_info->species(),
    -driver  => 'mysql'
);

if (!$ensdba) {
    $logger->logdie("connecting to Ensembl DB");
}

my $enssa = $ensdba->get_SliceAdaptor()
    || $logger->logdie("getting Ensembl SliceAdaptor");


my $cl_levels = $cla->fetch_levels;

my $cl_hash = $cla->fetch_conservation_level_hash;
if (!$cl_hash) {
    $logger->logdie("fetching conservation levels");
}

$logger->info("Conservation levels:");
foreach my $level (@$cl_levels) {
    my $cl = $cl_hash->{$level};
    $logger->info(
        sprintf("\t%d\t%0.2f\n",
            $cl->level,
            $cl->min_conservation
        )
    );
}

my $thresh_level1 = $thla->fetch_by_level(1);
if (!$thresh_level1) {
    $logger->logdie("fetching threshold level 1");
}

my $min_threshold = $thresh_level1->threshold();
$min_threshold = $min_threshold * 100 . '%';

$logger->info("min. TFBS threshold: $min_threshold");

my $sr_level1 = $srla->fetch_by_level(1);
if (!$sr_level1) {
    $logger->logdie("fetching search region level 1");
}

my $max_upstream_bp     = $sr_level1->upstream_bp();
my $max_downstream_bp   = $sr_level1->downstream_bp();

$logger->info(
    "max. upstream/downstream bp: $max_upstream_bp/$max_downstream_bp"
);

my $where;
if ($start_gid) {
    $where = "gene_id >= $start_gid";
    if ($end_gid) {
    	$where .= " and gene_id <= $end_gid";
    }
} elsif ($end_gid) {
    $where = "gene_id <= $end_gid";
}

my $gids = $ga->fetch_gene_ids($where);
if (!$gids) {
    $logger->logdie("fetching gene IDs");
}

my $collections = JASPAR_COLLECTIONS;
if ($collection_str) {
	my @collections = split ",", $collection_str;
	$collections = \@collections;
}

my $tax_groups = JASPAR_TAX_GROUPS;
if ($tax_group_str) {
	my @tax_groups = split ",", $tax_group_str;
	$tax_groups = \@tax_groups;
}

my $matrix_set = $jdbh->get_MatrixSet(
    -collection => $collections,
    -tax_group  => $tax_groups,
    -matrixtype => 'PWM',
    -min_ic	=> $min_ic
);

if (!$matrix_set) {
    $logger->logdie("fetching JASPAR matrix set");
}

#my $tf_id_map = fetch_tf_id_mapping();

open(OTFH, ">$out_tfbs_file") || $logger->logdie("opening output TFBS file");

my $min_conservation_level = 1;
my %cons_lvl_cons_regs;
foreach my $gid (@$gids) {
    $logger->info("processing gene $gid");

    my $gene = $ga->fetch_by_id($gid);
    if (!$gene) {
        $logger->logdie("fetching gene");
    }

    foreach my $level (@$cl_levels) {
        if (!$gene->fetch_conserved_regions($level)) {
            if ($level == $min_conservation_level) {
                $logger->warn("no conserved regions");
                next;
            }
        }
    }

    #if (!$gene->sequence()) {
    #    $logger->error("no gene sequence");
    #    next;
    #}

    #my $seq  = $gene->sequence()->seq();
    #my $mseq = $gene->sequence()->masked_seq();

    my $slice = $enssa->fetch_by_region(
        "chromosome", $gene->chr(), $gene->start(), $gene->end()
    );

    if (!$slice) {
        $logger->logdie(
            sprintf(
                "fetching slice chr%s:%d-%d from Ensembl", 
                $gene->chr(), $gene->start(), $gene->end()
            )
        );
    }

    my $biolseq = Bio::LocatableSeq->new(
        -id         => sprintf("chr%s:%d-%d",
                            $gene->chr, $gene->start, $gene->end),
        -start      => $gene->start,
        -end        => $gene->end,
        -strand     => 1,
        #-seq        => $mseq
        -seq        => $slice->get_repeatmasked_seq()->seq()
    );

    my $phca = ORCA::Analysis::PhastCons->new(
        -seq            => $biolseq,
        -chr            => $gene->chr,
        -start          => $gene->start,
        -end            => $gene->end,
        -exons          => $gene->exons
    );

    $phca->conserved_regions($gene->conserved_regions($min_conservation_level));

    my $search_regs = $gene->promoter_search_regions(
        $max_upstream_bp, $max_downstream_bp
    );

    #
    # Convert search regions to sequence based coords
    # XXX this is already done by the OPOSSUM::Gene::promoter_search_regions()
    # method
    #
    #foreach my $sr (@$search_regs) {
    #    $sr->start($sr->start - $gene->start + 1);
    #    $sr->end($sr->end - $gene->start + 1);
    #}

    my @cons_tfbss;
    foreach my $sr (@$search_regs) {
        my $tfbss = $phca->compute_conserved_tfbss(
            -matrix_set                 => $matrix_set,
            -min_tfbs_score             => $min_threshold,
            -min_tfbs_cr_overlap        => MIN_TFBS_CR_OVERLAP,
            -filter_overlapping_tfbss   => FILTER_OVERLAPPING_TFBSS,
            -start                      => $sr->start,
            -end                        => $sr->end
        );

        push @cons_tfbss, @$tfbss if $tfbss;
    }

    if (!@cons_tfbss) {
        $logger->warn("no conserved TFBSs found");
        next;
    }

    #@cons_tfbss = @{filter_conserved_tfbss_by_search_regions(
    #    \@cons_tfbss, $search_regs
    #)};
    #
    #if (!$cons_tfbss) {
    #    $logger->warn("no conserved TFBSs withing promoter search regions");
    #    next;
    #}

    # Convert the TFBS::Site list to an OPOSSUM::ConservedTFBS list
    #my $op_tfbss = sites_to_opossum_tfbss(
    #    $cons_tfbss, $gid, $tf_id_map
    #);
    #
    # For efficiency, don't convert TFBS::Sites to OPOSSUM:ConservedTFBSs.
    # Just add tag values for gene ID and TF ID.
    foreach my $tfbs (@cons_tfbss) {
        $tfbs->add_tag_value('gene_id', $gid);

        my $ext_tf_id = $tfbs->pattern()->ID();

        #my $tf_id;
        #if ($tf_id_map->{$ext_tf_id}) {
        #    $tf_id = $tf_id_map->{$ext_tf_id}->id();
        #}

        #$logger->logdie(
        #    "determining oPOSSUM TF ID for external TF ID $ext_tf_id"
        #) if !$tf_id;

        #$tfbs->add_tag_value('tf_id', $tf_id);
    }

    tfbss_set_max_conservation($gene, $cl_levels, \@cons_tfbss);

    write_conserved_tfbss(\*OTFH, $gid, \@cons_tfbss);
}
close(OTFH);

my $end_time = time;
$localtime = localtime($end_time);
my $elapsed_secs = $end_time - $start_time;

$logger->info("compute_conserved_tfbss completed on $localtime");
$logger->info("Elapsed time (s): $elapsed_secs");

exit;

#
# Return only those conserved TFBSs which are contained completely within
# the promoter search regions.
# 
sub filter_conserved_tfbss_by_search_regions
{
    my ($tfbss, $search_regs) = @_;

    my @filtered_tfbss;

    foreach my $tfbs (@$tfbss) {
        foreach my $sr (@$search_regs) {
            if ($tfbs->start >= $sr->start && $tfbs->end <= $sr->end) {
                push @filtered_tfbss, $tfbs;
                last;
            }
        }
    }

    return @filtered_tfbss ? \@filtered_tfbss : undef;
}

=head3
sub fetch_tf_id_mapping
{
    my $tf_info_set = $tfia->fetch_tf_info_set();

    my @tf_ids = $tf_info_set->tf_ids;

    my %map;
    foreach my $id (@tf_ids) {
        my $tf_info = $tf_info_set->get_tf_info($id);

        $map{$tf_info->external_id} = $tf_info;
    }

    return %map ? \%map : undef;
}
=cut

sub write_conserved_tfbss
{
    my ($fh, $gid, $tfbss) = @_;

    foreach my $tfbs (@$tfbss) {
        printf $fh "%d\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%s\t%d\t%.3f\n",
            $gid,
            # XXX need map of TF ID
            #($tfbs->get_tag_values('tf_id'))[0],
            $tfbs->pattern->ID,
            $tfbs->start,
            $tfbs->end,
            $tfbs->strand,
            $tfbs->score,
            $tfbs->rel_score,
            $tfbs->seq->seq,
            # XXX need map of TF conservation level and conservation
            ($tfbs->get_tag_values('conservation_level'))[0],
            ($tfbs->get_tag_values('conservation'))[0]
    }
}

#
# Convert a list of TFBS::Site objects to a list of OPOSSUM::ConservedTFBS
# objects
#
sub sites_to_opossum_tfbss
{
    my ($sites, $gid, $tf_id_map) = @_;

    return if !$sites || !$sites->[0];

    my @tfbss;
    foreach my $site (@$sites) {
        my $tf_id;
        if ($tf_id_map) {
            $tf_id = $tf_id_map->{$site->pattern->ID}->id;
        } else {
            $tf_id = $site->pattern->ID || $site->pattern->name;
        }

        if (!$tf_id) {
            $logger->logdie("determining TF ID");
        }

        push @tfbss, OPOSSUM::ConservedTFBS->new(
			-gene_pair_id   => $gid,
			-tf_id          => $tf_id,
			-start          => $site->start,
			-end            => $site->end,
			-strand         => $site->strand,
			-score          => $site->score,
			-rel_score      => $site->rel_score,
			-seq            => $site->seq->seq
        );
    }

    return @tfbss ? \@tfbss : undef;
}

sub tfbss_set_max_conservation
{
    my ($gene, $levels, $tfbss) = @_;

    # Reverse sort so we can check conservation from highest to lowest
    my @rlevels = reverse @$levels;

    foreach my $tfbs (@$tfbss) {
        my $max_level = 0;
        my $max_conservation = 0;
        foreach my $level (@rlevels) {
            my $crs = $gene->conserved_regions($level);
            next if !$crs;

            foreach my $cr (@$crs) {
                #
                # Conserved regions are sorted so speed up search by exiting
                # loop early if we have already searched past end of TFBS
                #
                last if $cr->start > $tfbs->end - MIN_TFBS_CR_OVERLAP + 1;
                next if $cr->end < $tfbs->start + MIN_TFBS_CR_OVERLAP - 1;

                #
                # A TFBS could overlap more than one conserved region at
                # the same conservation level, so explicitly check if
                # the current conservation score is greater than any
                # previously set one.
                #
                if ($cr->conservation > $max_conservation) {
                    $max_level = $level;
                    $max_conservation = $cr->conservation;
                }
            }

            #
            # Exit loop early if TFBS overlaps a conserved region at this
            # level (any remaining conservation levels are at a lower level
            # of conservation than this one
            #
            last if $max_level;
        }

        if ($max_level) {
            $tfbs->add_tag_value('conservation_level', $max_level);

            #
            # Must remove conservation tag that was added during the phastCons
            # analysis and then re-add it, otherwise conservation score is just
            # pushed onto the existing conservation tag array
            #
            $tfbs->remove_tag('conservation');
            $tfbs->add_tag_value('conservation', $max_conservation);
        } else {
            $logger->error("computing max. conservation for "
                        . $tfbs->pattern()->name());
        }
    }
}
