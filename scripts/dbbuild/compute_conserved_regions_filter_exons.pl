#!/usr/bin/perl -w

=head1 NAME

compute_conserved_regions_filter_exons.pl

=head1 SYNOPSIS

  compute_conserved_regions_filter_exons.pl -d opossum_db_name -h opossum_db_host
    -ud ucsc_db_name -ut ucsc_track_name [-s start] [-e end] [-fe filter_all_exons]
    -o cons_regs_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -d opossum_db_name   = Name of oPOSSUM db we are working on
   -h opossum_db_host   = Host name of the ORCA db
   -s start             = Start computing with this gene id
   -e end               = Finish computing with this gene id
   -fe filter_all_exons	= filter all exons within the promoter search region
   -ud ucsc_db_name     = UCSC DB name
   -ut ucsc_track_name  = UCSC track name
   -o cons_regs_file    = Ouput conserved regions file (for import into
                          oPOSSUM conserved_regions table using mysqlimport)
   -l log_file          = Name of log file to which processing and error
                          messages are written.
                          (Default = compute_conserved_regions.log)

=head1 DESCRIPTION

This is the oPOSSUM_2010 script for computing conserved regions at each
conservation level using phastCons conservation.

=head1 ALGORITHM

Foreach gene in the oPOSSUM database, read the corresponding gene, exon
and sequence information and use the ORCA::Analysis::PhastCons module to
compute the conserved regions at each of the conservation levels. The
conserved regions are filtered by promoter search regions based around
the maximum upstream and downstream sequence (search region level 1).
Write out the conserved regions to the specified output file.

=head1 MODIFICATIONS

Andrew Kwon
- added filter exon option, to exclude exons from conserved regions
- added GC content calculation

=head1 AUTHOR

  David Arenillas, Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca, tjkwon@cmmt.ubc.ca

=cut

use strict;

#
# Use most current (development) libs.
# Comment out to use installed libs.
#

# on Watson / cluster
#use lib '/home/dave/devel/oPOSSUM_2010/lib';
#use lib "/raid2/tjkwon/oPOSSUM_2010/lib";
#use lib '/home/dave/devel/ORCAtk/lib';
#use lib '/raid2/local/src/ensembl-57/ensembl/modules';

# on burgundy
use lib "/space/devel/oPOSSUM3/lib";
use lib '/usr/local/src/ensembl-57/ensembl/modules';
use lib '/space/devel/ORCAtk/lib';

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use Bio::LocatableSeq;
use OPOSSUM::DBSQL::DBAdaptor;
use ORCA::Analysis::PhastCons;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use constant DEBUG              => 1;

use constant LOG_FILE 		    => 'compute_conserved_regions.log';

use constant OPOSSUM_DB_USER	=> "opossum_r";

use constant MIN_CONSERVED_REGION_LENGTH    => 20;

# Local Ensembl database settings
use constant ENSEMBL_DB_HOST    	=> 'vm2.cmmt.ubc.ca';
use constant ENSEMBL_DB_USER    	=> 'ensembl_r';
use constant ENSEMBL_DB_PASS    	=> '';

use constant ENSEMBL_DB_NAME_HUMAN	=> 'homo_sapiens_core_57_37b';
use constant ENSEMBL_DB_NAME_MOUSE	=> 'mus_musculus_core_57_37j';
use constant ENSEMBL_DB_NAME_FLY	=> 'drosophila_melanogaster_core_57_513b';
use constant ENSEMBL_DB_NAME_WORM	=> 'caenorhabditis_elegans_core_57_200a';
use constant ENSEMBL_DB_NAME_YEAST	=> 'saccaromyces_cerevisiae_core_57_1j';

use constant HUMAN => 'human';
use constant MOUSE => 'mouse';
use constant FLY => 'fly';
use constant FRUITFLY => 'fruitfly';
use constant WORM => 'worm';
use constant NEMATODE => 'nematode';
use constant YEAST => 'yeast';

my $log_file = LOG_FILE;
my $opossum_db_name;
my $opossum_db_host;
my $ucsc_db_name;
my $ucsc_track_name;
my $start_gid;
my $end_gid;
my $filter_all_exons;
my $out_regions_file;
my $out_tfbs_file;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
    'ud=s'  => \$ucsc_db_name,
    'ut=s'  => \$ucsc_track_name,
    's=i'	=> \$start_gid,
    'e=i'	=> \$end_gid,
	'fe=i'	=> \$filter_all_exons,
    'o=s'	=> \$out_regions_file,
    't=s'	=> \$out_tfbs_file,
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

if (!$ucsc_db_name) {
    pod2usage(
        -msg        => "No UCSC DB name specified",
        -verbose    => 1
    );
}

if (!$ucsc_track_name) {
    pod2usage(
        -msg        => "No UCSC track name specified",
        -verbose    => 1
    );
}

if (!$out_regions_file) {
    pod2usage(
        -msg        => "No output conserved regions file specified",
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


my $opdba = OPOSSUM::DBSQL::DBAdaptor->new(
    -host            => $opossum_db_host,
    -dbname          => $opossum_db_name,
    -user            => OPOSSUM_DB_USER,
    -password        => undef
);

if (!$opdba) {
    $logger->logdie($DBI::errstr);
}

my $dbia = $opdba->get_DBInfoAdaptor;
if (!$dbia) {
    $logger->logdie("getting DBInfoAdaptor");
}

my $db_info = $dbia->fetch_db_info();
if (!$db_info) {
    $logger->logdie("fetching DB info");
}

my $ga = $opdba->get_GeneAdaptor;
if (!$ga) {
    $logger->logdie("getting GeneAdaptor");
}

my $exa = $opdba->get_ExonAdaptor;
if (!$exa) {
    $logger->logdie("getting ExonAdaptor");
}

my $cla = $opdba->get_ConservationLevelAdaptor;
if (!$cla) {
    $logger->logdie("getting ConservationLevelAdaptor");
}

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

my $srla = $opdba->get_SearchRegionLevelAdaptor;
if (!$srla) {
    $logger->logdie("getting SearchRegionLevelAdaptor");
}

my $srl_levels = $srla->fetch_levels;

my $srl_hash = $srla->fetch_search_region_level_hash;
if (!$srl_hash) {
    $logger->logdie("fetching search region levels");
}

$logger->info("Search region levels:");
foreach my $level (@$srl_levels) {
    my $srl = $srl_hash->{$level};
    $logger->info(
        sprintf("\t%d\t%d\t%d\n",
            $srl->level,
            $srl->upstream_bp,
            $srl->downstream_bp
        )
    );
}

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

my $species = $db_info->species();

my $max_upstream_bp = $db_info->max_upstream_bp;
my $max_downstream_bp = $max_upstream_bp;

#
# Special case for yeast, search region goes to 3' end of gene. Setting
# max_downstream_bp to undef causes the oPOSSUM::Gene::promoter_search_regions
# method to search to end of gene.
#
if ($species eq 'yeast') {
    $max_downstream_bp = undef;
}

# select the right ENSEMBL DB name
my $ens_db_name;
if ($species eq HUMAN) {
	$ens_db_name = ENSEMBL_DB_NAME_HUMAN;
} elsif ($species eq MOUSE) {
	$ens_db_name = ENSEMBL_DB_NAME_MOUSE;
} elsif ($species eq FLY or $species eq FRUITFLY) {
	$ens_db_name = ENSEMBL_DB_NAME_FLY;
} elsif ($species eq WORM or $species eq NEMATODE) {
	$ens_db_name = ENSEMBL_DB_NAME_WORM;
} elsif ($species eq YEAST) {
	$ens_db_name = ENSEMBL_DB_NAME_YEAST;
} else {
	$logger->logdie("Unknown species name $species. "
					. "Options are: human, mouse, fly, worm and yeast.");
}

# connect to ENSEMBL DB
my $ens_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host	    => ENSEMBL_DB_HOST,
    -user	    => ENSEMBL_DB_USER,
    -pass	    => ENSEMBL_DB_PASS,
    -dbname	    => $ens_db_name,
    -species	=> $species,
    -driver	    => 'mysql'
);

if (!$ens_db) {
    $logger->logdie("connecting to Ensembl $species core DB $ens_db_name");
}

my $ens_ga = $ens_db->get_GeneAdaptor || $logger->logdie("getting GeneAdaptor");
my $ens_sa = $ens_db->get_SliceAdaptor || $logger->logdie("getting SliceAdaptor");

open(OCRFH, ">$out_regions_file")
    || $logger->logdie("opening output conserved regions file");

foreach my $gid (@$gids) {
    $logger->info("processing gene $gid");

    my $opgene = $ga->fetch_by_gene_id($gid);
    if (!$opgene) {
    	$logger->logdie("fetching gene");
    }
	
	# opossum gene starts and ends have upstream and downstream offset built in
	# - for downstream, up to 3' end of the gene only
	my $sr_start = $opgene->start;
	my $sr_end = $opgene->end;
	my $slice = $ens_sa->fetch_by_region('chromosome', $opgene->chr, $sr_start, $sr_end);
	my $search_seq = $slice->seq;
    
	#my $opseq = $seqa->fetch_by_gene_id($gid);
    #if (!$opseq) {
    #	$logger->logdie("fetching sequence");
    #}

    #my $biolseq = Bio::LocatableSeq->new(
    #    -primary_id => $gid,
    #    -display_id => sprintf("chr%s:%d-%d",
    #                        $opgene->chr, $opgene->start, $opgene->end),
    #    -start      => $opgene->start,
    #    -end        => $opgene->end,
    #    -seq        => $opgene->sequence->masked_seq
    #);

    my $search_regs = $opgene->promoter_search_regions(
        $max_upstream_bp, $max_downstream_bp
    );

    #
    # Convert search regions to sequence (1-based) coords
    # Changed OPOSSUM::Gene::PromoterSearchRegions routine to return
    # sequence based coords.
    #
    #foreach my $sr (@$search_regs) {
    #    $sr->start($sr->start - $opgene->start + 1);
    #    $sr->end($sr->end - $opgene->start + 1);
    #}

    #
    # For yeast filter out all exons (including upstream gene exons) from
    # conserved regions. For all other species only filter exons of this
    # gene.
    #
    # NOTE: in the yeast case this only includes exons from genes actually in
    # the oPOSSUM DB (maybe not all genes) but oPOSSUM should be including
    # all 'KNOWN' genes with assigned symbols so this is probably OK.
    #
	# Andrew: change this so that you can do so for any species
	#
    my $exons;
    if ($filter_all_exons) {
        $exons = $exa->fetch_by_region(
            $opgene->chr, $opgene->start, $opgene->end
        );
    } else {
        $exons = $opgene->exons();
    }

    #
    # XXX this could be made more efficient. We don't need to retrieve
    # phastCons scores for the whole gene region, just the portion within
    # the maximal search region.
    #
    my $phca = ORCA::Analysis::PhastCons->new(
        -db             => $ucsc_db_name,
        -track          => $ucsc_track_name,
        #-seq            => $biolseq,
        -chr            => $opgene->chr,
        -start          => $opgene->start,
        -end            => $opgene->end,
        -exons          => $exons
    );

    foreach my $clevel (@$cl_levels) {
        my $cl = $cl_hash->{$clevel};

        my $conservation_level  = $cl->level;
        my $min_conservation    = $cl->min_conservation;

        my $cons_regs = $phca->compute_conserved_regions(
            -min_conservation   => $min_conservation,
            -filter_exons       => 1,
            -min_length         => MIN_CONSERVED_REGION_LENGTH,
        );

        if (!$cons_regs || !$cons_regs->[0]) {
            if ($clevel == 1) {
                $logger->warn("No conserved regions");
            } else {
                $logger->info("No level $clevel conserved regions");
            }
            last;
        }

        $cons_regs = filter_conserved_regions_by_search_regions(
            $cons_regs, $search_regs
        );

        if (!$cons_regs || !$cons_regs->[0]) {
            if ($clevel == 1) {
                $logger->warn(
                    "No conserved regions within promoter search regions"
                );
            } else {
                $logger->info("No level $clevel conserved regions within"
                    . " promoter search regions"
                );
            }
            last;
        }
        
		write_conserved_regions(\*OCRFH, $opgene->id, $clevel, $cons_regs, $search_seq);
    }
}
close(OCRFH);

my $end_time = time;
$localtime = localtime($end_time);
my $elapsed_secs = $end_time - $start_time;

$logger->info("oPOSSUM update completed on $localtime");
$logger->info("Elapsed time (s): $elapsed_secs");

exit;

#
# Return only those conserved regions which are contained within or which
# overlap the promoter search regions.
# 
# DO NOT truncate conserved regions at search region boundaries.
#
sub filter_conserved_regions_by_search_regions
{
    my ($cons_regs, $search_regs) = @_;

    my @filtered_regs;

    foreach my $cr (@$cons_regs) {
        foreach my $sr (@$search_regs) {
            if ($cr->end >= $sr->start && $cr->start <= $sr->end) {
                #
                # XXX Should we truncate conserved regions which overlap the
                # search region at the start/end of the search regions???
                #
                push @filtered_regs, $cr;
                last;
            }
        }
    }

    return @filtered_regs ? \@filtered_regs : undef;
}

sub write_conserved_regions
{
    my ($fh, $gid, $level, $crs, $search_seq) = @_;

    foreach my $cr (@$crs) {
		
		#
		# calculate the GC content
		#
		my $reg_length = $cr->end - $cr->start + 1;
		my $cr_start = $cr->start;
		my $cr_end = $cr->end;
		my $cr_score = $cr->score;
		my $cr_seq = substr($search_seq, $cr_start, $reg_length);
		
		$logger->debug("cl level: $level\t$cr_start-$cr_end");
		$logger->debug("$cr_seq");
		
		my $gc_content = calculate_gc_content($cr_seq);
		$logger->debug("GC: $gc_content");
		
	    printf $fh "%d\t%d\t%d\t%d\t%.5f\t%.5f\n",
            $gid, $level, $cr->start, $cr->end, $cr->score, $gc_content;
    }
}

sub calculate_gc_content
{
	my ($seq) = @_;
	
	my %count = (
		'A' => 0,
		'C' => 0,
		'G' => 0,
		'T' => 0
	);
	
	for (my $i = 0; $i < length($seq); $i++)
	{
		my $nt = substr($seq, $i, 1);
		$count{$nt}++;
	}
	
	my $gc = $count{'G'} + $count{'C'};
	my $gc_content = $gc / length($seq);
	
	return $gc_content;
}
