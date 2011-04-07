#!/usr/local/bin/perl -w

#
# quick & dirty script to look for those tfbs cluster hit sum_lengths
# that add up to more than the conserved region length.
# Probably due to the 'creeping' clusters: TFBS that sits at the edge 
# of the conserved regions, and keep extending out of bounds
# when added up, can be significant
#


use strict;

use lib '/space/devel/oPOSSUM3/cgi-bin';
use oPossumWebOpt;
use oPossumGeneWebOpt;

use lib OPOSSUM_LIB_PATH;

use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor;

my $db_name = 'oPOSSUM_2010_worm';
my $clevel = 3;
my $slevel = 3;
my $tlevel = 2;
my $outfile;
GetOptions(
	-db		=> \$db_name,
	-cl		=> \$clevel,
	-srl	=> \$slevel,
	-tl		=> \$tlevel,
	-o		=> \$outfile
);

if (!$outfile) {
	die "Must provide outfile name\n";
}

my $dba = OPOSSUM::DBSQL::DBAdaptor->new(
    -host     => OPOSSUM_DB_HOST,
    -dbname   => $db_name,
    -user     => OPOSSUM_DB_USER,
    -password => OPOSSUM_DB_PASS
);

die "Can't connect to $db_name\n" if !$dba;

my $ga = $dba->get_GeneAdaptor;
die "Can't get gene adaptor\n" if !$ga;
my $cra = $dba->get_ConservedRegionAdaptor;
die "Can't get conserved region adaptor\n" if !$cra;
my $aca = $dba->get_AnalysisClusterCountsAdaptor;
die "Can't get CountsAdaptor\n" if !$aca;
my $srla = $dba->get_SearchRegionLevelAdaptor;
die "Can't get search region level adaptor\n" if !$srla;
my $tcca = $dba->get_TFBSClusterCountAdaptor;
die "Can't get tfbs cluster counts adaptor\n" if !$tcca;

my $srl_levels = $srla->fetch_levels();
die("getting search region levels\n") if !$srl_levels;

my $srl_hash = $srla->fetch_search_region_level_hash();
die("getting search region levels\n") if !$srl_hash;

open (OUTPUT, ">$outfile") or die "Can't open $outfile\n";

my $cids = $tcca->fetch_cluster_ids;
my $gids = $ga->fetch_gene_ids;

foreach my $gid (@$gids)
{
    my $gene = $ga->fetch_by_id($gid);
    die "fetching gene $gid\n" if !$gene;

    my $srl = $srl_hash->{$slevel};

    my $upstream_bp   = $srl->upstream_bp();
    my $downstream_bp = $srl->downstream_bp();

    my $length = $cra->fetch_length_by_upstream_downstream(
        $gid, $clevel, $upstream_bp, $downstream_bp
    );
    
    my $counts = $tcca->fetch_gene_tfbs_cluster_counts(
        -gene_id => $gid,
        #-cluster_ids => $cids,
        -conservation_level => $clevel,
        -threshold_level => $tlevel,
        -search_region_level => $slevel
    );
    
    foreach my $count (@$counts) {
        if ($count->sum_length > $length) {
            print OUTPUT "Gene $gid Cluster " . $count->cluster_id
            . ": crl = $length\tsuml = "
            . $count->sum_length . "\n";
        }
    }

}

