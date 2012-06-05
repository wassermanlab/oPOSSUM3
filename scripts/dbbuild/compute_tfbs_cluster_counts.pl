#!/usr/bin/perl -w

=head1 NAME

compute_tfbs_cluster_counts.pl

=head1 SYNOPSIS

  compute_tfbs_cluster_counts.pl
        -h opossum_db_host -d opossum_db_name
        -cd cluster_db_name -ch cluster_db_host[-s start] [-e end]
        [-id cluster_id] -o out_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

  -h opossum_db_host  = Host name of oPOSSUM DB to process
  -d opossum_db_name  = Name of oPOSSUM DB to process
  -cd cluster_db_name = Name of TFBS_cluster DB to process
  -ch cluster_db_host = Host name of TFBS_cluster DB to process
  -s start            = Starting gene ID
  -e end              = Ending gene ID
  -cid cluster_id     = Specific cluster ID to be counted. If this is
                        specified, only this cluster hits will be counted.
  -o out_file         = Output TFBS cluster counts file
  -l log_file         = Name of log file to which processing and error
                        messages are written
                        (Default = compute_tfbs_cluster_counts.log)

=head1 DESCRIPTION

For each gene, for each TF cluster, at each level of conservation, promoter search
region and matrix score it counts the number of TFBS clusters from the conserved_tfbss
table and write out these counts to a file for importing into the tfbs_cluster_counts
table.

=head1 AUTHOR

  Andrew Kwon, adapting code by David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

#use lib "/space/devel/oPOSSUM3/lib";
#use lib "/space/devel/TFBS_Cluster/lib";
#use lib "/raid2/tjkwon/oPOSSUM_2010/lib";
#use lib "/raid2/tjkwon/TFBS_Cluster/lib";
#use lib '/home/dave/devel/TFBS';
use lib "/apps/oPOSSUM3/lib";

use strict;

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::ConservedTFBS;
use TFBSCluster::DBSQL::DBAdaptor;

use constant DEBUG     => 1;
use constant UPDATE_DB => 0;
use constant LOG_FILE  => 'compute_tfbs_cluster_counts.log';

use constant OPOSSUM_DB_HOST => 'opossum.cmmt.ubc.ca';

use constant CLUSTER_DB_NAME => 'TFBS_cluster';
use constant CLUSTER_DB_HOST => 'opossum.cmmt.ubc.ca';

my $opossum_db_name;
my $opossum_db_host;
my $cluster_db_name = CLUSTER_DB_NAME;
my $cluster_db_host = CLUSTER_DB_HOST;
my $start_gid;
my $end_gid;
my $cid;
my $out_file;
my $log_file;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
    'cd=s'  => \$cluster_db_name,
    'ch=s'  => \$cluster_db_host,
    'cid=i' => \$cid,
    's=i'   => \$start_gid,
    'e=i'   => \$end_gid,
    'o=s'   => \$out_file,
    'l=s'   => \$log_file
);

if (!$opossum_db_name) {
    pod2usage(
        -msg        => "No oPOSSUM DB name specified",
        -verbose    => 1
    );
}

if (!$opossum_db_host) {
#    pod2usage(
#        -msg        => "No oPOSSUM DB host specified",
#        -verbose    => 1
#    );
    $opossum_db_host = OPOSSUM_DB_HOST;
}

if (!$out_file) {
    pod2usage(
        -msg        => "No output TFBS cluster counts file specified",
        -verbose    => 1
    );
}

if (!$log_file) {
    $log_file = LOG_FILE;
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
my $appender = Log::Log4perl::Appender->new(
    "Log::Dispatch::File",
    filename => $log_file,
    mode     => "write"
);
#my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("%M:%L %p: %m%n");
$appender->layout($layout);
$logger->add_appender($appender);

my $start_time = time;
my $localtime  = localtime($start_time);

$logger->info("compute_tfbs_cluster_counts started $localtime\n");

my $opdba = OPOSSUM::DBSQL::DBAdaptor->new(
    -dbname   => $opossum_db_name,
    -host     => $opossum_db_host,
    -user     => 'opossum_r',
    -password => ''
);

if (!$opdba) {
    $logger->logdie("could not connect to oPOSSUM database\n" . $DBI::errstr);
}

my $cldba = TFBSCluster::DBSQL::DBAdaptor->new(
    -dbname   => $cluster_db_name,
    -host     => $cluster_db_host,
    -user     => 'opossum_r',
    -password => ''
);

if (!$cldba) {
    $logger->logdie("could not connect to TFBS_cluster database\n" . $DBI::errstr);
}

#
# Get some OPOSSUM adaptors up front
#
my $ga = $opdba->get_GeneAdaptor
    || $logger->logdie("getting GeneAdaptor");

my $srla = $opdba->get_SearchRegionLevelAdaptor()
    || $logger->logdie("getting SearchRegionLevelAdaptor");

my $cla = $opdba->get_ConservationLevelAdaptor()
    || $logger->logdie("getting ConservationLevelAdaptor");

my $tla = $opdba->get_ThresholdLevelAdaptor()
    || $logger->logdie("getting ThresholdLevelAdaptor");

my $ctfbsa = $opdba->get_ConservedTFBSAdaptor()
    || $logger->logdie("getting ConservedTFBSAdaptor");

my $cl_levels = $cla->fetch_levels()
    || $logger->logdie("getting conservation levels");

my $thl_levels = $tla->fetch_levels()
    || $logger->logdie("getting threshold levels");

my $srl_levels = $srla->fetch_levels()
    || $logger->logdie("getting search region levels");

my $cl_hash = $cla->fetch_conservation_level_hash()
    || $logger->logdie("getting conservation levels hash");

my $thl_hash = $tla->fetch_threshold_level_hash()
    || $logger->logdie("getting threshold levels hash");

my $srl_hash = $srla->fetch_search_region_level_hash()
    || $logger->logdie("getting search region levels hash");

#
# Get some TFBSCluster adaptors up front
#
my $tca = $cldba->get_TFClusterAdaptor
    || $logger->logdie("getting TFClusterAdaptor");

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

my $tfclusters;
if ($cid) {
    push @$tfclusters, $tca->fetch_by_cluster_id($cid);
} else {
    $tfclusters = $tca->fetch_all;
}

# first, go through all tfs and create a mapping btw tfid and clusterid
# I will be keeping arrays for each hash-key
my %tf_to_cluster_ids;
foreach my $cluster (@$tfclusters)
{
    my $tfids = $cluster->tf_ids;
    foreach my $tfid (@$tfids)
    {
        $tf_to_cluster_ids{$tfid} = $cluster->id;
    }
}

open(OCFH, ">$out_file") || $logger->logdie("opening output file $out_file");
foreach my $gid (@$gids)
{
    $logger->info("Processing gene $gid");

    foreach my $cl_level (@$cl_levels) {
        $logger->debug("Conservation level $cl_level");

        my $cl = $cl_hash->{$cl_level};

        foreach my $srl_level (@$srl_levels) {
            $logger->debug("Search region level $srl_level");

            my $srl = $srl_hash->{$srl_level};
            my $upstream_bp     = $srl->upstream_bp();
            my $downstream_bp   = $srl->downstream_bp();

            foreach my $thl_level (@$thl_levels) {
                $logger->debug("Threshold level $thl_level");

                my $thl = $thl_hash->{$thl_level};
                my $threshold = $thl->threshold();
                
                my $cons_tfbs_set;
                
                foreach my $cluster (@$tfclusters)
                {
                    my $cid = $cluster->id;
                    my $tf_ids = $cluster->tf_ids;
                    $cons_tfbs_set = OPOSSUM::ConservedTFBSSet->new();
                    foreach my $tfid (@$tf_ids)
                    {
                        my $cons_set = $ctfbsa->fetch_set_by_gene_tf_id(
                            -gene_id            => $gid,
                            -tf_id              => $tfid,
                            -conservation_level => $cl_level,
                            -threshold          => $threshold,
                            -upstream_bp        => $upstream_bp,
                            -downstream_bp      => $downstream_bp
                        );
                        $cons_tfbs_set->add_tf_site_set($cons_set);
                    }
                    
                    my $cons_tfbss = $cons_tfbs_set->tf_sites('start');
                    next if !$cons_tfbss or scalar @$cons_tfbss == 0;
                    
                    # merge the sites by cluster
                    my $merged_sites = merge_cluster_sites($cons_tfbss, $cid);
                    my $count = $merged_sites ? scalar(@$merged_sites) : 0;
                    
                    my $total_length = 0;
                    foreach my $site (@$merged_sites) {
                        $total_length += length($site->seq);
                    }

                    #$logger->debug(
                    #    sprintf "$gid\t%d\t$cl_level\t$thl_level"
                    #        ."\t$srl_level\t%d\t%d\n",
                    #            $cid, $count, $total_length
                    #);

                    printf OCFH "$gid\t%s\t$cl_level\t$thl_level"
                        . "\t$srl_level\t%d\t%d\n",
                            $cid, $count, $total_length;
                }
            }
        }
    }
}
close(OCFH);

my $end_time = time;
$localtime = localtime($end_time);
my $elapsed_secs = $end_time - $start_time;

$logger->info("compute_tfbs_cluster_counts completed $localtime");
$logger->info("Elapsed time (s): $elapsed_secs");

###########################

sub merge_cluster_sites
{
    my ($tfsites, $cluster_id) = @_;
    
    if (!defined $tfsites) {
        $logger->error("merge_cluster_sites: No tf sites provided");
        return;
    }
    
    my @merged_sites;
    push @merged_sites, $$tfsites[0];
    for (my $i = 1; $i < scalar(@$tfsites); $i++)
    {        
        my $tfsite = $$tfsites[$i];
        my $prevsite = $merged_sites[$#merged_sites];
        $prevsite->id($cluster_id);
        
        # if overlap, keep the max score
        # merge the two sites
        if (overlap($prevsite, $tfsite))
        {
            if ($prevsite->end < $tfsite->end) {
                
                # merge the sequences
                # first, check the strands of the sites
                # if negative, reverse complement
                # I should only do this if they are overlapping
				
                if ($prevsite->strand != $tfsite->strand) {
                    if ($prevsite->strand == -1) {
                        my $seq = Bio::Seq->new(-seq => $prevsite->seq);
                        $prevsite->seq($seq->revcom->seq);
                    } else {
                        my $seq = Bio::Seq->new(-seq => $tfsite->seq);
                        $tfsite->seq($seq->revcom->seq);
                    }
				}
				
				my $ext_seq = substr($tfsite->seq, $prevsite->end - $tfsite->start + 1);               
				$prevsite->end($tfsite->end);                
                $prevsite->seq($prevsite->seq . $ext_seq);
            }

            if ($tfsite->rel_score > $prevsite->rel_score) {
                    $prevsite->rel_score($tfsite->rel_score);
            }
            if ($tfsite->score > $prevsite->score) {
                    $prevsite->score($tfsite->score);
            }

        } else {
            $tfsite->id($cluster_id);
            push @merged_sites, $tfsite;
        }
    }
    
    return \@merged_sites;
}

sub overlap
{
    my ($tf1, $tf2) = @_;
    
    if (($tf1->start <= $tf2->start and $tf1->end > $tf2->start)
        or ($tf2->start <= $tf1->start and $tf2->end > $tf1->start))
    {
        return 1;
    }
    return 0;
}

