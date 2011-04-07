#!/usr/local/bin/perl -w

=head1 NAME

compute_conserved_region_lengths.pl

=head1 SYNOPSIS

  compute_conserved_region_lengths.pl -d opossum_db_name -h opossum_db_host
    [-s start] [-e end] -o out_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -d opossum_db_name   = Name of ORCA db we are working on.
   -h opossum_db_host   = Host name of the ORCA db.
   -s start             = Start computing with this gene id.
   -e end               = Finish computing with this gene id.
   -o out_file          = Ouput conserved region lengths file (for import
                          into oPOSSUM conserved_region_lengths table
                          using mysqlimport).
   -l log_file          = Name of log file to which processing and error
                          messages are written.
                          (Default = compute_conserved_region_lengths.log)

=head1 DESCRIPTION

This is the oPOSSUM_2010 script for computing conserved TFBSs. The genes,
promoters, and conserved_regions tables must already be populated. All lengths
are stored including 0-lengths if there are no conserved regions which fall
into promoter search regions at a particular conservation and search region
level.

=head1 ALGORITHM

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use lib '/space/devel/oPOSSUM_2010/lib';

use strict;

use Getopt::Long;
use Pod::Usage;
use OPOSSUM::DBSQL::DBAdaptor;
use Log::Log4perl qw(get_logger :levels);

use constant DEBUG	    => 0;
use constant LOG_FILE   => 'compute_conserved_region_lengths.log';

use constant OPOSSUM_DB_USER	=> "opossum_r";
use constant OPOSSUM_DB_PASS	=> "";

my $opossum_db_name;
my $opossum_db_host;
my $start_gid;
my $end_gid;
my $out_file;
my $log_file;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
    's=i'   => \$start_gid,
    'e=i'   => \$end_gid,
    'o=s'   => \$out_file,
    'l=s'   => \$log_file
);

if (!$opossum_db_name) {
    pod2usage(
        -msg        => "Please profile an oPOSSUM DB name",
        -verbose    => 1
    );
}

if (!$opossum_db_host) {
    pod2usage(
        -msg        => "Please profile an oPOSSUM DB host name",
        -verbose    => 1
    );
}

if (!$out_file) {
    pod2usage(
        -msg        => "Please profile an output file name",
        -verbose    => 1
    );
}

$log_file = LOG_FILE if !$log_file;

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

$logger->info("compute_conserved_region_lengths started on $localtime\n");

my $opdbh = OPOSSUM::DBSQL::DBAdaptor->new(
    -dbname		=> $opossum_db_name,
    -host		=> $opossum_db_host,
    -user		=> OPOSSUM_DB_USER,
    -password	=> OPOSSUM_DB_PASS
);

if (!$opdbh) {
    $logger->logdie("could not connect to oPOSSUM DB\n");
}

my $ga = $opdbh->get_GeneAdaptor;
$logger->logdie("getting GeneAdaptor") if !$ga;

my $cra = $opdbh->get_ConservedRegionAdaptor;
$logger->logdie("getting ConservedRegionAdaptor") if !$cra;

my $srla = $opdbh->get_SearchRegionLevelAdaptor;
$logger->logdie("getting SearchRegionLevelAdaptor") if !$srla;

my $cla = $opdbh->get_ConservationLevelAdaptor;
$logger->logdie("getting ConservationLevelAdaptor") if !$cla;

my $cl_levels = $cla->fetch_levels();
$logger->logdie("getting conservation levels") if !$cl_levels;

my $srl_levels = $srla->fetch_levels();
$logger->logdie("getting search region levels") if !$srl_levels;

my $cl_hash = $cla->fetch_conservation_level_hash();
$logger->logdie("getting conservation levels") if !$cl_hash;

my $srl_hash = $srla->fetch_search_region_level_hash();
$logger->logdie("getting search region levels") if !$srl_hash;

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

open(OFH, ">$out_file")
    || $logger->logdie("opening output conserved region lengths file\n");

foreach my $gid (@$gids) {
    $logger->info("Processing gene $gid");

    my $gene = $ga->fetch_by_id($gid);
    if (!$gene) {
        $logger->logdie("fetching gene $gid");
    }

    #if (!$gene->fetch_promoters) {
    #    $logger->warn("No promoters for gene $gid");
    #    next;
    #}

    #
    # Initial check just to make sure gene has some conserved regions before
    # going into more complicated processing below.
    #
    if (!$gene->fetch_conserved_regions(1)) {
        $logger->warn("No conserved regions for gene $gid");
        next;
    }

    foreach my $cl_level (@$cl_levels) {
        $logger->debug("Conservation level $cl_level");

        foreach my $srl_level (@$srl_levels) {
            $logger->debug("Search region level $srl_level");

            my $srl = $srl_hash->{$srl_level};

            my $upstream_bp   = $srl->upstream_bp();
            my $downstream_bp = $srl->downstream_bp();

            my $length = $cra->fetch_length_by_upstream_downstream(
                $gid, $cl_level, $upstream_bp, $downstream_bp
            );

            printf OFH "$gid\t$cl_level\t$srl_level\t$length\n";
	    }
    }
}
close(OFH);
