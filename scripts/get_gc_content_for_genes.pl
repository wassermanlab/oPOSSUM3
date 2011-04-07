#!/usr/local/bin/perl -w

=head1 NAME

get_gc_content_for_genes.pl

=head1 SYNOPSIS

  get_gc_content_for_genes.pl -d opossum_db_name -h opossum_db_host
    -i input_genes.txt -t gene_id_type -o out_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -d opossum_db_name   = Name of ORCA db we are working on.
   -h opossum_db_host   = Host name of the ORCA db.
   -n num genes
   -i in_file			= Input gene list
   -t gene_id_type		= Input gene id type
   -o out_file          = Ouput conserved region gc content file
   -l log_file          = Name of log file to which processing and error
                          messages are written.
                          (Default = get_gc_content_for_genes.log)

=head1 DESCRIPTION

Given the list of input genes, for each search region level, compute the 
overall GC content of the conserved regions.
Report back with GC content for each search region level, and the overall
GC content.

=head1 ALGORITHM

=head1 AUTHOR

  Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca

=cut

use strict;

use lib '/space/devel/oPOSSUM3/scripts/standalone';

use oPossumGeneInclude;

use Getopt::Long;
use Pod::Usage;
use OPOSSUM::DBSQL::DBAdaptor;
use Log::Log4perl qw(get_logger :levels);

use constant DEBUG	    => 0;
use constant LOG_FILE   => 'get_gc_content_for_genes.log';

my $opossum_db_name;
my $opossum_db_host;
my $in_file;
my $gene_id_type;
my $num_genes;
my $out_file;
my $log_file;
my $has_operon;
my $cl_level;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
	'i=s'	=> \$in_file,
	'n=i'	=> \$num_genes,
	't=s'	=> \$gene_id_type,
    'o=s'   => \$out_file,
	'op'	=> \$has_operon,
    'l=s'   => \$log_file,
	'cl=i'	=> \$cl_level
);

if (!$opossum_db_name) {
    pod2usage(
        -msg        => "Please profile an oPOSSUM DB name",
        -verbose    => 1
    );
}

if ((!$in_file and !$num_genes) or ($in_file and $num_genes)) {
    pod2usage(
        -msg        => "Please profile one of input file name or number of genes",
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
$gene_id_type = DFLT_GENE_ID_TYPE if !$gene_id_type;
$cl_level = DFLT_CONSERVATION_LEVEL;

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
    -host		=> OPOSSUM_DB_HOST,
    -user		=> OPOSSUM_DB_USER,
    -password	=> OPOSSUM_DB_PASS
);

if (!$opdbh) {
    $logger->logdie("could not connect to oPOSSUM DB\n");
}

my $ga = $opdbh->get_GeneAdaptor;
$logger->logdie("getting GeneAdaptor") if !$ga;

my $oa = $opdbh->get_OperonAdaptor;
$logger->logdie("getting OperonAdaptor") if !$oa;

my $cra = $opdbh->get_ConservedRegionAdaptor;
$logger->logdie("getting ConservedRegionAdaptor") if !$cra;

my $srla = $opdbh->get_SearchRegionLevelAdaptor;
$logger->logdie("getting SearchRegionLevelAdaptor") if !$srla;

my $crla = $opdbh->get_ConservedRegionLengthAdaptor;
$logger->logdie("getting ConservedRegionLengthAdaptor") if !$crla;

my $srl_levels = $srla->fetch_levels();
$logger->logdie("getting search region levels") if !$srl_levels;

my $srl_hash = $srla->fetch_search_region_level_hash();
$logger->logdie("getting search region levels") if !$srl_hash;

my $gids;
my $included_gene_ids;
my $missing_gene_ids;
my $gid_gene_ids;
my $operon_first_gids;
my $operon_unique_gids;

if ($in_file)
{
	my $gene_ids = read_ids_file($in_file, $logger);
	unless ($gene_ids) {
		$logger->logdie("No target gene IDs read from file $in_file");
	}

	($gids, $included_gene_ids, $missing_gene_ids, 
		$gid_gene_ids, $operon_first_gids, $operon_unique_gids
	) = fetch_opossum_gene_ids(
		$ga, $oa, $gene_id_type, $gene_ids, $has_operon);
} else {
	($gids, $operon_first_gids, $operon_unique_gids) = fetch_random_opossum_gene_ids($ga, $oa, $num_genes, $has_operon);
}

if (!$gids) {
    $logger->logdie("fetching gene IDs");
}


open(OFH, ">$out_file")
    || $logger->logdie("opening output gc content file\n");

foreach my $srl (@$srl_levels) {
    $logger->debug("Search region level $srl");

	my $sr_up = $srl_hash->{$srl}->upstream_bp();
	my $sr_down = $srl_hash->{$srl}->downstream_bp();
		
	my $gc_content = fetch_cr_gc_content(
		$ga, $cra, $gids, $cl_level, $sr_up, $sr_down);
		
	my $cr_length = $crla->fetch_total_length(
		-conservation_level		=> $cl_level,
		-search_region_level	=> $srl,
		-gene_ids				=> $gids,
		-operon_gene_ids		=> $operon_first_gids,
		-has_operon				=> $has_operon);
	
    printf OFH "$srl\t$gc_content\t$cr_length\n";
}
close(OFH);
