#!/usr/local/bin/perl -w

=head1 NAME

compute_conserved_region_gc_content.pl (deprecated)

=head1 SYNOPSIS

  compute_conserved_region_gc_content.pl -d opossum_db_name -h opossum_db_host
    [-s start] [-e end] -o out_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -d opossum_db_name   = Name of ORCA db we are working on.
   -h opossum_db_host   = Host name of the ORCA db.
   -s start             = Start computing with this gene id.
   -e end               = Finish computing with this gene id.
   -i conserved_tfbss_file = conserved_tfbss.txt file with missing gc_content column
   -o out_file          = Ouput conserved region gc_content file (for import
                          into oPOSSUM conserved_region_gc_content table
                          using mysqlimport).
   -l log_file          = Name of log file to which processing and error
                          messages are written.
                          (Default = compute_conserved_region_gc_content.log)

=head1 DESCRIPTION

This is the oPOSSUM_2010 script for computing the GC content of conserved regions.
The genes, promoters, and conserved_regions tables must already be populated.

=head1 ALGORITHM

=head1 NOTE

Should use compute_conserved_regions_filter_exons.pl intead.

=head1 AUTHOR

  Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca

=cut

use lib "/space/devel/oPOSSUM_2010_cluster/lib";
use lib "/raid2/tjkwon/oPOSSUM_2010/lib";
use lib '/usr/local/src/ensembl-57/ensembl/modules';
use lib '/raid2/local/src/ensembl-57/ensembl/modules';

use strict;

use Getopt::Long;
use Pod::Usage;
use OPOSSUM::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Log::Log4perl qw(get_logger :levels);

use constant DEBUG	    => 0;
use constant LOG_FILE   => 'compute_conserved_region_gc_content.log';

use constant OPOSSUM_DB_USER	=> "opossum_r";
use constant OPOSSUM_DB_PASS	=> "";

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
use constant WORM => 'worm';
use constant YEAST => 'yeast';

my $opossum_db_name;
my $opossum_db_host = "vm5.cmmt.ubc.ca";
my $species;
my $start_gid;
my $end_gid;
my $in_file;
my $out_file;
my $log_file;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
	'sp=s'	=> \$species,
    's=i'   => \$start_gid,
    'e=i'   => \$end_gid,
	'i=s'	=> \$in_file,
    'o=s'   => \$out_file,
    'l=s'   => \$log_file
);

if (!$species) {
    pod2usage(
        -msg        => "Please provide the species name for analysis",
        -verbose    => 1
    );
}
$species = lc($species);

if (!$opossum_db_name) {
	$opossum_db_name = "oPOSSUM_2010_$species";
}

#if (!$opossum_db_host) {
#    pod2usage(
#        -msg        => "Please profile an oPOSSUM DB host name",
#        -verbose    => 1
#    );
#}

if ($in_file and ($start_gid or $end_gid)) {
	pod2usage(
		-msg		=> "No gene range allowed if conserved_tfbss.txt given",
		-verbose	=> 1
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

$logger->info("compute_conserved_region_gc_content started on $localtime\n");

# select the right ENSEMBL DB name
my $ens_db_name;
if ($species eq HUMAN) {
	$ens_db_name = ENSEMBL_DB_NAME_HUMAN;
} elsif ($species eq MOUSE) {
	$ens_db_name = ENSEMBL_DB_NAME_MOUSE;
} elsif ($species eq FLY) {
	$ens_db_name = ENSEMBL_DB_NAME_FLY;
} elsif ($species eq WORM) {
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

# connect to OPOSSUM DB
my $opdbh = OPOSSUM::DBSQL::DBAdaptor->new(
    -dbname		=> $opossum_db_name,
    -host		=> $opossum_db_host,
    -user		=> OPOSSUM_DB_USER,
    -password	=> OPOSSUM_DB_PASS
);

if (!$opdbh) {
    $logger->logdie("could not connect to OPOSSUM DB\n");
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

if ($in_file)
{
	parse_in_file($in_file, $out_file);
	exit;
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

open(OFH, ">$out_file")
    || $logger->logdie("opening output conserved region gc_content file\n");

foreach my $gid (@$gids)
{
    $logger->info("Processing gene $gid");

    my $gene = $ga->fetch_by_id($gid);
    if (!$gene) {
        $logger->logdie("fetching gene $gid");
    }

	# opossum gene starts and ends have upstream and downstream offset built in
	# - for downstream, up to 3' end of the gene only
	my $sr_start = $gene->start;
	my $sr_end = $gene->end;
	my $slice = $ens_sa->fetch_by_region('chromosome', $gene->chr, $sr_start, $sr_end);
	my $search_seq = $slice->seq;
	
	foreach my $clevel (@$cl_levels)
	{
		my $crs = $gene->fetch_conserved_regions($clevel);
		if (!$crs) {
			$logger->warn("No conserved regions for gene $gid at level $clevel");
			next;
		}
		
		foreach my $cr (@$crs)
		{
	#		my $reg_start = $sr_start + $cr->start - 1;
	#		my $reg_end = $sr_start + $cr->end;
			my $reg_length = $cr->end - $cr->start + 1;
			
			# finish this part
			my $cl_level = $cr->conservation_level;
			my $cr_start = $cr->start;
			my $cr_end = $cr->end;
			my $cr_score = $cr->conservation;
			my $cr_seq = substr($search_seq, $cr_start, $reg_length);
			$logger->info("cl level: $cl_level\t$cr_start-$cr_end");
			$logger->info("$cr_seq");
			my $gc_content = calculate_gc_content($cr_seq);
			printf OFH "%d\t%d\t%d\t%d\t%.3f\t%.3f\n",
				$gid, $cl_level, $cr_start, $cr_end, $cr_score, $gc_content;
		}
	}    
}

close(OFH);

sub parse_in_file
{
	my ($in_file, $out_file) = @_;

	open(IFH, "$in_file")
		|| $logger->logdie("opening input conserved_tfbss.file\n");
	open(OFH, ">$out_file")
	    || $logger->logdie("opening output conserved region gc_content file\n");
	
	my $gene;
	my $slice;
	my $search_seq;
	my $sr_start;
	my $sr_end;
	while (my $line = <IFH>)
	{
		chomp($line);
		my (@row) = split "\t", $line;
		
		if (!defined $gene or $gene->id != $row[0])
		{
			$gene = $ga->fetch_by_id($row[0]);
			if (!$gene) {
			    $logger->logdie("fetching gene " . $row[0]);
			}
			# opossum gene starts and ends have upstream and downstream offset built in
			# - for downstream, up to 3' end of the gene only
			$sr_start = $gene->start;
			$sr_end = $gene->end;
			$slice = $ens_sa->fetch_by_region('chromosome', $gene->chr, $sr_start, $sr_end);
			$search_seq = $slice->seq;
		}
		
		my $reg_start = $sr_start + $row[2] - 1;
		my $reg_end = $sr_start + $row[3];
		my $reg_length = $row[3] - $row[2] + 1;
		
		# finish this part
		my $cl_level = $row[1];
		my $cr_seq = substr($search_seq, $reg_start, $reg_length);
		my $gc_content = calculate_gc_content($cr_seq);
		printf OFH "%s\t%.3f\n", $line, $gc_content;
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
