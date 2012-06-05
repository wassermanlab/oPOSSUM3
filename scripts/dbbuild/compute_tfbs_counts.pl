#!/usr/bin/perl -w

=head1 NAME

compute_tfbs_counts.pl

=head1 SYNOPSIS

  compute_tfbs_counts.pl -d db_name -h db_host [-s start] [-e end]
                         [-t tf_id_type] -o out_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

  -d db_name        = Name of oPOSSUM DB to process
  -h db_host        = Host name of oPOSSUM DB to process
  -s start          = Starting gene pair ID
  -e end            = Ending gene pair ID
  -i tf_id          = Specific TF ID for which to compute counts
  -t tf_id_type     = The TF ID 'type'. Should be the first two characters
                      in the TF ID to identify which collection it's from,
                      e.g. 'MA', 'PB', 'MF' etc.
                      This option should probably be replaced with a
                      collection name option.
  -o out_file       = Output TFBS counts file
  -l log_file       = Name of log file to which processing and error
                      messages are written
                      (Default = compute_tfbs_counts.log)

=head1 DESCRIPTION

For each gene, for each TF, at each level of conservation, promoter search
region and matrix score it counts the number of TFBSs from the conserved_tfbss
table and write out these counts to a file for importing into the tfbs_counts
table.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

#use lib '/space/devel/oPOSSUM3/lib';
use lib '/apps/oPOSSUM3/lib';

use strict;

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use OPOSSUM::DBSQL::DBAdaptor;

use constant DEBUG     => 0;
use constant UPDATE_DB => 1;
use constant LOG_FILE  => 'compute_tfbs_counts.log';

my $opossum_db_name;
my $opossum_db_host;
my $start_gid;
my $end_gid;
my $tf_id;
my $tf_id_type;
my $out_file;
my $log_file;
GetOptions(
    'd=s' => \$opossum_db_name,
    'h=s' => \$opossum_db_host,
    's=i' => \$start_gid,
    'e=i' => \$end_gid,
    'i=s' => \$tf_id,
    't=s' => \$tf_id_type,
    'o=s' => \$out_file,
    'l=s' => \$log_file
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

if (!$out_file) {
    pod2usage(
        -msg        => "No output TFBS counts file specified",
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

$logger->info("compute_tfbs_counts started $localtime\n");

my $opdba = OPOSSUM::DBSQL::DBAdaptor->new(
    -dbname   => $opossum_db_name,
    -host     => $opossum_db_host,
    -user     => 'opossum_r',
    -password => ''
);

if (!$opdba) {
    $logger->logdie("could not connect to oPOSSUM database\n" . $DBI::errstr);
}

#
# Get some adaptors up front
#
my $dbia = $opdba->get_DBInfoAdaptor
    || $logger->logdie("getting DBInfoAdaptor");

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

my $db_info = $dbia->fetch_db_info() || $logger->logdie("fetch DB info");

my $species = $db_info->species();

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

open(OCFH, ">$out_file") || $logger->logdie("opening output file $out_file");
foreach my $gid (@$gids) {
    $logger->info("Processing gene $gid");

    foreach my $cl_level (@$cl_levels) {
        $logger->debug("Conservation level $cl_level");

        my $cl = $cl_hash->{$cl_level};

        foreach my $srl_level (@$srl_levels) {
            $logger->debug("Search region level $srl_level");

            my $srl = $srl_hash->{$srl_level};
            my $upstream_bp     = $srl->upstream_bp();
            my $downstream_bp   = $srl->downstream_bp();

            #
            # Special case for yeast. Set downstream bp to undef which results
            # in search regions being extended to 3' end of gene.
            #
            #if ($species eq 'yeast') {
            #    $downstream_bp = undef;
            #}

            foreach my $thl_level (@$thl_levels) {
                $logger->debug("Threshold level $thl_level");

                my $thl = $thl_hash->{$thl_level};
                my $threshold = $thl->threshold();


                my @tf_ids;
                if ($tf_id) {
                    push @tf_ids, $tf_id;
                }

                my $tfbs_counts = $ctfbsa->fetch_gene_tfbs_counts(
                    -gene_id            => $gid,
                    -tf_ids             => @tf_ids ? \@tf_ids : undef,
                    -conservation_level => $cl_level,
                    -threshold          => $threshold,
                    -upstream_bp        => $upstream_bp,
                    -downstream_bp      => $downstream_bp
                );

                if ($tfbs_counts) {
                    foreach my $count (@$tfbs_counts) {
                        next if defined $tf_id_type
                            && $count->tf_id() !~ /^$tf_id_type/;

                        #$logger->debug(
                        #    sprintf "$gid\t%d\t$cl_level\t$thl_level"
                        #        ."\t$srl_level\t%d\n",
                        #            $count->tf_id(),
                        #            $count->count()
                        #);

                        printf OCFH "$gid\t%s\t$cl_level\t$thl_level"
                            . "\t$srl_level\t%d\n",
                                $count->tf_id(),
                                $count->count();
                    }
                }
            }
        }
    }
}
close(OCFH);

my $end_time = time;
$localtime = localtime($end_time);
my $elapsed_secs = $end_time - $start_time;

$logger->info("compute_tfbs_counts completed $localtime");
$logger->info("Elapsed time (s): $elapsed_secs");
