#!/usr/local/bin/perl -w

=head1 NAME

compute_conserved_regions.pl

=head1 SYNOPSIS

  compute_conserved_regions.pl -d opossum_db_name -h opossum_db_host
    -ud ucsc_db_name -ut ucsc_track_name [-s start] [-e end]
    -o cons_regs_file [-l log_file]

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -d opossum_db_name   = Name of oPOSSUM db we are working on
   -h opossum_db_host   = Host name of the ORCA db
   -s start             = Start computing with this gene id
   -e end               = Finish computing with this gene id
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

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

#
# Use most current (development) libs.
# Comment out to use installed libs.
#
use lib '/space/devel/oPOSSUM_2010/lib';
use lib '/space/devel/ORCAtk/lib';

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use Bio::LocatableSeq;
use OPOSSUM::DBSQL::DBAdaptor;
use ORCA::Analysis::PhastCons;

use constant DEBUG              => 1;

use constant LOG_FILE 		    => 'compute_conserved_tfbss.log';

use constant OPOSSUM_DB_USER	=> "opossum_r";

use constant MIN_CONSERVED_REGION_LENGTH    => 20;

my $log_file = LOG_FILE;
my $opossum_db_name;
my $opossum_db_host;
my $ucsc_db_name;
my $ucsc_track_name;
my $start_gid;
my $end_gid;
my $out_regions_file;
my $out_tfbs_file;
GetOptions(
    'd=s'   => \$opossum_db_name,
    'h=s'   => \$opossum_db_host,
    'ud=s'  => \$ucsc_db_name,
    'ut=s'  => \$ucsc_track_name,
    's=i'	=> \$start_gid,
    'e=i'	=> \$end_gid,
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

my $ga = $opdba->get_GeneAdaptor;
if (!$ga) {
    $logger->error("getting GeneAdaptor");
}

my $cla = $opdba->get_ConservationLevelAdaptor;
if (!$cla) {
    $logger->error("getting ConservationLevelAdaptor");
}

my $cl_levels = $cla->fetch_levels;

my $cl_hash = $cla->fetch_conservation_level_hash;
if (!$cl_hash) {
    $logger->error("fetching conservation levels");
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
    $logger->error("getting SearchRegionLevelAdaptor");
}

my $srl_levels = $srla->fetch_levels;

my $srl_hash = $srla->fetch_search_region_level_hash;
if (!$srl_hash) {
    $logger->error("fetching search region levels");
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

my $max_upstream_bp = $srl_hash->{1}->upstream_bp;
my $max_downstream_bp = $srl_hash->{1}->downstream_bp;

open(OCRFH, ">$out_regions_file")
    || $logger->logdie("opening output conserved regions file");

foreach my $gid (@$gids) {
    $logger->info("processing gene $gid");

    my $opgene = $ga->fetch_by_gene_id($gid);
    if (!$opgene) {
    	$logger->logdie("fetching gene");
    }

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

    #my $exons = $exa->fetch_by_gene_id($gid);

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

    my $phca = ORCA::Analysis::PhastCons->new(
        -db             => $ucsc_db_name,
        -track          => $ucsc_track_name,
        #-seq            => $biolseq,
        -chr            => $opgene->chr,
        -start          => $opgene->start,
        -end            => $opgene->end,
        -exons          => $opgene->exons
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

        write_conserved_regions(\*OCRFH, $opgene->id, $clevel, $cons_regs);
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
    my ($fh, $gid, $level, $crs) = @_;

    foreach my $cr (@$crs) {
	    printf $fh "%d\t%d\t%d\t%d\t%.3f\n",
            $gid, $level, $cr->start, $cr->end, $cr->score;
    }
}
