#!/usr/local/bin/perl -w

=head1 NAME

fetch_worm_gene_info.pl

=head1 SYNOPSIS

  fetch_worm_gene_info.pl -i gene_ids_file -s species -d ensembl_db_name
        -u upstream_bp -og out_genes_file [-op out_promoters_file]
        [-ox out_exons_file] [-os out_seqs_file]
        [-l log_file]

=head1 ARGUMENTS

  -i gene_ids_file       = Name of input WormMart export file
  -d ensembl_db          = Name of Ensembl DB from which to extract genes/
                           exons/promoters/sequences
  -u upstream_bp         = Number of upstream bp to include with the gene
                           coordinated and associated sequence retrieved.
                           This must be AT LEAST equal to the maximum value
                           of the upstream_bp column in the
                           search_region_levels table.
                           Default = 10000
  -og out_genes_file     = Ouput genes file
  -ox out_exons_file     = Output exons file
  -o_pr out_promoters_file = Output promoters file
  -o_op out_operons_file = Output operons file
  -os out_seqs_file      = Output sequences file
  -l log_file            = Name of log file to which processing and error
                           messages are written.
                           (Default = fetch_gene_info_<input_file>.log)

=head1 DESCRIPTION

Given a list of input WormMart export file with gene and operon information,
connect to the specified Ensembl database and extract the corresponding
gene/exon/promoter/seq information needed to populate the OPOSSUM
genes/exons/promoters/sequences tables. Output the information to the
corresponding tab delimited files suitable for loading into OPOSSUM using mysqlimport.

=head1 AUTHOR

  Andrew Kwon, adapting script by David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

# Ensembl API lib
use lib '/usr/local/src/ensembl-57/ensembl/modules';
#use lib "/space/devel/oPOSSUM_2010/lib";
use lib "/space/devel/oPOSSUM_2010_cluster/lib";

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use Log::Dispatch::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use OPOSSUM::Gene;
#use OPOSSUM::Sequence;
#use OPOSSUM::Promoter;
#use OPOSSUM::Exon;

use constant DEBUG          => 0;
use constant UPSTREAM_BP    => 1500;

# Local Ensembl database settings
use constant ENSEMBL_DB_HOST    	=> 'vm2.cmmt.ubc.ca';
use constant ENSEMBL_DB_USER    	=> 'ensembl_r';
use constant ENSEMBL_DB_PASS    	=> '';
use constant ENSEMBL_DEFAULT_DB_NAME => 'caenorhabditis_elegans_core_57_200a';

use constant WORM => 'worm';

my $in_file; # input file from WormMart with gene and operon info
my $species = WORM;
my $ens_db_name = ENSEMBL_DEFAULT_DB_NAME;
my $upstream_bp;
my $out_genes_file;
my $out_operons_file;
my $out_promoters_file;
my $out_exons_file;
my $out_seqs_file;
my $log_file;
GetOptions(
    'i=s'	 => \$in_file,
    'u=i'    => \$upstream_bp,
    'd=s'	 => \$ens_db_name,
    'og=s'   => \$out_genes_file,
    'o_op=s' => \$out_operons_file,
    'o_pr=s' => \$out_promoters_file,
    'ox=s'   => \$out_exons_file,
    'os=s'   => \$out_seqs_file,
    'l=s'	 => \$log_file
);

if (!$in_file) {
    pod2usage(
        -msg        => "Please specify an input gene stable IDs file",
        -verbose    => 1
    );
}

#if (!$ens_db_name) {
#    pod2usage(
#        -msg        => "Please specify the Ensembl DB name",
#        -verbose    => 1
#    );
#}

if (!$out_genes_file) {
    pod2usage(
        -msg        => "Please specify the output genes file name",
        -verbose    => 1
    );
}

if (!$out_operons_file) {
    pod2usage(
        -msg        => "Please specify the output operons file name",
        -verbose    => 1
    );
}
$upstream_bp = UPSTREAM_BP if !defined $upstream_bp;

if (!$log_file) {
    $log_file = $in_file;
    $log_file =~ s/.*\///;
    $log_file =~ s/\..*$//;
    $log_file = "fetch_gene_info_$log_file.log";
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
my $layout = Log::Log4perl::Layout::PatternLayout->new("%p: %m%n");

$appender->layout($layout);
$logger->add_appender($appender);

my $start_time = time;
my $localtime = localtime($start_time);

$logger->info("fetch_gene_info started at $localtime\n");

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

my $ga = $ens_db->get_GeneAdaptor || $logger->logdie("getting GeneAdaptor");
my $sa = $ens_db->get_SliceAdaptor || $logger->logdie("getting SliceAdaptor");

open(IN, "$in_file")
    || $logger->logdie("opening input gene stable IDs file $in_file - $!");

open(OGFH, ">$out_genes_file")
    || $logger->logdie("opening output genes file $out_genes_file");

open(O_OPFH, ">$out_operons_file")
    || $logger->logdie("opening output operons file $out_operons_file");

if ($out_seqs_file) {
    open(OSFH, ">$out_seqs_file")
        || $logger->logdie("opening output sequences file $out_seqs_file");
}

if ($out_promoters_file) {
    open(O_PRFH, ">$out_promoters_file")
        || $logger->logdie("opening output promoters file $out_promoters_file");
}

if ($out_exons_file) {
    open(OXFH, ">$out_exons_file")
        || $logger->logdie("opening output exons file $out_exons_file");
}

# go through gene list
# organize operon list first
<IN>; # skip header

my %operon_list;
my %gene_operon_map;
my @gene_list;

my $gid = 1;
my %incl_list;
while (my $line = <IN>)
{
    chomp $line;
    my ($wb_id, $gene_symbol, $seq_id, $op_symbol, $op_chr, $op_start, $op_end)
        = split "\t", $line;
    
    next if !$seq_id; # only keep those that have a valid SeqID
    next if $incl_list{$seq_id};
    
    push @gene_list, "$seq_id\t$wb_id";
    
    if ($op_symbol)
    {
        push @{$operon_list{$op_symbol}}, $gid;
    }
    
    $incl_list{$seq_id} = $gid;
    $gid++;
}

my $opid = 1;
foreach my $operon_symbol (keys %operon_list)
{
    my $gene_ids = $operon_list{$operon_symbol};
    #print "$opid\t$operon_symbol\t" . join(',', @$gene_ids) . "\n";

    foreach my $gid (@$gene_ids)
    {
        print O_OPFH "$opid\t$operon_symbol\t$gid\n";
        $gene_operon_map{$gid} = $opid;
    }

    $opid++;
}

# now onto the genes
for (my $i = 0; $i < scalar(@gene_list); $i++)
{
    my ($seq_id, $wb_id) = split "\t", $gene_list[$i];
    my $gene_id = $i + 1;
    
    $logger->info(sprintf("Processing gene %7d\t%s", $gene_id, $seq_id));

    my $ens_gene = $ga->fetch_by_stable_id($seq_id, 1);

    if (!$ens_gene) {
        $logger->error("fetching Ensembl $species gene $seq_id");
        next;
    }

    # Nah...including everything
    # XXX skip mitochondrial genes for now
    #if ($ens_gene->slice->seq_region_name eq 'MT') {
    #    $logger->warn("mitochondrial gene - not processing for now");
    #    next;
    #}

    my $gene_name       = $ens_gene->external_name() || '';
    my $gene_biotype    = $ens_gene->biotype() || '';
    my $gene_chr        = $ens_gene->slice->seq_region_name();
    my $gene_strand     = $ens_gene->strand();
    #my $gene_desc       = $ens_gene->description() || '';

    #
    # Store gene start/end coordinates to include the upstream sequence
    #
    my $gene_start;
    my $gene_end;
    my $gene_tss;
    if ($ens_gene->strand == 1) {
        $gene_tss   = $ens_gene->start;
        $gene_end   = $ens_gene->end;
        $gene_start = $gene_tss - $upstream_bp;
    } elsif ($ens_gene->strand == -1) {
        $gene_start = $ens_gene->start;
        $gene_tss   = $ens_gene->end;
        $gene_end   = $gene_tss + $upstream_bp;
    }

    #my $gene = OPOSSUM::Gene->new(
    #    -id             => $gene_id,
    #    -ensembl_id     => $seq_id,
    #    -symbol         => $ens_gene->external_name,
    #    -description    => $ens_gene->description,
    #    -chr            => $ens_gene->slice->seq_region_name,
    #    -start          => $gene_start,
    #    -end            => $gene_end,
    #    -strand         => $ens_gene->strand,
    #);
    #
    #write_gene(\*OGFH, $gene);
    my $operon_id = $gene_operon_map{$gene_id};
    if (!$operon_id) {
        $operon_id = "";
    }
    print OGFH "$gene_id\t$seq_id\t$gene_name\t$gene_biotype\t$gene_chr"
              . "\t$gene_start\t$gene_end\t$gene_tss\t$gene_strand";
#              . "\t$operon_id\n";


    if ($out_seqs_file) {
        my $slice = $sa->fetch_by_region(
            'chromosome', $gene_chr, $gene_start, $gene_end
        );

        #my $seq = OPOSSUM::Sequence->new(
        #    -gene_id    => $gene_id,
        #    -seq        => $slice->seq,
        #    -masked_seq => $slice->get_repeatmasked_seq->seq
        #);
        #
        #write_sequence(\*OSFH, $seq);
        printf OSFH "$gene_id\t%s\t%s\n",
            $slice->seq,
            $slice->get_repeatmasked_seq->seq;
    }

    if ($out_promoters_file) {
        my $transcripts = $ens_gene->get_all_Transcripts();

        my %tss_included;
        foreach my $transcript (@$transcripts) {
            my $tss;
            if ($gene_strand == 1) {
                $tss = $transcript->start; 
            } elsif ($gene_strand == -1) {
                $tss = $transcript->end; 
            }

            # Avoid duplicate TSSs
            next if $tss_included{$tss};
            $tss_included{$tss} = 1;

            #my $promoter = OPOSSUM::Promoter->new(
            #    -gene_id                => $gene_id,
            #    -tss                    => $tss,
            #    -ensembl_transcript_id  => $transcript->stable_id || ''
            #);
            #
            #write_promoter(\*OPFH, $promoter);
            printf O_PRFH "$gene_id\t$tss\t%s\n", $transcript->stable_id || '';
        }
    }

    if ($out_exons_file) {
        my $ens_exons = $ens_gene->get_all_Exons();

        my $exons = sort_and_combine_exons($ens_exons);

        foreach my $exon (@$exons) {
            printf OXFH "$gene_id\t%d\t%d\n", $exon->start, $exon->end;
        }
    }
}
close(IN);
close(OGFH);
close(OPFH);
close(OSFH) if $out_seqs_file;
close(OPFH) if $out_promoters_file;
close(OXFH) if $out_exons_file;

my $end_time = time;
$localtime = localtime($end_time);
my $elapsed_secs = $end_time - $start_time;

$logger->info("fetch_gene_info completed at $localtime\n");
$logger->info("Elapsed time (s): $elapsed_secs");

exit;

#
# This routine sorts exons by start position and combines any overlapping
# exons (from different transcripts) into single exons.
#
# It assumes that the exons parameter passed in is a listref of some sort of
# objects with start and end methods. 
#
sub sort_and_combine_exons
{
    my ($exons) = @_;

    return if !$exons || !$exons->[0];

    my @temp_exons = sort {$a->start <=> $b->start} @$exons;
    my $num_exons = scalar @temp_exons;
    for (my $i = 0; $i < $num_exons; $i++) {
        my $exon1 = $temp_exons[$i] if exists($temp_exons[$i]);
        if ($exon1) {
            for (my $j = $i+1; $j < $num_exons; $j++) {
                my $exon2 = $temp_exons[$j] if exists ($temp_exons[$j]);
                if ($exon2) {
                    if (combine_exon($exon1, $exon2)) {
                        if ($exon2->start < $exon1->start) {
                            $exon1->start($exon2->start);
                        }

                        if ($exon2->end > $exon1->end) {
                            $exon1->end($exon2->end);
                        }
                        delete $temp_exons[$j];
                    } else {
                        last;
                    }
                }
            }
        }
    }

    my @unique_exons;
    foreach my $exon (@temp_exons) {
        if (defined $exon) {
            push @unique_exons, $exon;
        }
    }

    return @unique_exons ? \@unique_exons : undef;
}

sub combine_exon
{
    my ($exon1, $exon2) = @_;

    my $combine = 1;
    $combine = 0 if $exon1->start > $exon2->end + 1
        || $exon1->end < $exon2->start - 1;

    return $combine;
}


#
# Define the promoter search regions based on the amount of upstream/downstream
# sequence around each promoter TSS. Truncate promoter regions at 3' end of
# gene. Combine regions into unique (non-overlapping) regions.
#
sub define_search_regions
{
    my ($promoters, $upstream_bp, $downstream_bp, $strand) = @_;

    my @search_regions;
    foreach my $p (@$promoters) {
        my $ptss      = $p->tss;
        my $pstart    = $p->start;
        my $pend      = $p->end;

        my $tss_start;
        my $tss_end;
        if ($strand == 1) {
            if (defined $upstream_bp) {
                $tss_start = $ptss - $upstream_bp;
                $tss_start = $pstart if $pstart > $tss_start;
            } else {
                $tss_start = $pstart;
            }

            if (defined $downstream_bp) {
                $tss_end = $ptss + $downstream_bp - 1;
                $tss_end = $pend if $pend < $tss_end;
            } else {
                $tss_end = $pend;
            }
        } elsif ($strand == -1) {
            if (defined $upstream_bp) {
                $tss_end = $ptss + $upstream_bp;
                $tss_end = $pend if $pend < $tss_end;
            } else {
                $tss_end = $pend;
            }

            if (defined $downstream_bp) {
                $tss_start = $ptss - $downstream_bp + 1;
                $tss_start = $pstart if $pstart > $tss_start;
            } else {
                $tss_start = $pstart;
            }
        } else {
            $logger->logdie("no strand specified");
            return;
        }

        push @search_regions, Bio::SeqFeature::Generic->new(
            -start	=> $tss_start,
            -end	=> $tss_end
        );
    }

    return combine_search_regions(\@search_regions);
}

#
# Combine any promoter search regions which may overlap 
#
sub combine_search_regions
{
    my ($regs) = @_;

    return if !$regs || !$regs->[0];

    @$regs = sort {$a->start <=> $b->start} @$regs;

    my $num_regs = scalar @$regs;
    for (my $i = 0; $i < $num_regs; $i++) {
        my $reg1 = $regs->[$i] if exists($regs->[$i]);
        if ($reg1) {
            for (my $j = $i+1; $j < $num_regs; $j++) {
                my $reg2 = $regs->[$j] if exists ($regs->[$j]);
                if ($reg2) {
                    if (feature_combine($reg1, $reg2)) {
                        if ($reg2->start < $reg1->start) {
                            $reg1->start($reg2->start);
                        }

                        if ($reg2->end > $reg1->end) {
                            $reg1->end($reg2->end);
                        }
                        delete $regs->[$j];
                    } else {
                        last;
                    }
                }
            }
        }
    }

    my @unique_regs;
    foreach my $reg (@$regs) {
        if (defined $reg) {
            push @unique_regs, $reg;
        }
    }

    return @unique_regs ? \@unique_regs : undef;
}

sub write_gene
{
    my ($fh, $gene) = @_;

    printf $fh "%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n",
        $gene->id,
        $gene->seq_id,
        $gene->symbol,
        $gene->biotype,
        $gene->chr,
        $gene->start,
        $gene->end,
        $gene->tss,
        $gene->strand,
        $gene->operon_id;
}

sub write_sequence
{
    my ($fh, $seq) = @_;

    print $fh "%d\t%s\t%s\n",
        $seq->gene_id,
        $seq->seq,
        $seq->masked_seq;
}

sub write_promoter
{
    my ($fh, $promoter) = @_;

    print $fh "%d\t%d\t%s\n",
        $promoter->gene_id,
        $promoter->tss,
        $promoter->ensembl_id;
}
