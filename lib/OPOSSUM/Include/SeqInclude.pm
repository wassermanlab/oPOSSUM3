=head1 NAME

 OPOSSUM::Include::SeqInclude.pm

=head1 SYNOPSIS

=head1 DESCRIPTION

  Contains all options and routines that are common to all the seq-based
  analysis scripts.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

use OPOSSUM::Opt::SeqOpt;
use OPOSSUM::Include::BaseInclude;


sub read_seqs
{
    my ($file, %job_args) = @_;

    my $seqIO = Bio::SeqIO->new(-file => $file, -format => 'fasta');

    my @seqs;
    while (my $seq = $seqIO->next_seq()) {
        push @seqs, $seq;
    }
    $seqIO->close();

    return @seqs ? \@seqs : undef;
}

sub check_seqs
{
    my ($seqs, $t_or_bg, %job_args) = @_;

    foreach my $seq (@$seqs) {
        unless ($seq->seq) {
            my $msg = sprintf(
                  "Poorly formatted $t_or_bg sequences. Sequence ID %s has no"
                . " actual sequence associated with it.", $seq->display_id
            );

            fatal($msg, %job_args);
        }
    }
}

sub id_seqs
{
    my ($seqs) = @_;
    
    my @seq_ids; # seq0 .. seqN
    my %seq_id_display_ids; # mapping btw seq_id and original display_id
    my %seq_id_seqs;
    my $seq_num = 0;
    foreach my $seq (@$seqs) {
        my $seq_id = "seq$seq_num";
    
        push @seq_ids, $seq_id;
    
        my $display_id = $seq->display_id();
        if ($display_id) {
            if ($display_id =~ /^(chr\w+):(\d+)-(\d+)/) {
                $display_id = "$1:$2-$3";
            }
        } else {
            $display_id = $seq_id;
        }
    
        $seq_id_display_ids{$seq_id} = $display_id;
        $seq_id_seqs{$seq_id} = $seq;
    
        $seq_num++;
    }
    
    return (\@seq_ids, \%seq_id_seqs, \%seq_id_display_ids);
}

sub read_peak_pos
{
    my ($file, $seq_ids, %job_args) = @_;
    
    my $text = read_file($file);
    my $ids = parse_id_text($text);
    
    if (scalar @$ids != scalar @$seq_ids) {
        fatal("Number of max peak positions do not match the number of sequences",
              %job_args);
    }
    
    my %seq_peaks;
    for (my $i = 0; $i < scalar @$ids; $i++) {
        $seq_peaks{$$seq_ids[$i]} = $$ids[$i];
    }
    
    return %seq_peaks ? \%seq_peaks : undef;
}


sub calc_seqs_total_gc_content
{
    my ($seqs) = @_;

    my %count = (
        'A' => 0,
        'C' => 0,
        'G' => 0,
        'T' => 0,
        'N' => 0
    );

    my $total_length = 0;
    foreach my $seq (@$seqs) {
        my @nts = split //, $seq->seq();

        foreach my $nt (@nts) {
            $count{uc $nt}++;
        }
    }

    my $gc_count    = $count{'G'} + $count{'C'};
    my $total_count = $gc_count + $count{'A'} + $count{'T'};
    my $gc_content  = $gc_count / $total_count;

    return ($gc_content, $total_count);
}



1;
