#!/usr/local/bin/perl -w

use strict;

use Getopt::Long;

my $num_genes;
my $in_file;
my $out_file;
GetOptions(
    'i=s'   => \$in_file,
    'o=s'   => \$out_file,
    'n=i'   => \$num_genes
);

die "No input genes file provided\n" if !$in_file;
die "No number of random genes provided\n" if !$num_genes;

my $genes = read_genes($in_file);

my $random_genes = pick_random_genes($genes, $num_genes);
die "Error computing random genes\n" if !$random_genes;

open(OFH, ">$out_file") || die "Error opening output genes file $out_file - $!\n";

foreach my $gene (@$random_genes) {
    print OFH "$gene\n";
}
close(OFH);

exit;

sub read_genes
{
    my ($file) = @_;

    open(FH, $file) || die "Could not open input genes file $file - $!\n";

    my @genes;
    while (my $line = <FH>) {
        chomp $line;

        # XXX assumes no header line

        push @genes, $line;
    }
    close(FH);

    return @genes ? \@genes : undef;
}

sub pick_random_genes
{
    my ($genes, $ngenes) = @_;

    my $total = scalar @$genes;

    my @rand_genes;
    my %rand_idxs;
    my $picked = 0;
    while ($picked < $ngenes) {
        my $rand_idx = int(rand($total));
        if (!$rand_idxs{$rand_idx}) {
            $rand_idxs{$rand_idx} = 1;
            my $gene = $genes->[$rand_idx];
            push @rand_genes, $gene;
            $picked++;
        }
    }

    return @rand_genes ? \@rand_genes : undef;
}