#!/usr/local/bin/perl -w

=head3

WS200 contains a fair number of incorrect annotation of operon-gene parings.

=cut

use strict;

use Getopt::Long;
use Pod::Usage;

my $operon_file = "/g-sp2/data/oPOSSUM_2010/build/worm/ws200/operons.txt";
my $gene_file = "/g-sp2/data/oPOSSUM_2010/build/worm/ws200/genes.txt";
my $error_file = "$ENV{HOME}/Worm/ws200/Wormbase200_all.txt.incorrect_gene_operon_matches";
my $out_file = "$operon_file.fixed";

GetOptions(
	'g=s'	=> \$gene_file,
    'op=s'   => \$operon_file,
	'e=s'	=> \$error_file,
    'o=s'  => \$out_file,
);

open (GENE, "$gene_file") or die "Can't open $gene_file\n";
open (OPERON, "$operon_file") or die "Can't open $operon_file\n";
open (ERROR, "$error_file") or die "Can't open $error_file\n";
open (OUTPUT, ">$out_file") or die "Can't write to $out_file\n";

my %gene_list;
while (my $line = <GENE>)
{
	chomp($line);
	my ($gid, $ensid, $symbol) = split "\t", $line;
	$gene_list{$gid} = $ensid;
}

my %error_list;
while (my $line = <ERROR>)
{
	chomp($line);
	my ($ensid, $opname) = split "\t", $line;
	$error_list{$ensid} = $opname;
}

while (my $line = <OPERON>)
{
	chomp($line);
	my ($opid, $opname, $geneid) = split "\t", $line;
	my $ensid = $gene_list{$geneid};
	next if !$ensid;
	if (!$error_list{$ensid}) {
		print OUTPUT "$line\n";
	}
}

