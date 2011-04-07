
=head1 NAME

OPOSSUM_CLUSTER::DBSQL::Analysis::PhastConsAdaptor - Adaptor for MySQL queries to
retrieve information from the oPOSSUM database as an
ORCA::Analysis::PhastCons object.

=head1 SYNOPSIS

$caa = $db_adaptor->get_PhastConsAdaptor();

=head1 DESCRIPTION

Read necessary information from the alignments, exons and conserved_regions
tables of the oPOSSUM DB to build a PhastCons object.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM_CLUSTER::DBSQL::Analysis::PhastConsAdaptor;

use strict;

use Carp;

use Bio::LocatableSeq;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use OPOSSUM_CLUSTER::DBSQL::BaseAdaptor;
use ORCA::Analysis::PhastCons;

use vars '@ISA';
@ISA = qw(OPOSSUM_CLUSTER::DBSQL::BaseAdaptor);

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_by_gene_id

 Title    : fetch_by_gene_id
 Usage    : $phca = $phcaa->fetch_by_gene_id(
                $gene_id, $cons_level, $phc_db, $phc_track
            );
 Function : Fetch a PhastCons analysis object by it's associated
            Gene ID. If a conservation level is given, the conserved
            regions stored at this level are retrieved and attached to
            the object. Otherwise the conserved regions must be re-computed
            by calling compute_conserved_regions.
 Returns  : An ORCA::Analysis::PhastCons object.
 Args	  : A gene ID,
            Optionally a conservation level.
            Optionally a phastCons DB name
            Optionally a phastCons track (table) name

=cut

sub fetch_by_gene_id
{
    my ($self, $gid, $cons_level, $phc_db, $phc_track) = @_;

    #
    # Fetch actual Gene object
    #
    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "ERROR getting GeneAdaptor\n";
        return;
    }

    my $gene = $ga->fetch_by_gene_id($gid);
    if (!$gene) {
        carp "Could not fetch gene $gid\n";
        return;
    }

    #
    # Fetch associated sequence
    #
    #my $seqa = $self->db->get_SequenceAdaptor;
    #if (!$seqa) {
    #    carp "ERROR getting SequenceAdaptor\n";
    #    return;
    #}

    #my $op_seq = $seqa->fetch_by_gene_id($gid);
    #if (!$op_aln) {
    #    carp "Could not fetch sequences for gene $gid\n";
    #    return;
    #}
    $gene->fetch_sequence();

    #
    # Fetch associated exons
    #
    #my $exa = $self->db->get_ExonAdaptor;
    #if (!$exa) {
    #    carp "ERROR getting ExonAdaptor\n";
    #    return;
    #}

    #my $op_exons = $exa->fetch_by_gene_id($gid);
    #if (!$op_exons) {
    #    carp "Could not fetch exons for Gene $gid\n";
    #    return;
    #}
    $gene->fetch_exons();

    my ($seq, $mseq) = _create_bio_seqs($gene);

    #my $bio_exons = _opossum_to_bio_exons($op_exons);

    my $phca = ORCA::Analysis::PhastCons->new(
        -db    => $phc_db,
        -track => $phc_track,
        -seq   => $mseq,
        -chr   => $gene->chr,
        -start => $gene->start,
        -end   => $gene->end,
        -exons => $gene->exons
    );

    #
    # If conservation level parameter is specified then also fetch conserved
    # regions (otherwise they must be recomputed).
    #
    if ($cons_level) {
        my $cra = $self->db->get_ConservedRegionAdaptor;
        if (!$cra) {
            carp "ERROR getting ConservedRegionAdaptor\n";
            return;
        }

        my $op_crs = $cra->fetch_list_by_gene_id($gid, $cons_level);

        if (!$op_crs) {
            carp "ERROR fetching conserved regions for Gene $gid at"
                . " conservation level $cons_level\n";
        }
        return;

        my $cr_feats = _opossum_conserved_regions_to_features($op_crs);

        $phca->conserved_regions($cr_feats);
    }

    return $phca;
}

sub _create_bio_seqs
{
    my ($gene) = @_;

    my $seq = Bio::LocatableSeq->new(
        -primay_id => $gene->id,
        -display_id =>
            sprintf("chr%s:%d-%d", $gene->chr, $gene->start, $gene->end),
        -alphabet => "dna",
        -seq      => $gene->sequence->seq,
        -start    => $gene->start,
        -end      => $gene->end,
        -strand   => 1
    );

    my $mseq = Bio::LocatableSeq->new(
        -primary_id => $gene->id,
        -display_id =>
            sprintf("chr%s:%d-%d", $gene->chr, $gene->start, $gene->end),
        -alphabet => "dna",
        -seq      => $gene->sequence->masked_seq,
        -start    => $gene->start,
        -end      => $gene->end,
        -strand   => 1
    );

    return ($seq, $mseq);
}

sub _opossum_to_bio_exons
{
    my ($op_exons) = @_;

    my @bio_exons;
    foreach my $op_exon (@$op_exons) {
        push @bio_exons, Bio::SeqFeature::Gene::Exon->new(
            -primary_tag => 'exon',
            -source_tag  => 'oPOSSUM',
            -start       => $op_exon->rel_start,
            -end         => $op_exon->rel_end,
            -strand      => 1
        );
    }

    return @bio_exons ? \@bio_exons : undef;
}

sub _opossum_conserved_regions_to_features
{
    my ($op_crs) = @_;

    return if !$op_crs;

    my @feats;
    foreach my $op_cr (@$op_crs) {
        push @feats, Bio::SeqFeature::Generic->new(
            -source_tag => 'oPOSSUM',
            -start      => $op_cr->start,
            -end        => $op_cr->end,
            -strand     => 1,
            -score      => $op_cr->conservation
        );
    }

    return @feats ? \@feats : undef;
}
