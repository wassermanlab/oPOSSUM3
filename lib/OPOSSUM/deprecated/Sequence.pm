
=head1 NAME

OPOSSUM::Sequence - Sequence object (sequences DB record)

=head1 DESCRIPTION

A Sequence object models a record retrieved from the sequences table
of the oPOSSUM DB. It contains the raw and masked sequences and the
gene ID which the sequences are related to. The sequences no longer contain
positional information which is contained in the gene object.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Sequence;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);

=head2 new

 Title   : new
 Usage   : $seq = OPOSSUM::Sequence->new(
			    -gene_id	=> 1,
			    -seq		=> 'ACTGCTGAAAA...',
			    -masked_seq	=> '...NNNTCAGAGT');

 Function: Construct a new Sequence object
 Returns : a new OPOSSUM::Sequence object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 gene_id

 Title   : gene_id
 Usage   : $gene_id = $seq->gene_id()
           or $seq->gene_id('123');

 Function: Get/set the ID of the Gene associated with this Sequence.
 Returns : A numeric Gene ID
 Args    : None or a numeric Gene ID

=cut

sub gene_id
{
    my ($self, $gene_id) = @_;

    if ($gene_id) {
        $self->{-gene_id} = $gene_id;
    }
    return $self->{-gene_id};
}

#
# No longer store positional information in the Sequence object. This is
# stored in the Gene object to which this Sequence object is related.
#
#=head2 chr
#
# Title   : chr
# Usage   : $chr = $seq->chr() or $seq->chr('9');
#
# Function: Get/set the species 1 chromosome name of the aligned sequences.
# Returns : A string
# Args    : None or a chromosome name
#
#=cut
#
#sub chr
#{
#    my ($self, $chr) = @_;
#
#    if ($chr) {
#        $self->{-chr} = $chr;
#    }
#
#    return $self->{-chr};
#}
#
#=head2 start
#
# Title   : start
# Usage   : $start = $seq->start() or $seq->start($start);
#
# Function: Get/set the species 1 sequence start position
# Returns : An integer
# Args    : None or a new gene start site 
#
#=cut
#
#sub start
#{
#    my ($self, $start) = @_;
#
#    if ($start) {
#        $self->{-start} = $start;
#    }
#
#    return $self->{-start};
#}
#
#=head2 end
#
# Title   : end
# Usage   : $end = $seq->end() or $seq->end($end);
#
# Function: Get/set the species 1 sequence end position
# Returns : An integer
# Args    : None or a new gene end site 
#
#=cut
#
#sub end
#{
#    my ($self, $end) = @_;
#
#    if ($end) {
#        $self->{-end} = $end;
#    }
#
#    return $self->{-end};
#}

=head2 seq

 Title   : seq
 Usage   : $seq = $seq->seq() or $seq->seq($seq);

 Function: Get/set the species 1 sequence
 Returns : The species 1 sequence
 Args    : None or a new species 1 sequence

=cut

sub seq
{
    my ($self, $seq) = @_;

    if ($seq) {
        $self->{-seq} = $seq;
    }

    return $self->{-seq};
}

=head2 masked_seq

 Title   : masked_seq
 Usage   : $seq = $seq->masked_seq() or $seq->masked_seq($seq);

 Function: Get/set the masked species 1 sequence
 Returns : The species 1 masked sequence
 Args    : None or a new species 1 masked sequence

=cut

sub masked_seq
{
    my ($self, $seq) = @_;

    if ($seq) {
        $self->{-masked_seq} = $seq;
    }

    return $self->{-masked_seq};
}

=head2 gene

 Title   : gene
 Usage   : $gene = $seq->gene()
           or $seq->gene($gene);

 Function: Get/set the Gene object associated with this sequence
 Returns : A Gene object
 Args    : None or a new Gene object 

=cut

sub gene
{
    my ($self, $gene) = @_;

    if ($gene) {
        if ($gene->isa("OPOSSUM::Gene")) {
            $self->{-gene} = $gene;
        } else {
            carp "not an OPOSSUM::Gene";
            return undef;
        }
    } elsif (!$self->{-gene}) {
        $self->fetch_gene();
    }

    return $self->{-gene};
}

=head2 fetch_gene

 Title   : fetch_gene
 Usage   : $gene = $seq->fetch_gene()

 Function: Get the Gene object associated with this sequence.
 Returns : A Gene object
 Args    : None

=cut

sub fetch_gene
{
    my $self = shift;

    if (!$self->adaptor()) {
        carp "no adaptor defined trying to get Gene";
        return;
    }

    my $ga = $self->adaptor()->db()->fetch_GeneAdaptor();
    if (!$ga) {
        carp "could not get GeneAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_id($self->gene_id());

    $self->{-gene} = $gene;
}

1;
