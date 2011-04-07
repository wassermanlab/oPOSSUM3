=head1 NAME

OPOSSUM::Exon - Exon object (models oPOSSUM DB exons table record)

=head1 DESCRIPTION

A Exon object models a record retrieved from the exons table of the oPOSSUM DB.
It contains the positional information of the exon.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Exon;

use strict;

use Carp;
use OPOSSUM::DBObject;
use Bio::SeqFeature::Generic;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $exon = OPOSSUM::Exon->new(
               -gene_id     => 1,
               -start       => 63025275,
               -end         => 63025474
           );

 Function: Construct a new Exon object
 Returns : a new OPOSSUM::Exon object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        %args
    }, ref $class || $class;

    return $self;
}

=head2 gene_id

 Title   : gene_id
 Usage   : $geneid = $exon->gene_id() or $exon->gene_id($geneid);

 Function: Get/set the ID of the Gene object associated with this
           Exon.
 Returns : A numeric ID
 Args    : None or a numeric ID

=cut

sub gene_id
{
    my ($self, $id) = @_;

    if ($id) {
        $self->{-gene_id} = $id;
    }

    return $self->{-gene_id}
}

=head2 start

 Title   : start
 Usage   : $start = $exon->start() or $exon->start($start);

 Function: Get/set the chromosomal start position of this exon
 Returns : An integer
 Args    : None or a new exon start position 

=cut

sub start
{
    my ($self, $start) = @_;

    if ($start) {
        $self->{-start} = $start;
    }

    return $self->{-start};
}

=head2 end

 Title   : end
 Usage   : $end = $exon->end() or $exon->end($end);

 Function: Get/set the chromosomal end position of this exon
 Returns : An integer
 Args    : None or a new exon end position 

=cut

sub end
{
    my ($self, $end) = @_;

    if ($end) {
        $self->{-end} = $end;
    }

    return $self->{-end};
}

=head2 biotype

 Title   : biotype
 Usage   : $biotype = $exon->biotype() or $exon->biotype($biotype);

 Function: Get/set the biotype of this exon (coding vs. UTR)
 Returns : A string
 Args    : None or a new exon biotype string

=cut

sub biotype
{
    my ($self, $biotype) = @_;

    if ($biotype) {
        $self->{-biotype} = $biotype;
    }

    return $self->{-biotype};
}

=head2 gene

 Title   : gene
 Usage   : $gene = $ex->gene() or $ex->gene($gene);

 Function: Get/set the Gene that this exon is associated with.
 Returns : An OPOSSUM::Gene object.
 Args    : None or an OPOSSUM::Gene object.

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

 Function: Get the Gene object associated with this exon.
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
