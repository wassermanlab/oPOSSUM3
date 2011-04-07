
=head1 NAME

OPOSSUM::Promoter - Promoter object (promoters DB record)

=head1 DESCRIPTION

A Promoter object models a record retrieved from the promoters table of the
oPOSSUM DB. It contains the gene ID of the Gene that this promoter is
associated with, the TSS and the Ensembl transcript stable ID from which this
TSS was taken. It no longer contains chr, start/end.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Promoter;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $promoter = OPOSSUM::Promoter->new(
		    -id				        => 1,
		    -gene_id		        => 1,
		    -tss			        => 103070274,
		    -ensembl_transcript_id	=> 'ENST00000369658');

 Function: Construct a new Promoter object
 Returns : a new OPOSSUM::Promoter object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

#
# No longer use a unique ID for the Promoter. We only need to know it's
# association with a Gene via the gene ID.
#
#=head2 id
#
# Title   : id
# Usage   : $id = $promoter->id() or $promoter->id('33992');
#
# Function: Get/set the ID of the Promoter. This should be a unique
#           identifier for this object within the implementation. If
#           the Promoter object was read from the oPOSSUM database,
#           this should be set to the value in the promoter_id column.
# Returns : A numeric ID
# Args    : None or a numeric ID
#
#=cut
#
#sub id
#{
#    my ($self, $id) = @_;
#
#    if ($id) {
#        $self->{-id} = $id;
#    }
#
#    return $self->{-id};
#}
#
#sub promoter_id
#{
#    my ($self, $id) = @_;
#
#    return $self->id($id);
#}

=head2 gene_id

 Title   : gene_id
 Usage   : $gpid = $promoter->gene_id() or $promoter->gene_id($gpid);

 Function: Get/set the ID of the Gene object associated with this
           Promoter.
 Returns : A numeric ID
 Args    : None or a numeric ID

=cut

sub gene_id
{
    my ($self, $id) = @_;

    if ($id) {
        $self->{-gene_id} = $id;
    }

    return $self->{-gene_id};
}

=head2 chr

 Title   : chr
 Usage   : $chr = $promoter->chr() or $promoter->chr($chr);

 Function: Get/set the chromosome name of the gene.
 Returns : A string
 Args    : None or a chromosome name

=cut

#
# No longer store chr, start, end. We just need the TSS. Start and end of
# promoter regions are computed values using the TSS and the amount
# of upstream/downstream sequence specified at the application level.
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
#=head2 strand
#
# Title   : strand
# Usage   : $strand = $promoter->strand() or $promoter->strand(1);
#
# Function: Get/set the strand of the gene from which this
#           Promoter comes.
# Returns : 1 or -1
# Args    : None or a new strand value
#
#=cut
#
#sub strand
#{
#    my ($self, $strand) = @_;
#
#    if ($strand) {
#        $self->{-strand} = $strand;
#    }
#
#    return $self->{-strand};
#}
#
#=head2 start
#
# Title   : start
# Usage   : $start = $promoter->start() or $promoter->start($start);
#
# Function: Get/set the chromosomal start position of this promoter
#           region
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
# Usage   : $end = $promoter->end() or $promoter->end($end);
#
# Function: Get/set the chromosomal end position of this promoter
# 	   region
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

=head2 tss

 Title   : tss
 Usage   : $tss = $promoter->tss() or $promoter->tss($tss);

 Function: Get/set the chromosomal position of the transcription
           start site
 Returns : An integer
 Args    : None or a new transcription start site 

=cut

sub tss
{
    my ($self, $tss) = @_;

    if ($tss) {
        $self->{-tss} = $tss;
    }

    return $self->{-tss};
}

=head2 ensembl_transcript_id

 Title   : ensembl_transcript_id
 Usage   : $id = $promoter->ensembl_transcript_id()
           or $promoter->ensembl_transcript_id($id);

 Function: Get/set the Ensembl transcript ID which the TSS of this promoter
           is based on
 Returns : An Ensembl transcript ID
 Args    : None or an Ensembl transcript ID

=cut

sub ensembl_transcript_id
{
    my ($self, $id) = @_;

    if ($id) {
        $self->{-ensembl_transcript_id} = $id;
    }

    return $self->{-ensembl_transcript_id};
}

=head2 gene

 Title   : gene
 Usage   : $gene = $promoter->gene() or $promoter->gene($gene);

 Function: Get/set the gene associated with this promoter
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

 Function: Fetch the Gene object associated with this promoter from the
           database.
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
