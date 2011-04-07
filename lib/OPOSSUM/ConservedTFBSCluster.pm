=head1 NAME

OPOSSUM::ConservedTFBS - ConservedTFBS object (conserved_tfbss DB record)

=head1 DESCRIPTION

A ConservedTFBS object models a record retrieved from the conserved_tfbss
table of the oPOSSUM DB. The ConservedTFBS object contains the start and
end positions of the TFBS site well as the PSSM score and level of
conservation.

=head1 MODIFICATIONS

 2010/03/04
 - This class now inherits from TFBS::Site

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut
package OPOSSUM::ConservedTFBS;

use strict;
use Carp;
use OPOSSUM::DBObject;
use TFBS::Site;

use vars qw(@ISA);

@ISA = qw(TFBS::Site OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $tfbs = OPOSSUM::ConservedTFBS->new(
                -id                 => 1,
                -gene_id            => 1,
                -start              => 646823928,
                -end                => 646823937,
                -strand             => 1,
                -score              => 1.897,
                -rel_score          => 0.765,
                -seq                => 'CCAAGGATAG',
                -conservation_level => 1,
                -conservation       => 0.832
            );

 Function: Construct a new ConservedTFBS object
 Returns : a new OPOSSUM::ConservedTFBS object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        %args
    }, ref $class || $class;

    return $self;
}

=head2 param

 Title    : param
 Usage    : $value = $ctfs->param($param)
            or $ctfs->param($param, $value);
 Function : Get/set the value of a parameter
 Returns  : Value of the named parameter
 Args     : [1] name of a parameter
            [2] on set, the value of the parameter

=cut

sub param
{
    my ($self, $param, $value) = @_;

    if ($param) {
        if (defined $value) {
            $self->{_params}->{$param} = $value;
        }
        return $self->{_params}->{$param};
    }
    return keys %{$self->{_params}};
}

=head2 id

 Title   : id
 Usage   : $id = $ctfs->id() or $ctfs->id($id);

 Function: Get/set the ID of the TFBS profile (matrix) associated with
           this conserved TFBS.
 Returns : The TFBS profile ID.
 Args    : None or a new ID.

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
        $self->{-id} = $id;
    }
    return $self->{-id};
}

=head2 tf_id

 Title   : tf_id
 Usage   : $id = $ctfs->tf_id() or $ctfs->tf_id($id);

 Function: Synonymous to the 'id' method.

=cut

sub tf_id
{
    my ($self, $id) = @_;

    return $self->id($id);
}

=head2 gene_id

 Title   : gene_id
 Usage   : $gene_id = $ctfs->gene_id()
           or $ctfs->gene_id($gene_id);

 Function: Get/set the ID of the Gene object associated with this
           conserved TFBS.
 Returns : A string.
 Args    : None or a new Gene ID.

=cut

sub gene_id
{
    my ($self, $gene_id) = @_;

    if (defined $gene_id) {
        $self->{-gene_id} = $gene_id;
    }
    return $self->{-gene_id};
}

=head2 start

 Title   : start
 Usage   : $start = $ctfs->start() or $ctfs->start($start);

 Function: Get/set the start position of this conserved TF site
 Returns : An integer.
 Args    : None or a new start position.

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
 Usage   : $end = $ctfs->end() or $ctfs->end($end);

 Function: Get/set the end position of this conserved TF site
 Returns : An integer.
 Args    : None or a new end position.

=cut

sub end
{
    my ($self, $end) = @_;

    if ($end) {
        $self->{-end} = $end;
    }
    return $self->{-end};
}

=head2 strand

 Title   : strand
 Usage   : $strand = $ctfs->strand() or $ctfs->strand($strand);

 Function: Get/set the strand of this conserved TFBS
 Returns : 1 or -1.
 Args    : None or a new strand.

=cut

sub strand
{
    my ($self, $strand) = @_;

    if ($strand) {
        $self->{-strand} = $strand;
    }
    return $self->{-strand};
}

=head2 seq

 Title   : seq
 Usage   : $seq = $ctfs->seq() or $ctfs->seq($seq);

 Function: Get/set the sequence of this conserved TF site
 Returns : A string.
 Args    : None or a new sequence.

=cut

sub seq
{
    my ($self, $seq) = @_;

    if ($seq) {
        $self->{-seq} = $seq;
    }
    return $self->{-seq};
}

=head2 score

 Title   : score
 Usage   : $score = $ctfs->score() or $ctfs->score($score);

 Function: Get/set the matrix score of this conserved TFBS
 Returns : A real number.
 Args    : None or a new score.

=cut

sub score
{
    my ($self, $score) = @_;

    if ($score) {
        $self->{-score} = $score;
    }
    return $self->{-score};
}

=head2 rel_score

 Title   : rel_score
 Usage   : $score = $ctfs->rel_score() or $ctfs->rel_score($score);

 Function: Get/set the matrix relative score of this conserved TFBS
 Returns : A real number.
 Args    : None or a new relative score.

=cut

sub rel_score
{
    my ($self, $score) = @_;

    if ($score) {
        $self->{-rel_score} = $score;
    }
    return $self->{-rel_score};
}

=head2 conservation_level

 Title   : conservation_level
 Usage   : $level = $ctfs->conservation_level()
           or $ctfs->conservation_level($level);

 Function: Get/set the conservation level of the region that this conserved
           TF site falls within.
 Returns : An integer.
 Args    : None or a new conservation level.

=cut

sub conservation_level
{
    my ($self, $conservation_level) = @_;

    if (defined $conservation_level) {
        $self->{-conservation_level} = $conservation_level;
    }
    return $self->{-conservation_level};
}

=head2 conservation

 Title   : conservation
 Usage   : $cons = $ctfs->conservation() or $ctfs->conservation($cons);

 Function: Get/set the conservation score of the region that this conserved
           TF site falls within.
 Returns : A real number.
 Args    : None or a new conservation score.

=cut

sub conservation
{
    my ($self, $conservation) = @_;

    if (defined $conservation) {
        $self->{-conservation} = $conservation;
    }
    return $self->{-conservation};
}

=head2 tf_info

 Title   : tf_info
 Usage   : $tf_info = $ctfs->tf_info() or $ctfs->tf_info($tf_info);

 Function: Get/set the TFInfo object associated with this conserved TFBS
 Returns : An OPOSSUM::TFInfo object.
 Args    : None or a new OPOSSUM::TFInfo object.

=cut

sub tf_info
{
    my ($self, $info) = @_;

    if (defined $info) {
        if ($info->isa("OPOSSUM::TFInfo")) {
            $self->{-info} = $info;
        } else {
            carp "not an OPOSSUM::TFInfo";
            return undef;
        }
    }

    return $self->{-info};
}

=head2 gene

 Title   : gene
 Usage   : $gp = $ctfs->gene() or $ctfs->gene($gp);

 Function: Get/set the Gene object associated with this conserved
           TFBS.
 Returns : An OPOSSUM::Gene object.
 Args    : None or a new OPOSSUM::Gene object.

=cut

sub gene
{
    my ($self, $gp) = @_;

    if ($gp) {
        if ($gp->isa("OPOSSUM::Gene")) {
            $self->{-gene} = $gp;
        } else {
            carp "not an OPOSSUM::Gene";
            return undef;
        }
    }

    return $self->{-gene};
}

1;
