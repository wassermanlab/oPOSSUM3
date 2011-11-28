=head1 NAME

OPOSSUM::ConservedRegion - ConservedRegion object (conserved_regions DB record)

=head1 DESCRIPTION

A ConservedRegion object models a record retrieved from the conserved_regions
table of the oPOSSUM DB. The ConservedRegion object contains the start and
end positions of the conserved region on the gene sequence for a given level
of conservation.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut
package OPOSSUM::ConservedRegion;

use strict;
use Carp;
use OPOSSUM::DBObject;
use Bio::SeqFeature::Generic;

use vars qw(@ISA);

@ISA = qw(Bio::SeqFeature::Generic OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $cr = OPOSSUM::ConservedRegion->new(
		       -gene_id             => 1,
		       -conservation_level  => 1,
		       -start               => 1565,
		       -end                 => 1697,
               -conservation        => 0.76,
			   -gc_content			=> 0.48
           );

 Function: Construct a new ConservedRegion object
 Returns : a new OPOSSUM::ConservedRegion object

=cut

sub new
{
    my ($class, %args) = @_;

    my $obj = Bio::SeqFeature::Generic->new(%args);

    my $self = bless $obj, ref $class || $class;

    $self->adaptor($args{-adaptor});
    $self->gene_id($args{-gene_id});
    $self->conservation_level($args{-conservation_level} || $args{-level});
    $self->score($args{-conservation} || $args{-score});
	$self->gc_content($args{-gc_content});

    return $self;
}

=head2 gene_id

 Title   : gene_id
 Usage   : $gpid = $cr->gene_id() or $cr->gene_id($gpid);

 Function: Get/set the ID of the Gene that this conserved region
           is associated with.
 Returns : An integer Gene ID.
 Args    : None or an OPOSSUM::Gene ID.

=cut

sub gene_id
{
    my ($self, $gpid) = @_;

    if (defined $gpid) {
        $self->{-gene_id} = $gpid;
    }

    return $self->{-gene_id};
}

=head2 conservation_level

 Title   : conservation_level
 Usage   : $level = $cr->conservation_level()
           or $cr->conservation_level($level)

 Function: Get/set the conservation level
 Args	 : None or an integer conservation level
 Returns : Integer conservation level

=cut

sub conservation_level
{
    my ($self, $level) = @_;

    if (defined $level) {
        $self->{-conservation_level} = $level;
    }

    return $self->{-conservation_level};
}

=head2 level

 Title   : level
 Usage   : Deprecated - use conservation_level instead

=cut

sub level
{
    my ($self, $level) = @_;

    carp "deprecated method 'level'; please use 'conservation_level' instead\n";

    return $self->conservation_level($level);
}

#=head2 start
#
# Title   : start
# Usage   : $start = $cr->start() or $cr->start($start_pos);
#
# Function: Get/set the start position of the conserved region on the
#           gene sequence.
# Returns : An integer.
# Args    : None or an integer start position.
#
#=cut
#
#sub start
#{
#    my ($self, $start) = @_;
#
#    if (defined $start) {
#    	$self->{-start} = $start;
#    }
#
#    return $self->{-start};
#}
#
#=head2 end
#
# Title   : end
# Usage   : $end = $cr->end() or $cr->end($end_pos);
#
# Function: Get/set the end position of the conserved region on the
#           gene sequence.
# Returns : An integer.
# Args    : None or an integer end position.
#
#=cut
#
#sub end
#{
#    my ($self, $end) = @_;
#
#    if (defined $end) {
#    	$self->{-end} = $end;
#    }
#
#    return $self->{-end};
#}

=head2 conservation

 Title   : conservation
 Usage   : $conservation = $cr->conservation()
           or $cr->conservation($conservation);

 Function: Get/set the conservation score (percent identity) of the
           conserved region.
 Returns : A real number.
 Args    : None or a conservation score (percent identity).

=cut

sub conservation
{
    my $self = shift;

    if (@_) {
        my $score = shift;
        $self->score($score);
    }

    return $self->score;
}

=head2 gc_content

 Title   : gc_content
 Usage   : $gc_content = $cr->gc_content()
           or $cr->gc_content($gc_content);

 Function: Get/set the gc_content score (percent identity) of the
           conserved region.
 Returns : A real number.
 Args    : None or a gc_content score (percent identity).

=cut

sub gc_content
{
    my $self = shift;

    if (@_) {
        my $val = shift;
        $self->{-gc_content} = $val;
    }

    return $self->{-gc_content};
}

#=head2 length
#
# Title   : length
# Usage   : $length = $cr->length();
#
# Function: Get the length of the conserved region on the species 1
#           promoter sequence.
# Returns : An integer.
# Args    : None.
#
#=cut
#
#sub length
#{
#    my $self = shift;
#
#    return $self->end - $self->start + 1;
#}

=head2 gene

 Title   : gene
 Usage   : $gene = $cr->gene() or $cr->gene($gp);

 Function: Get/set the Gene that this conserved region is associated
           with.
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

 Function: Get the Gene object associated with this conserved region.
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
