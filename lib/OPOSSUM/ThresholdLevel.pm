=head1 NAME

OPOSSUM::ThresholdLevel - ThresholdLevel object (threshold_levels DB record)

=head1 DESCRIPTION

A ThresholdLevel object models a record retrieved from the threshold_levels
table of the oPOSSUM DB. This table stores the various matrix match thresholds
at which a TFBS is considered to be a hit.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut
package OPOSSUM::ThresholdLevel;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $tl = OPOSSUM::ThresholdLevel->new(
                    -level		=> 1,
                    -threshold	=> 75.0);

 Function: ThresholdLevel object
 Returns : a new OPOSSUM::ThresholdLevel object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
		    %args
    }, ref $class || $class;

    return $self;
}

=head2 level

 Title   : level
 Usage   : $level = $tl->level() or $tl->level($level);

 Function: Get/set the threshold level.
 Returns : An integer.
 Args    : None or an integer.

=cut

sub level
{
    my ($self, $level) = @_;

    if (defined $level) {
        $self->{-level} = $level;
    }
    return $self->{-level};
}

sub threshold_level
{
    my ($self, $level) = @_;

    return $self->level($level);
}

=head2 threshold

 Title   : threshold
 Usage   : $thresh = $cl->threshold() or $cl->threshold($thresh);

 Function: Get/set the matrix score threshold corresponding to this level.
 Returns : A real number.
 Args    : None or a real number.

=cut

sub threshold
{
    my ($self, $threshold) = @_;

    if (defined $threshold) {
        $self->{-threshold} = $threshold;
    }
    return $self->{-threshold};
}

1;
