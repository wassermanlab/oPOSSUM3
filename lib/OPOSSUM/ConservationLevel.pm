=head1 NAME

OPOSSUM::ConservationLevel - ConservationLevel object (conservation_levels
DB record)

=head1 DESCRIPTION

A ConservationLevel object models a record retrieved from the
conservation_levels table of the oPOSSUM DB. This table stores the various
levels at which conserved regions were computed, stored and searched for
TFBSs within the region around the TSSs of the promoter pairs. Conservation
was computed at different levels.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ConservationLevel;

use strict;
use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $cl = OPOSSUM::ConservationLevel->new(
			    -level              => 1,
			    -min_conservation   => 0.7
           );

 Function: Construct a new ConservationLevel object
 Returns : a new OPOSSUM::ConservationLevel object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
		    %args
    }, ref $class || $class;

    #
    # Explicitly call min_conservaton setter to check format (% string or
    # real number).
    #
    $self->min_conservation($args{-min_conservation});

    return $self;
}

=head2 level

 Title   : level
 Usage   : $level = $cl->level() or $cl->level(1);

 Function: Get/set the conservation level.
 Returns : An integer.
 Args    : None or an integer.

=cut

sub level
{
    my ($self, $level) = @_;

    if ($level) {
        $self->{-level} = $level;
    }

    return $self->{-level}
}

sub conservation_level
{
    my ($self, $level) = @_;

    return $self->level($level);
}

=head2 min_conservation

 Title   : min_conservation
 Usage   : $mp = $cl->min_conservation() or $cl->min_conservation(1);

 Function: Get/set the minimum conservation (percentage identity) of
           the conserved regions associated with this conservation level.
 Returns : A real number.
 Args    : None or a real number.

=cut

sub min_conservation
{
    my ($self, $cons) = @_;

    if (defined $cons) {
        if ($cons =~ /(\d+)%/) {
            $cons = $1 / 100;
        }
        
        if ($cons < 0 || $cons > 1) {
            carp "Minimum conservation $cons out of range. Please specify"
                . " real number in the range 0 to 1 or a percentage string"
                . " e.g. '70%'\n";
        }
        $self->{-min_conservation} = $cons;
    }

    return $self->{-min_conservation}
}

1;
