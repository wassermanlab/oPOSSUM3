=head1 NAME

OPOSSUM::ConservationLevelSet.pm - module to hold a set of ConservationLevel
objects

=head1 DESCRIPTION

Implements a set of ConservationLevel objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ConservationLevelSet;

use strict;

use Carp;
use OPOSSUM::ConservationLevel;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title    : new
 Usage    : $cls = OPOSSUM::ConservationLevelSet->new();
 Function : Create a new OPOSSUM::ConservationLevelSet object.
 Returns  : A OPOSSUM::ConservationLevelSet object.
 Args     : None

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
    			_cl_array	=> [],
			%args
		    }, ref $class || $class;

    return $self;
}

=head2 num_levels

 Title    : num_levels
 Usage    : $num = $cls->num_levels();
 Function : Return the number of conservation levels in the set
 Returns  : An integer
 Args     : None

=cut

sub num_levels
{
    return @{$_[0]->{_cl_array}} ? scalar @{$_[0]->{_cl_array}} : 0;
}

=head2 get_levels

 Title    : get_levels
 Usage    : @levels = $cls->get_levels();
 Function : Return a list of conservation levels
 Returns  : A list of integer levels
 Args     : None

=cut

sub get_levels
{
    my $self = shift;

    my @levels;
    foreach my $cl_level (@{$self->{_cl_array}}) {
    	push @levels, $cl_level->level;
    }

    return sort @levels;
}

=head2 add_conservation_level

 Title    : add_conservation_level
 Usage    : $cls->add_conservation_level($cl);
 Function : Add a new ConservationLevel object to the set
 Returns  : Nothing
 Args     : A OPOSSUM::ConservationLevel object

=cut

sub add_conservation_level
{
    my ($self, $new_cl) = @_;

    return if !defined $new_cl;

    if (!$new_cl->isa("OPOSSUM::ConservationLevel")) {
    	carp "not a OPOSSUM::ConservationLevel";
	return;
    }

    foreach my $cl (@{$self->{_cl_array}}) {
    	if ($new_cl->level == $cl->level) {
	    carp "conservation level already exists";
	    return;
	}
    }

    push @{$self->{_cl_array}}, $new_cl;

    return $new_cl;
}

=head2 get_conservation_level

 Title    : get_conservation_level
 Usage    : $cl = $cls->get_conservation_level($idx);
 Function : Get a ConservationLevel object from the set by it's index
 Returns  : A OPOSSUM::ConservationLevel object
 Args     : The index of the ConservationLevel object in the set

=cut

sub get_conservation_level
{
    my ($self, $idx) = @_;

    return if !defined $idx;

    return if $idx > $self->num_levels;

    return $self->{_cl_array}->[$idx];
}

=head2 get_conservation_level_by_level

 Title    : get_conservation_level_by_level
 Usage    : $cl = $cls->get_conservation_level_by_level($level);
 Function : Get a ConservationLevel object from the set by it's level
 	    number
 Returns  : A OPOSSUM::ConservationLevel object
 Args     : The level number of the ConservationLevel object in the set

=cut

sub get_conservation_level_by_level
{
    my ($self, $level) = @_;

    return if !$level;

    foreach my $cl (@{$self->{_cl_array}}) {
    	if ($cl->level == $level) {
	    return $cl;
	}
    }

    return undef;
}

=head2 sort

 Title    : sort
 Usage    : $cls->sort($reverse);
 Function : Sort the list of conservation levels in the set by their
	    level numbers
 Returns  : Nothing
 Args     : [1] Optionally a boolean which, if true, reverses the sort
	    order

=cut

sub sort
{
    my ($self, $reverse) = @_;

    return if !@{$self->{_cl_array}};

    if ($reverse) {
	$self->{_cl_array} = [sort {$b->level <=> $a->level}
    					@{$self->{_cl_array}}];
    } else {
	$self->{_cl_array} = [sort {$a->level <=> $b->level}
    					@{$self->{_cl_array}}];
    }
}

1;
