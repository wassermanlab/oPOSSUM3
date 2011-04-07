=head1 NAME

OPOSSUM::ThresholdLevelSet.pm - module to hold a set of ThresholdLevel objects

=head1 DESCRIPTION

Implements a set of ThresholdLevel objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ThresholdLevelSet;

use strict;

use Carp;
use OPOSSUM::ThresholdLevel;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title    : new
 Usage    : $tls = OPOSSUM::ThresholdLevelSet->new();
 Function : Create a new OPOSSUM::ThresholdLevelSet object.
 Returns  : An OPOSSUM::ThresholdLevelSet object.
 Args     : None

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        _tl_array	=> [],
        %args
    }, ref $class || $class;

    return $self;
}

=head2 num_levels

 Title    : num_levels
 Usage    : $num = $tls->num_levels();
 Function : Return the number of threshold levels in the set
 Returns  : An integer
 Args     : None

=cut

sub num_levels
{
    return @{$_[0]->{_tl_array}} ? scalar @{$_[0]->{_tl_array}} : 0;
}

=head2 get_levels

 Title    : get_levels
 Usage    : @levels = $cls->get_levels();
 Function : Return a list of threshold levels
 Returns  : A list of integer levels
 Args     : None

=cut

sub get_levels
{
    my $self = shift;

    my @levels;
    foreach my $tl_level (@{$self->{_tl_array}}) {
        push @levels, $tl_level->level;
    }

    return sort @levels;
}


=head2 add_threshold_level

 Title    : add_threshold_level
 Usage    : $tls->add_threshold_level($tl);
 Function : Add a new ThresholdLevel object to the set
 Returns  : Nothing
 Args     : An OPOSSUM::ThresholdLevel object

=cut

sub add_threshold_level
{
    my ($self, $new_tl) = @_;

    return if !defined $new_tl;

    if (!$new_tl->isa("OPOSSUM::ThresholdLevel")) {
    	carp "not an OPOSSUM::ThresholdLevel";
        return;
    }

    foreach my $tl (@{$self->{_tl_array}}) {
    	if ($new_tl->level == $tl->level) {
            carp "threshold level already exists";
            return;
        }
    }

    push @{$self->{_tl_array}}, $new_tl;

    return $new_tl;
}

=head2 get_threshold_level

 Title    : get_threshold_level
 Usage    : $tl = $cls->get_threshold_level($idx);
 Function : Get a ThresholdLevel object from the set by it's index
 Returns  : An OPOSSUM::ThresholdLevel object
 Args     : The index of the ThresholdLevel object in the set

=cut

sub get_threshold_level
{
    my ($self, $idx) = @_;

    return if !defined $idx;

    return if $idx > $self->num_levels;

    return $self->{_tl_array}->[$idx];
}

=head2 get_threshold_level_by_level

 Title    : get_threshold_level_by_level
 Usage    : $tl = $cls->get_threshold_level_by_level($level);
 Function : Get a ThresholdLevel object from the set by it's level
            number
 Returns  : An OPOSSUM::ThresholdLevel object
 Args     : The level number of the ThresholdLevel object in the set

=cut

sub get_threshold_level_by_level
{
    my ($self, $level) = @_;

    return if !$level;

    foreach my $tl (@{$self->{_tl_array}}) {
    	if ($tl->level == $level) {
            return $tl;
        }
    }

    return undef;
}

=head2 sort

 Title    : sort
 Usage    : $tls->sort($reverse);
 Function : Sort the list of threshold levels in the set by their
            level numbers
 Returns  : Nothing
 Args     : [1] Optionally a boolean which, if true, reverses the sort
            order

=cut

sub sort
{
    my ($self, $reverse) = @_;

    return if !@{$self->{_tl_array}};

    if ($reverse) {
        $self->{_tl_array} = [sort {$b->level <=> $a->level}
                              @{$self->{_tl_array}}];
    } else {
        $self->{_tl_array} = [sort {$a->level <=> $b->level}
                              @{$self->{_tl_array}}];
    }
}

1;
