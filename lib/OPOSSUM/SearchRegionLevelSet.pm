=head1 NAME

OPOSSUM::SearchRegionLevelSet.pm - module to hold a set of SearchRegionLevel
objects

=head1 DESCRIPTION

Implements a set of SearchRegionLevel objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::SearchRegionLevelSet;

use strict;

use Carp;

use OPOSSUM::SearchRegionLevel;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title    : new
 Usage    : $srls = OPOSSUM::SearchRegionLevelSet->new();
 Function : Create a new OPOSSUM::SearchRegionLevelSet object.
 Returns  : An OPOSSUM::SearchRegionLevelSet object.
 Args     : None

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
    			_srl_array	=> [],
			%args
		    }, ref $class || $class;

    return $self;
}

=head2 num_levels

 Title    : num_levels
 Usage    : $num = $srls->num_levels();
 Function : Return the number of search region levels in the set
 Returns  : An integer
 Args     : None

=cut

sub num_levels
{
    return @{$_[0]->{_srl_array}} ? scalar @{$_[0]->{_srl_array}} : 0;
}

=head2 get_levels

 Title    : get_levels
 Usage    : @levels = $srls->get_levels();
 Function : Return a list of search region levels
 Returns  : A list of integer levels
 Args     : None

=cut

sub get_levels
{
    my $self = shift;

    my @levels;
    foreach my $srl_level (@{$self->{_srl_array}}) {
	push @levels, $srl_level->level;
    }

    return sort @levels;
}


=head2 add_search_region_level

 Title    : add_search_region_level
 Usage    : $srls->add_search_region_level($srl);
 Function : Add a new SearchRegionLevel object to the set
 Returns  : Nothing
 Args     : An OPOSSUM::SearchRegionLevel object

=cut

sub add_search_region_level
{
    my ($self, $new_srl) = @_;

    return if !defined $new_srl;

    if (!$new_srl->isa("OPOSSUM::SearchRegionLevel")) {
    	carp "not an OPOSSUM::SearchRegionLevel";
	return;
    }

    foreach my $srl (@{$self->{_srl_array}}) {
    	if ($new_srl->level == $srl->level) {
	    carp "search_region level already exists";
	    return;
	}
    }

    push @{$self->{_srl_array}}, $new_srl;

    return $new_srl;
}

=head2 get_search_region_level

 Title    : get_search_region_level
 Usage    : $srl = $cls->get_search_region_level($idx);
 Function : Get a SearchRegionLevel object from the set by it's index
 Returns  : An OPOSSUM::SearchRegionLevel object
 Args     : The index of the SearchRegionLevel object in the set

=cut

sub get_search_region_level
{
    my ($self, $idx) = @_;

    return if !defined $idx;

    return if $idx > $self->num_levels;

    return $self->{_srl_array}->[$idx];
}

=head2 get_search_region_level_by_level

 Title    : get_search_region_level_by_level
 Usage    : $srl = $cls->get_search_region_level_by_level($level);
 Function : Get a SearchRegionLevel object from the set by it's level
 	    number
 Returns  : An OPOSSUM::SearchRegionLevel object
 Args     : The level number of the SearchRegionLevel object in the set

=cut


sub get_search_region_level_by_level
{
    my ($self, $level) = @_;

    return if !$level;

    foreach my $srl (@{$self->{_srl_array}}) {
    	if ($srl->level == $level) {
	    return $srl;
	}
    }

    return undef;
}

=head2 sort

 Title    : sort
 Usage    : $srls->sort($reverse);
 Function : Sort the list of search region levels in the set by their
	    level numbers
 Returns  : Nothing
 Args     : [1] Optionally a boolean which, if true, reverses the sort
	    order

=cut

sub sort
{
    my ($self, $reverse) = @_;

    return if !@{$self->{_srl_array}};

    if ($reverse) {
	$self->{_srl_array} = [sort {$b->level <=> $a->level}
    					@{$self->{_srl_array}}];
    } else {
	$self->{_srl_array} = [sort {$a->level <=> $b->level}
    					@{$self->{_srl_array}}];
    }
}

1;
