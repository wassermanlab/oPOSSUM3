=head1 NAME

ConservedRegionSet.pm - module to hold a set of ConservedRegion objects

=head1 DESCRIPTION

Implements a set of ConservedRegion objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=head1 MODIFICATIONS

  AK 2010/09/03
  - added method total_gc_content
  - modified to reflect the addition of gc_content to schema
  
=cut

package OPOSSUM::ConservedRegionSet;

use strict;

use Carp;
use Bio::SeqFeature::Generic;
use OPOSSUM::ConservedRegion;


=head2 new

 Title    : new
 Usage    : $crs = OPOSSUM::ConservedRegionSet->new();
 Function : Create a new OPOSSUM::ConservedRegionSet object.
 Returns  : An OPOSSUM::ConservedRegionSet object.
 Args     : None

=cut

sub new
{
    my ($class) = @_;

    my $self = bless {
        _params		=> {},
        _region_array	=> []
    }, ref $class || $class;

    return $self;
}

=head2 param

 Title    : param
 Usage    : $value = $crs->param($param)
            or $crs->param($param, $value);
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

=head2 size

 Title    : size
 Usage    : $size = $crs->size();
 Function : Return the size of the set (number of ConservedRegion objects)
 Returns  : An integer
 Args     : None

=cut

sub size
{
    my $self = shift;

    return $self->{_region_array} ? scalar @{$self->{_region_array}} : 0;
}

=head2 total_length

 Title    : total_length
 Usage    : $length = $crs->total_length();
 Function : Return the total length of the conserved regions in the set
 Returns  : An integer
 Args     : None

=cut

sub total_length
{
    my $self = shift;

    my $length = 0;
    foreach my $cr (@{$self->{_region_array}}) {
        $length += $cr->length;
    }

    return $length;
}

=head2 total_gc_content

 Title    : total_gc_content
 Usage    : $gc_content = $crs->total_gc_content();
 Function : Return the total gc_content of the conserved regions in the set
 Returns  : A float
 Args     : None

=cut

sub total_gc_content
{
    my $self = shift;

    my $sum_gc_length = 0;
	my $sum_cr_length = 0;
    foreach my $cr (@{$self->{_region_array}}) {
        my $cr_length = $cr->end - $cr->start + 1;
		my $gc_length = $cr->gc_content * $cr_length;
		$sum_gc_length += $gc_length;
		$sum_cr_length += $cr_length;
    }

	my $gc_content = $sum_gc_length / $sum_cr_length;
	
    return $gc_content;
}

=head2 conserved_regions

 Title    : conserved_regions
 Usage    : $regions = $crs->conserved_regions($sort_by);
 Function : Return the (optionally sorted) list of conserved regions
 	    in the set
 Returns  : A listref of ConservedRegion objects
 Args     : [1] Optionally a ConservedRegion field nane to sort by

=cut

sub conserved_regions
{
    my ($self, $sort_by) = @_;

    my @crs = @{$self->{_region_array}};
    if ($sort_by) {
        if ($sort_by eq 'start') {
            @crs = sort {   $a->start <=> $b->start
                         || $a->end <=> $b->end
            } @crs;
        } elsif ($sort_by eq 'end') {
            @crs = sort {   $a->end <=> $b->end
                         || $a->start <=> $b->start
            } @crs;
        } elsif ($sort_by eq 'conservation') {
            @crs = sort {   $a->conservation <=> $b->conservation
                         || $a->start <=> $b->start
            } @crs;
        } elsif ($sort_by eq 'gc_content') {
			@crs = sort {	$a->gc_content <=> $b->gc_content
						 || $a->start <=> $b->start
			} @crs;
		} else {
            carp "unrecognized sort field";
        }
    }

    return @crs ? \@crs : undef;
}

=head2 add_conserved_region

 Title    : add_conserved_region
 Usage    : $crs->add_conserved_region($cr);
 Function : Add a new ConservedRegion object to the set
 Returns  : Nothing
 Args     : An OPOSSUM::ConservedRegion object

=cut

sub add_conserved_region
{
    my ($self, $cr) = @_;

    return if !$cr;
    if (!ref $cr || !$cr->isa("OPOSSUM::ConservedRegion")) {
        carp "conserved region is not an OPOSSUM::ConservedRegion object";
        return;
    }

    my $level = $self->param('conservation_level');
    if ($level && $cr->level != $level) {
    	carp "region conservation level does not match set";
    	return;
    }
	
    push @{$self->{_region_array}}, $cr;
}

=head2 as_features

 Title    : as_features
 Usage    : $fps = $crs->as_features($sort_by);
 Function : Return the conserved region set as 
            Bio::SeqFeature::Generic
 Returns  : A reference to a list of Bio::SeqFeature::Generic objects
 Args     : [1] Optionally a ConservedRegion field nane to sort by

=cut

sub as_features
{
    my ($self, $sort_by) = @_;

    my @features;
    foreach my $cr ($self->conserved_regions($sort_by)) {
        push @features, Bio::SeqFeature::Generic->new(
			-source_tag => "oPOSSUM",
			-start		=> $cr->start,
			-end		=> $cr->end,
			-score      => $cr->conservation,
			-tag		=> { gc_content => $cr->gc_content },
        );
    }

    return @features ? \@features : undef;
}

1;
