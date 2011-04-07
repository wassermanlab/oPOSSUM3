=head1 NAME

ConservedRegionLengthSet.pm - module to hold a set of ConservedRegionLength
objects

=head1 DESCRIPTION

Implements a set of ConservedRegionLength objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ConservedRegionLengthSet;

use strict;

use Carp;
use OPOSSUM::ConservedRegionLength;

=head2 new

 Title    : new
 Usage    : $crls = OPOSSUM::ConservedRegionLengthSet->new();
 Function : Create a new OPOSSUM::ConservedRegionLengthSet object.
 Returns  : An OPOSSUM::ConservedRegionLengthSet object.
 Args     : None

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        _params  => {},
        _crl_set => {}
        },
        ref $class || $class;

    return $self;
}

=head2 param

 Title    : param
 Usage    : $value = $crls->param($param)
            or $crls->param($param, $value);
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
 Usage    : $size = $crls->size();
 Function : Return the size of the set (number of ConservedRegionLength
            objects)
 Returns  : An integer
 Args     : None

=cut

sub size
{
    return scalar keys %{$_[0]->{_crl_set}} || 0;
}

=head2 get_ids

 Title    : get_ids
 Usage    : $ids = $crls->get_ids();
 Function : Return the list of IDs of the conserved region lengths in
            the set
 Returns  : An ref to a list of IDs
 Args     : None

=cut

sub get_ids
{
    return keys %{$_[0]->{_crl_set}};
}

=head2 total_conserved_region_length

 Title    : total_conserved_region_length
 Usage    : $length = $crls->total_conserved_region_length();
 Function : Return the total length of all the conserved region length
            objects in the set.
 Returns  : An integer
 Args     : None

=cut

sub total_conserved_region_length
{
    my ($self) = @_;

    my $crl_set = $self->{_crl_set};
    my $cr_len  = 0;
    foreach my $id (keys %$crl_set) {
        $cr_len += $crl_set->{$id}->length;
    }
    return $cr_len;
}

=head2 get_conserved_region_length

 Title    : get_conserved_region_length
 Usage    : $crl = $crls->get_conserved_region_length($id);
 Function : Return the ConservedRegionLength object by it's ID
 Returns  : A ConservedRegionLength object
 Args     : ID of the ConservedRegionLength object

=cut

sub get_conserved_region_length
{
    my ($self, $id) = @_;

    return if !$id;

    return $self->{_crl_set}->{$id};
}

=head2 add_conserved_region_length

 Title    : add_conserved_region_length
 Usage    : $crls->add_conserved_region_length($crl);
 Function : Add a new ConservedRegionLength object to the set
 Returns  : Nothing
 Args     : A ConservedRegionLength object

=cut

sub add_conserved_region_length
{
    my ($self, $crl) = @_;

    return if !$crl;
    if (!ref $crl || !$crl->isa("OPOSSUM::ConservedRegionLength")) {
        carp "not an OPOSSUM::ConservedRegionLength object";
        return;
    }

    $self->{_crl_set}->{$crl->gene_id} = $crl;
}

=head2 subset

 Title    : subset
 Usage    : $subset = $crls->subset($ids);
 Function : Return a subset of this set based on the provided IDs
 Returns  : An OPOSSUM::ConservedRegionLengthSet object
 Args     : A ref to a list of conserved region length IDs

=cut

sub subset
{
    my ($self, $ids) = @_;

    return if !$ids;

    my $subset = OPOSSUM::ConservedRegionLengthSet->new();
    foreach my $id (@$ids) {
        my $crl = $self->get_conserved_region_length($id);
        if ($crl) {
            $subset->add_conserved_region_length(
                OPOSSUM::ConservedRegionLength->new(
                    -gene_id             => $crl->gene_id,
                    -conservation_level  => $crl->conservation_level,
                    -search_region_level => $crl->search_region_level,
                    -length              => $crl->length
                )
            );
        }
    }

    my @pkeys = $self->param;
    foreach my $pkey (@pkeys) {
        $subset->param($pkey, $self->param($pkey));
    }

    return $subset;
}

1;
