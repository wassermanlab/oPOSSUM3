=head1 NAME

OPOSSUM::TFInfoSet.pm - module to hold a set of TFInfo objects

NOTE: This module is DEPRECATED! We will not longer store TF information
in the oPOSSUM database. Instead TF (matrix) information should be retrieved
directly from the appropriate JASPAR/PAZAR database.

=head1 DESCRIPTION

Implements a set of TFInfo objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut
package OPOSSUM::TFInfoSet;

use strict;

use Carp;
use OPOSSUM::TFInfo;

=head2 new

 Title    : new
 Usage    : $tfis = OPOSSUM::TFInfoSet->new();
 Function : Create a new OPOSSUM::TFInfoSet object.
 Returns  : An OPOSSUM::TFInfoSet object.
 Args     : None.

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        _params     => {},
        _tfi_hash    => {} 
    }, ref $class || $class;

    return $self;
}

=head2 param

 Title    : param
 Usage    : $value = $tfis->param($param)
 	    or $tfis->param($param, $value);
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
 Usage    : $size = $tfis->size();
 Function : Return the size of the set (number of TFInfo objects)
 Returns  : An integer
 Args     : None

=cut

sub size
{
    return scalar keys %{$_[0]->{_tfi_hash}} || 0;
}

=head2 tf_ids

 Title    : tf_ids
 Usage    : $ids = $tfis->tf_ids();
 Function : Return the IDs of the TFInfo objects in the set
 Returns  : A reference to a list of TF IDs
 Args     : None

=cut

sub tf_ids
{
    return keys %{$_[0]->{_tfi_hash}};
}

=head2 get_tf_info

 Title    : get_tf_info
 Usage    : $tfi = $tfis->get_tf_info($id);
 Function : Return a single TFInfo object from the set by it's ID
 Returns  : An OPOSSUM::TFInfo object
 Args     : ID of the TFInfo object

=cut

sub get_tf_info
{
    my ($self, $id) = @_;

    return $self->{_tfi_hash}->{$id};
}

=head2 add_tf_info

 Title    : add_tf_info
 Usage    : $tfis->add_tf_info($tfi);
 Function : Add a new TFInfo object to the set
 Returns  : Nothing
 Args     : An OPOSSUM::TFInfo object

=cut

sub add_tf_info
{
    my ($self, $tfi) = @_;

    return if !$tfi;
    if (!ref $tfi || !$tfi->isa("OPOSSUM::TFInfo")) {
        carp "not an OPOSSUM::TFInfo object";
        return;
    }

    $self->{_tfi_hash}->{$tfi->id} = $tfi;
}

=head2 get_tf_info_list

 Title    : get_tf_info_list
 Usage    : $tfi_list = $tfis->get_tf_info_list(
                -ids         => $ids,
                -sources     => $sources,
                -collections => $collections,
                -tax_groups  => $tax_groups,
                -min_ic      => $min_ic,
                -sort_by     => $sort_field
            );
 Function : Return a listref of TFInfo objects from the set based on
            arguments passed.
 Returns  : A listref OPOSSUM::TFInfo objects
 Args     : Optionally IDs of the TFInfo objects - overides all other
                arguments
            Optionally a source or list of sources
            Optionally a collection or list of collections
            Optionally a tax group or list of tax groups
            Optionally a minimum information content
            Optionally sort returned list by TFInfo field

=cut

sub get_tf_info_list
{
    my ($self, %args) = @_;

    my @tf_info_list;

    my $ids         = $args{-ids};
    my $sources     = $args{-sources};
    my $collections = $args{-collections};
    my $tax_groups  = $args{-tax_groups};
    my $min_ic      = $args{-min_ic};
    my $sort_by     = $args{-sort_by};

    if ($ids) {
        my @id_list;
        if (ref $ids eq 'ARRAY') {
            # -ids arg value is a listref of IDs
            @id_list = @$ids;
        } else {
            # -ids arg value is a single ID
            push @id_list, $ids;
        }

        foreach my $id (@id_list) {
            my $tfi = $self->get_tf_info($id);

            if ($tfi) {
                push @tf_info_list, OPOSSUM::TFInfo->new(
                    -id             => $tfi->id,
                    -source         => $tfi->source,
                    -collection     => $tfi->collection,
                    -external_id    => $tfi->external_id,
                    -name           => $tfi->name,
                    -class          => $tfi->class,
                    -family         => $tfi->family,
                    -tax_group      => $tfi->tax_group,
                    -width          => $tfi->width,
                    -ic             => $tfi->ic
                );
            }
        }
    } else {
        my %source_hash;
        if ($sources) {
            if (ref $sources eq 'ARRAY') {
                # -sources arg value is a listref of values
                foreach my $source (@$sources) {
                    $source_hash{$source} = 1;
                }
            } else {
                # -sources arg value is a single value
                $source_hash{$sources} = 1;
            }
        }

        my %collection_hash;
        if ($collections) {
            if (ref $collections eq 'ARRAY') {
                # -collections arg value is a listref of values
                foreach my $collection (@$collections) {
                    $collection_hash{$collection} = 1;
                }
            } else {
                # -collections arg value is a single value
                $collection_hash{$collections} = 1;
            }
        }

        my %tax_group_hash;
        if ($tax_groups) {
            if (ref $tax_groups eq 'ARRAY') {
                # -tax_groups arg value is a listref of values
                foreach my $tax_group (@$tax_groups) {
                    $tax_group_hash{$tax_group} = 1;
                }
            } else {
                # -tax_groups arg value is a single value
                $tax_group_hash{$tax_groups} = 1;
            }
        }

        $min_ic = 0 if !$min_ic;

        # Get all the TFs
        foreach my $id (@{$self->tf_ids()}) {
            my $tfi = $self->get_tf_info($id);

            if ($tfi) {
                my $include = 1;

                if (%source_hash) {
                    $include = 0 if !$source_hash{$tfi->source()};
                }

                if ($include && %collection_hash) {
                    $include = 0 if !$collection_hash{$tfi->collection()};
                }

                if ($include && %tax_group_hash) {
                    $include = 0 if !$tax_group_hash{$tfi->tax_group()};
                }

                if ($include && $tfi->ic() < $min_ic) {
                    $include = 0;
                }

                if ($include) {
                    push @tf_info_list, OPOSSUM::TFInfo->new(
                        -id             => $tfi->id,
                        -source         => $tfi->source,
                        -collection     => $tfi->collection,
                        -external_id    => $tfi->external_id,
                        -name           => $tfi->name,
                        -class          => $tfi->class,
                        -family         => $tfi->family,
                        -tax_group      => $tfi->tax_group,
                        -width          => $tfi->width,
                        -ic             => $tfi->ic
                    );
                }
            }
        }
    }

    if ($sort_by) {
        $sort_by = lc $sort_by;
        if ($sort_by eq 'id' || $sort_by eq 'tf_id' || $sort_by eq 'tfid') {
            @tf_info_list = sort {$a->id <=> $a->id} @tf_info_list;
        } elsif ($sort_by eq 'external_id') {
            @tf_info_list
                = sort {$a->external_id <=> $a->external_id} @tf_info_list;
        } elsif ($sort_by eq 'name') {
            @tf_info_list = sort {$a->name cmp $a->name} @tf_info_list;
        } else {
            carp "unknown sort field or sorting not supported for $sort_by\n";
        }
    }

    return @tf_info_list ? \@tf_info_list : undef;
}

=head2 get_tf_info_set

 Title    : get_tf_info_set
 Usage    : $tfi_set = $tfis->get_tf_info_set(
                -ids         => $ids,
                -sources     => $sources,
                -collections => $collections,
                -tax_groups  => $tax_groups,
                -min_ic      => $min_ic,
                -sort_by     => $sort_field
            );
 Function : Return a subset from this set based on arguments passed.
 Returns  : An OPOSSUM::TFInfoSet object
 Args     : Optionally IDs of the TFInfo objects - overides all other
                arguments
            Optionally a source or list of sources
            Optionally a collection or list of collections
            Optionally a tax group or list of tax groups
            Optionally a minimum information content

=cut

sub get_tf_info_set
{
    my ($self, %args) = @_;

    my $tfi_list = $self->get_tf_info_set(%args);

    return undef if !$tfi_list;

    my $tfi_set = OPOSSUM::TFInfoSet->new();
    foreach my $tfi (@$tfi_list) {
        $tfi_set->add_tf_info($tfi);
    }

    foreach my $arg (keys %args) {
        $arg =~ s/^-//;

        unless ($arg eq "sort_by" || $arg eq "ids") {
            $tfi_set->param($arg, $args{"-$arg"});
        }
    }

    return $tfi_set;
}


=head2 subset

 Title    : subset
 Usage    : $subset = $tfis->subset($ids);
 Function : Return a subset of this set based on the provided IDs
 Returns  : An OPOSSUM::TFInfoSet object
 Args     : A ref to a list of TF IDs

=cut

sub subset
{
    my ($self, $ids) = @_;

    return if !$ids;

    my $subset = OPOSSUM::TFInfoSet->new();
    foreach my $id (@$ids) {
    	my $tfi = $self->get_tf_info($id);
        if ($tfi) {
            $subset->add_tf_info(
	    		OPOSSUM::TFInfo->new(
					-id             => $tfi->id,
					-source         => $tfi->source,
					-collection     => $tfi->collection,
					-external_id    => $tfi->external_id,
					-name           => $tfi->name,
					-class          => $tfi->class,
					-family         => $tfi->family,
					-tax_group      => $tfi->tax_group,
					-width          => $tfi->width,
					-ic             => $tfi->ic
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
