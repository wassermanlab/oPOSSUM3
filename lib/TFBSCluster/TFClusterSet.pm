=head1 NAME

TFClusterSet.pm - module to hold a set of TFCluster objects

=head1 DESCRIPTION

Implements a set of TFCluster objects

=head1 AUTHOR

 Andrew Kwon (adapting the script by David Arenillas)
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut
package TFBSCluster::TFClusterSet;

use strict;

use TFBSCluster::TFCluster;

use Carp;

=head2 new

 Title	: new
 Usage	: $tfcs = TFBSCluster::TFClusterSet->new();
 Function : Create a new TFBSCluster::TFClusterSet object.
 Returns  : A TFBSCluster::TFClusterSet object.
 Args	 : None.

=cut

sub new
{
	my ($class, %args) = @_;

	my $self = bless {
			_params	=> {},
			_tfc_hash	=> {} 
		}, ref $class || $class;
    
    if ($args{-tf_cluster_list}) {
        $self->add_tf_cluster_list($args{-tf_cluster_list});
    }

    if ($args{-tf_cluster_set}) {
        $self->add_tf_cluster_set($args{-tf_cluster_set});
    }

	return $self;
}

=head2 size

 Title	: size
 Usage	: $size = $tfcs->size();
 Function : Return the size of the set (number of TFCluster objects)
 Returns  : An integer
 Args	 : None

=cut

sub size
{
	return scalar keys %{$_[0]->{_tfc_hash}} || 0;
}

=head2 ids

 Title    : ids
 Usage    : $ids = $tfcs->ids();
 Function : Return the ids of the TFCluster objects in the set
 Returns  : A list/listref of TFCluster IDs
 Args     : None

=cut

sub ids
{
    my @ids = keys %{$_[0]->{_tfc_hash}};

    #if (wantarray()) {
    #    return @ids;
    #} else {
        return @ids ? \@ids : undef;
    #}
}

=head2 cluster_ids

 Title	  : cluster_ids
 Usage	  : $ids = $tfcs->cluster_ids();
 Function : Alternate name for ids()
 Returns  : A list/listref of TFCluster IDs
 Args	  : None

=cut

sub cluster_ids
{
    return $_[0]->ids();
}

=head2 get_tf_classes

 Title	: get_tf_classes
 Usage	: $tff = $tfcs->get_tf_classes();
 Function : Returns the tf class names for the TFs in this cluster set
 Returns  : A listref of string (tf class names)
 Args	 : None

=cut

sub get_tf_classes
{
    my $self = shift;
    
    my %classes;
    my @classes;
    foreach my $cid (keys %{$self->{_tfc_hash}})
    {
		my $tfc = $self->{_tfc_hash}->{$cid};
		my $class = $tfc->class;
		if (!$class or $class eq "") {
			$class = "Unknown";
		}
		if (!$classes{$class}) {
			$classes{$class} = 1;
			push @classes, $class;
		}
    }
    
	my @sorted_classes = sort {$a cmp $b} @classes;
    return @sorted_classes ? \@sorted_classes : undef;
}

=head2 get_tf_families

 Title	: get_tf_families
 Usage	: $tff = $tfcs->get_tf_families();
 Function : Returns the tf family names for the TFs in this cluster set
 Returns  : A listref of string (tf family names)
 Args	 : None

=cut

sub get_tf_families
{
    my $self = shift;
    
    my %families;
    my @families;
    foreach my $cid (keys %{$self->{_tfc_hash}})
    {
		my $tfc = $self->{_tfc_hash}->{$cid};
		my $family = $tfc->family;
		if (!$family or $family eq "") {
			$family = "Unknown";
		}
		if (!$families{$family}) {
		    $families{$family} = 1;
		    push @families, $family;
		}
    }
    
	my @sorted_families = sort {$a cmp $b} @families;
    return @sorted_families ? \@sorted_families : undef;
}

=head2 get_tf_cluster

 Title	: get_tf_cluster
 Usage	: $tfc = $tfcs->get_tf_cluster($id);
 Function : Return a TFCluster object from the set by its ID
 Returns  : A TFBSCluster::TFCluster object
 Args	 : ID of the TFCluster object

=cut

sub get_tf_cluster
{
	my ($self, $id) = @_;

	return $self->{_tfc_hash}->{$id};
}

=head2 add_tf_cluster

 Title	: add_tf_cluster
 Usage	: $tfcs->add_tf_cluster($tfc);
 Function : Add a new TFCluster object to the set
 Returns  : Nothing
 Args	 : A TFBSCluster::TFCluster object

=cut

sub add_tf_cluster
{
	my ($self, $tfc) = @_;

	return if !$tfc;
	if (!ref $tfc || !$tfc->isa("TFBSCluster::TFCluster")) {
		carp "not a TFBSCluster::TFCluster object";
		return;
	}
	
    my $cid = $tfc->id();
    if ($self->{_tfc_hash}->{$cid}) {
        carp "TFCluster $cid is already in the set (IDs must be unique)";
        return;
    }
	
	$self->{_tfc_hash}->{$tfc->id} = $tfc;
}

=head2 get_tf_cluster_by_tf_id

 Title	: get_tf_cluster_by_tf_id
 Usage	: $cluster = $tfcs->get_tf_cluster_by_tf_id($tfid);
 Function: Get the cluster containing the specified TFInfo
 Returns: A TFBSCluster::TFCluster object
 Args	: A TF ID integer

=cut

sub get_tf_cluster_by_tf_id
{
	my ($self, $tf_id) = @_;
	
	return !$tf_id;
	
	foreach my $cid (keys %{$self->{_tfc_hash}})
	{
		my $tfc = $self->{_tfc_hash}->{$cid};
		return $tfc if $tfc->contains_tf_id($tf_id);
	}
	
	return;
}

=head2 add_tf_cluster_list

 Title    : add_tf_cluster_list
 Usage    : $cs->add_tf_cluster_list();
 Function : Add a list of TFCluster objects to the set
 Returns  : Nothing
 Args     : A listref of TFCluster objects

=cut

sub add_tf_cluster_list
{
    my ($self, $tf_cluster_list) = @_;

    return unless $tf_cluster_list && $tf_cluster_list->[0];
	
    unless ($tf_cluster_list->[0]->isa('TFBSCluster::TFCluster')) {
        carp("Not a TFBSCluster::TFCluster listref");
        return;
    }

    foreach my $tf_cluster (@$tf_cluster_list) {
        $self->add_tf_cluster($tf_cluster);
    }
}

=head2 add_tf_cluster_set

 Title    : add_tf_cluster_set
 Usage    : $tfs->add_tf_cluster_set();
 Function : Add a set of TFBSCluster::TFCluster objects to the set
 Returns  : Nothing
 Args     : A TFBSCluster::TFClusterSet object

=cut

sub add_tf_cluster_set
{
    my ($self, $tf_cluster_set) = @_;

    return unless $tf_cluster_set;;

    unless ($tf_cluster_set->isa('TFBSCluster::TFClusterSet')) {
        carp("Not an TFBSCluster::TFClusterSet");
        return;
    }

    return unless $tf_cluster_set->size > 0;

    my $clusters = $tf_cluster_set->get_tf_cluster_list;
    foreach my $tf_cluster (@$clusters) {
        $self->add_tf_cluster($tf_cluster);
    }
}

=head2 get_tf_cluster_list

 Title    : get_tf_cluster_list
 Usage    : $tf_cluster_list = $tfcs->get_tf_cluster_list($sort_field);
 Function : Return a list of the TFCluster objects in the set
 Returns  : A listref of TFBSCluster::TFCluster objects
 Args     : Optionally a field to sort TFCluster objects on (e.g. id)

=cut

sub get_tf_cluster_list
{
    my ($self, $sort_field) = @_;

    my $ids = $self->ids();

    return if !$ids;

    my @tf_cluster_list;
    foreach my $id (@$ids) {
        push @tf_cluster_list, $self->get_tf_cluster($id);
    }

    if ($sort_field) {
        if (uc $sort_field eq 'ID') {
            @tf_cluster_list = sort {$a->id() cmp $b->id()} @tf_cluster_list;
        } elsif (uc $sort_field eq 'NAME') {
			@tf_cluster_list = sort {$a->name() cmp $b->name()} @tf_cluster_list;
		}
    }

    return @tf_cluster_list ? \@tf_cluster_list : undef;
}

=head2 get_tf_cluster_set

 Title    : get_tf_cluster_set
 Usage    : $tf_cluster_set = $tfcs->get_tf_cluster_set();
 Function : Return a TFClusterSet from this set (copy)
 Returns  : A TFBSCluster::TFClusterSet object
 Args     : None.

=cut

sub get_tf_cluster_set
{
    my $self = shift;

    my $ids = $self->ids();

    return if !$ids;

    my $tf_cluster_set = TFBSCluster::TFClusterSet->new();

    foreach my $id (@$ids) {
        $tf_cluster_set->add_tf_cluster($self->get_tf_cluster($id));
    }

    return $tf_cluster_set;
}

=head2 subset

 Title	: subset
 Usage	: $subset = $tfcs->subset($ids);
 Function : Return a subset of this set based on the provided IDs
 Returns  : A TFBSCluster::TFClusterSet object
 Args	 : A ref to a list of TFCluster IDs

=cut

sub subset
{
	my ($self, %args) = @_;

    my $ids        = $args{-ids};
    my $tf_classes = $args{-tf_classes};

	my $subset = TFBSCluster::TFClusterSet->new();
	
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
            my $tfc = $self->get_tf_cluster($id);
            $subset->add_tf_cluster($tfc) if $tfc;
        }
    } else {
        my %tf_class_hash;
        if ($tf_classes) {
            if (ref $tf_classes eq 'ARRAY') {
                # -tf_class arg value is a listref of values
                foreach my $class (@$tf_classes) {
                    $tf_class_hash{$class} = 1;
                }
            } else {
                # -tf_classes arg value is a single value
                $tf_class_hash{$tf_classes} = 1;
            }
        }

        # Get all the TFs
        foreach my $id (@{$self->ids()}) {
            my $tfc = $self->get_tf_cluster($id);
            if ($tfc) {
                my $include = 1;
                if ($include && %tf_class_hash) {
					foreach my $class (@{$tfc->tf_classes}) {
						$include = 0
							if !$tf_class_hash{$class};
					}
                }
				
                if ($include) {
                    $subset->add_tf_cluster($tfc);
                }
            }
        }
    }

    #
    # First set subset params to same as this set
    #
    foreach my $param_name ($self->param()) {
        $subset->param($param_name, $self->param($param_name));
    }

    #
    # Overwrite subset params depending on args passed
    #
    foreach my $arg_key (keys %args) {
        if ($arg_key eq '-tf_class')
        {
            my $arg_name = $arg_key;
            $arg_name =~ s/^-//;
            
            $subset->param($arg_name, $args{$arg_key});
        }
    }

    return $subset;
}


=head2 param

 Title    : param
 Usage    : $value = $tfs->param($param)
            or $tfs->param($param, $value);
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

1;
