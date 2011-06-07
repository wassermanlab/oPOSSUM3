=head1 NAME

TFBSCluster::TFSet.pm - module to hold a set of TFs (TFBS::Matrix objects)

=head1 DESCRIPTION

Stores TFBS::Matrix objects. It differs from a TFBS::MatrixSet in that
it uses hash based methods to index the matrices by their matrix IDs for
more efficient access to individuel matrices, instead of having to use an
Iterator object and iterate through each matrix one by one.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package TFBSCluster::TFSet;

use strict;

use Carp;
use TFBS::Matrix;

=head2 new

 Title    : new
 Usage    : $tfs = TFBSCluster::TFSet->new();
 Function : Create a new TFBSCluster::TFSet object.
 Returns  : An TFBSCluster::TFSet object.
 Args     : None.

=cut

sub new
{
    my ($class, %args) = @_;

    if ($args{-matrix_list} && $args{-matrix_set}) {
        carp "Please provide either -matrix_list or -matrix_set arg,"
            . " but not both\n";
        return;
    }

    my $self = bless {
        _params     => {},
        _matrix_hash    => {} 
    }, ref $class || $class;

    if ($args{-matrix_list}) {
        $self->add_matrix_list($args{-matrix_list});
    }

    if ($args{-matrix_set}) {
        $self->add_matrix_set($args{-matrix_set});
    }

    return $self;
}

=head2 ids

 Title    : ids
 Usage    : $ids = $tfs->ids();
 Function : Return the IDs of the matrix objects in the set
 Returns  : A list/listref of matrix IDs
 Args     : None

=cut

sub ids
{
    my @ids = keys %{$_[0]->{_matrix_hash}};

    #if (wantarray()) {
    #    return @ids;
    #} else {
        return @ids ? \@ids : undef;
    #}
}

=head2 matrix_ids

 Title    : matrix_ids
 Usage    : $ids = $tfs->matrix_ids();
 Function : Alternate name for ids()

=cut

sub matrix_ids
{
    return $_[0]->ids();
}

=head2 tf_ids

 Title    : tf_ids
 Usage    : $ids = $tfs->tf_ids();
 Function : Alternate name for ids()

=cut

sub tf_ids
{
    return $_[0]->ids();
}

=head2 size

 Title    : size
 Usage    : $size = $tfs->size();
 Function : Return the size of the set (number of TF objects)
 Returns  : An integer
 Args     : None

=cut

sub size
{
    return scalar keys %{$_[0]->{_matrix_hash}} || 0;
}

=head2 add_matrix

 Title    : add_matrix
 Usage    : $tfs->add_matrix($matrix);
 Function : Add a new matrix to the set
 Returns  : TFBS::Matrix added
 Args     : A TFBS::Matrix object

=cut

sub add_matrix
{
    my ($self, $matrix) = @_;

    return if !$matrix;

    unless (ref $matrix && $matrix->isa("TFBS::Matrix")) {
        carp "add_matrix() argument is not a TFBS::Matrix object";
        return;
    }

    my $matrix_id = $matrix->ID();
    if ($self->{_matrix_hash}->{$matrix_id}) {
        carp "Matrix $matrix_id is already in the set (IDs must be unique)";
        return;
    }

    #
    # It's probably bad to add matrix to hash as a ref instead of copying
    # but there is no TFBS::Matrix copy function or other way to easily
    # copy TFBS::Matrix objects
    #
    $self->{_matrix_hash}->{$matrix_id} = $matrix;
}

=head2 get_matrix

 Title    : get_matrix
 Usage    : $matrix = $tfs->get_matrix($id);
 Function : Return a single matrix from the set by it's ID
 Returns  : A TFBS::Matrix object
 Args     : ID of the matrix

=cut

sub get_matrix
{
    my ($self, $matrix_id) = @_;

    return if !defined $matrix_id;

    return $self->{_matrix_hash}->{$matrix_id};
}

=head2 get_tf

 Title    : get_tf
 Usage    : Synonym for get_matrix()
 Function : Return a single TF from the set by it's ID
 Returns  : A TFBS::Matrix object
 Args     : ID of the matrix

=cut

sub get_tf
{
    my ($self, $matrix_id) = @_;

    return $self->get_matrix($matrix_id);
}

=head2 add_matrix_list

 Title    : add_matrix_list
 Usage    : $tfs->add_matrix_list();
 Function : Add a list of TFBS::Matrix objects to the set
 Returns  : Nothing
 Args     : A listref of TFBS::Matrix objects

=cut

sub add_matrix_list
{
    my ($self, $matrix_list) = @_;

    return unless $matrix_list && $matrix_list->[0];

    unless ($matrix_list->[0]->isa('TFBS::Matrix')) {
        carp("Not a TFBS::Matrix listref");
        return;
    }

    foreach my $matrix (@$matrix_list) {
        $self->add_tf($matrix);
    }
}

=head2 add_matrix_set

 Title    : add_matrix_set
 Usage    : $tfs->add_matrix_set();
 Function : Add a set of TFBS::Matrix objects to the set
 Returns  : Nothing
 Args     : A TFBS::MatrixSet object

=cut

sub add_matrix_set
{
    my ($self, $matrix_set) = @_;

    return unless $matrix_set;;

    unless ($matrix_set->isa('TFBS::MatrixSet')) {
        carp("Not a TFBS::MatrixSet");
        return;
    }

    return unless $matrix_set->size > 0;

    my $iter = $matrix_set->Iterator();
    while (my $matrix = $iter->next()) {
        $self->add_matrix($matrix);
    }
}

=head2 get_matrix_list

 Title    : get_matrix_list
 Usage    : $matrix_list = $tfs->get_matrix_list($sort_field);
 Function : Return a list of the matrix objects in the set
 Returns  : A listref TFBS::Matrix objects
 Args     : Optionally a field to sort matrices on (e.g. ID or name)

=cut

sub get_matrix_list
{
    my ($self, $sort_field) = @_;

    my $ids = $self->ids();

    return if !$ids;

    my @matrix_list;
    foreach my $id (@$ids) {
        push @matrix_list, $self->get_matrix($id);
    }

    if ($sort_field) {
        if (uc $sort_field eq 'ID') {
            @matrix_list = sort {$a->ID() cmp $b->ID()} @matrix_list;
        } elsif (uc $sort_field eq 'NAME') {
            @matrix_list = sort {$a->name() cmp $b->name()} @matrix_list;
        }
    }

    return @matrix_list ? \@matrix_list : undef;
}

=head2 get_matrix_set

 Title    : get_matrix_set
 Usage    : $matrix_set = $tfs->get_matrix_set();
 Function : Return a TFBS matrix set from this set
 Returns  : A TFBS::MatrixSet object
 Args     : None.

=cut

sub get_matrix_set
{
    my $self = shift;

    my $ids = $self->ids();

    return if !$ids;

    my $matrix_set = TFBS::MatrixSet->new();

    foreach my $id (@$ids) {
        $matrix_set->add_matrix($self->get_matrix($id));
    }

    return $matrix_set;
}

=head2 subset

 Title    : subset
 Usage    : $tf_subset = $tfs->subset(
                -ids         => $ids,
                -sources     => $sources,
                -collections => $collections,
                -tax_groups  => $tax_groups,
                -min_ic      => $min_ic,
            );
 Function : Return a subset of this set based on arguments passed
 Returns  : An TFBSCluster::TFSet object
 Args     : Optionally IDs of the matrices; overides all other
                arguments
            Optionally a collection or list of collections
            Optionally a tax group or list of tax groups
            Optionally a minimum information content

=cut

sub subset
{
    my ($self, %args) = @_;

    my $ids         = $args{-ids};
    my $collections = $args{-collections};
    my $tax_groups  = $args{-tax_groups};
    my $min_ic      = $args{-min_ic};

    my $subset = TFBSCluster::TFSet->new();

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
            my $matrix = $self->get_matrix($id);

            $subset->add_matrix($matrix) if $matrix;
        }
    } else {
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
        foreach my $id (@{$self->matrix_ids()}) {
            my $matrix = $self->get_matrix($id);

            if ($matrix) {
                my $include = 1;

                if ($include && %collection_hash) {
                    $include = 0
                        if !$collection_hash{$matrix->tag('collection')};
                }

                if ($include && %tax_group_hash) {
                    $include = 0
                        if !$tax_group_hash{$matrix->tag('tax_group')};
                }

                if ($include && $matrix->to_ICM->total_ic() < $min_ic) {
                    $include = 0;
                }

                if ($include) {
                    $subset->add_matrix($matrix);
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
        if (   $arg_key eq '-collections' || $arg_key eq '-tax_groups'
            || $arg_key eq '-min_ic'
        )
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
