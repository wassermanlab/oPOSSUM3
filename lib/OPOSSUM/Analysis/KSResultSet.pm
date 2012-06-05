=head1 NAME

OPOSSUM::KSResultSet.pm - module to hold the results of a KS analysis

=head1 DESCRIPTION

Implements a set of KSResult objects

=head1 AUTHOR

 Andrew Kwon, based on module by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::KSResultSet;

use strict;

use Carp;

use OPOSSUM::Analysis::KSResult;

=head2 new

 Title    : new
 Usage    : $frs = OPOSSUM::Analysis::KSResultSet->new();
 Function : Create a new OPOSSUM::Analysis::KSResultSet object.
 Returns  : An OPOSSUM::Analysis::KSResultSet object.
 Args     : None.

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        _result_hash => {},
        %args
    }, ref $class || $class;

    return $self;
}

=head2 size

 Title    : size
 Usage    : $size = $frs->size();
 Function : Return the size of the set
 Returns  : An integer
 Args     : None

=cut

sub size
{
    return scalar keys %{$_[0]->{_result_hash}} || 0;
}

=head2 num_results

 Title    : num_results
 Usage    : Synonym for size

=cut

sub num_results
{
    return $_[0]->size();
}

=head2 ids

 Title    : ids
 Usage    : $ids = $frs->ids();
 Function : Return IDs of all the results in the set
 Returns  : A list or listref of result (TF) IDs
 Args     : None

=cut

sub ids
{
    my @ids = keys %{$_[0]->{_result_hash}};

    if (wantarray()) {
        return @ids;
    } else {
        return @ids ? \@ids : undef;
    }
}

=head2 add_result

 Title    : add_result
 Usage    : $frs->add_result($result);
 Function : Add a new result to the set
 Returns  : The OPOSSUM::Analysis::KSResult object just added
 Args     : An OPOSSUM::Analysis::KSResult object

=cut

sub add_result
{
    my ($self, $result) = @_;

    return if !$result;

    unless ($result->isa("OPOSSUM::Analysis::KSResult")) {
        carp "Argument is not an OPOSSUM::Analysis::KSResult\n";
        return;
    }

    my $id = $result->id();
    if ($self->{_result_hash}->{$id}) {
        carp "Result with ID $id already exists in set\n";
        return;
    }

    $self->{_result_hash}->{$id} = $result;
}

=head2 get_result

 Title    : get_result
 Usage    : $result = $frs->get_result($id);
 Function : Return a result from the set by it's (TF) ID
 Returns  : An OPOSSUM::Analysis::KSResult object
 Args     : ID of the result (TF ID)

=cut

sub get_result
{
    my ($self, $id) = @_;

    return if !defined $id;

    my $result = $self->{_result_hash}->{$id};
    unless ($result) {
        carp "No KS result with ID $id; make sure you provided an ID"
            . " and not an index value\n";
        return;
    }

    return $result;
}

=head2 get_result_by_id

 Title    : get_result_by_id
 Usage    : Synonym for get_result()

=cut

sub get_result_by_id
{
    my ($self, $id) = @_;

    return $self->get_result($id);
}

=head2 get_list

 Title    : get_list
 Usage    : $fr_list = $frs->get_list($sort_field, $reverse);
 Function : Get a list of KS results
 Returns  : A listref of OPOSSUM::KSResult objects
 Args     : [1] a field to sort the list on
            [2] a boolean which, if true, causes the sort to be reversed

=cut

sub get_list
{
    my ($self, $field, $reverse) = @_;

    my $ids = $self->ids();

    my @list;
    foreach my $id (@$ids) {
        push @list, $self->get_result($id);
    }

    # Default sort on ID
    $field = 'id' if !$field;

    if ($reverse) {
        @list = sort {$b->$field cmp $a->$field} @list;
    } else {
        @list = sort {$a->$field cmp $b->$field} @list;
    }

    return @list ? \@list : undef;
}

=head2 param

 Title    : param
 Usage    : $value = $frs->param($param);
            $frs->param($param, $value);
 Function : Get/set a parameter value
 Returns  : A parameter value
 Args     : [1] A parameter name,
            [2] On set, a parameter value

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

=head2 sort_by

 Title    : sort_by
 Usage    : Deprecated

=cut

sub sort_by
{
    carp "sort_by() is deprecated; to get a sorted list of KS results"
        . " use the get_list() method\n";
}

1;
