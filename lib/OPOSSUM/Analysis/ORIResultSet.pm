=head1 NAME

OPOSSUM::ORIResultSet.pm - module to hold the results of an ORI analysis

=head1 DESCRIPTION

Implements a set of ORIResult objects

=head1 AUTHOR

 Shannan Ho Sui
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: shosui@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::ORIResultSet;

use strict;

use Carp;

use OPOSSUM::Analysis::ORIResult;

=head2 new

 Title    : new
 Usage    : $ors = OPOSSUM::Analysis::ORIResultSet->new();
 Function : Create a new OPOSSUM::Analysis::ORIResultSet object.
 Returns  : An OPOSSUM::Analysis::ORIResultSet object.
 Args     : None.

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
    			_result_array	=> [],
			%args
		    }, ref $class || $class;

    return $self;
}

=head2 param

 Title    : param
 Usage    : $value = $ors->param($param) or $ors->param($param, $value);
 Function : Get/set a parameter value
 Returns  : A parameter value
 Args     : [1] a parameter name,
 	    [2] on set, a parameter value

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

=head2 num_results

 Title    : num_results
 Usage    : $num = $ors->num_results();
 Function : Return the number of results in the set
 Returns  : An integer
 Args     : None

=cut

sub num_results
{
    return @{$_[0]->{_result_array}} ? scalar @{$_[0]->{_result_array}} : 0;
}

=head2 add_result

 Title    : add_result
 Usage    : $ors->add_result($result);
 Function : Add a new result to the set
 Returns  : Nothing
 Args     : An OPOSSUM::Analysis::ORIResult object

=cut

sub add_result
{
    my ($self, $new_result) = @_;

    return if !defined $new_result;

    if (!$new_result->isa("OPOSSUM::Analysis::ORIResult")) {
    	carp "not an OPOSSUM::Analysis::ORIResult";
	return;
    }

    my $new_id = $new_result->id;
    foreach my $result (@{$self->{_result_array}}) {
	if ($new_id eq $result->id) {
	    carp "result with ID $new_id already exists in set";
	    return;
	}
    }

    push @{$self->{_result_array}}, $new_result;

    return $new_result;
}

=head2 get_result

 Title    : get_result
 Usage    : $result = $ors->get_result($idx);
 Function : Return a result from the set by it's index in the set
 Returns  : An OPOSSUM::Analysis::ORIResult object
 Args     : Index of the result in the set

=cut

sub get_result
{
    my ($self, $idx) = @_;

    return if !defined $idx;

    return if $idx >= $self->num_results;

    return $self->{_result_array}->[$idx];
}

=head2 get_result_by_id

 Title    : get_result_by_id
 Usage    : $result = $ors->get_result_by_id($id);
 Function : Return a result from the set by it's unique ID
 Returns  : An OPOSSUM::Analysis::ORIResult object
 Args     : ID of the result in the set

=cut

sub get_result_by_id
{
    my ($self, $id) = @_;

    return if !defined $id;

    foreach my $result (@{$self->{_result_array}}) {
    	if ($result->id eq $id) {
	    return $result;
	}
    }
    return undef;
}

=head2 sort_by

 Title    : sort_by
 Usage    : $ors->sort_by($field, $reverse);
 Function : Sort the result set by the given field
 Returns  : Nothing
 Args     : [1] an OPOSSUM::Analysis::ORIResult field to sort the set on
	    [2] a boolean which, if true, causes the sort to be reversed

=cut

sub sort_by
{
    my ($self, $field, $reverse) = @_;

    return if !defined $field;

    return if !@{$self->{_result_array}};

    if ($reverse) {
	$self->{_result_array} = [sort {$b->$field <=> $a->$field}
    					@{$self->{_result_array}}];
    } else {
	$self->{_result_array} = [sort {$a->$field <=> $b->$field}
    					@{$self->{_result_array}}];
    }
}

1;
