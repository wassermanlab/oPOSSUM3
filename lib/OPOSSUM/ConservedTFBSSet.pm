=head1 NAME

ConservedTFBSSet.pm - module to hold a set of ConservedTFBS objects

=head1 DESCRIPTION

Implements a set of ConservedTFBS objects

=head1 MODIFICATION

Removed level function
Added param function
- Andrew on Aug. 11, 2010
Fixed add_tf_site_set
- Andrew on Aug. 27, 2010

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ConservedTFBSSet;

use strict;
use Carp;
use OPOSSUM::ConservedTFBS;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title    : new
 Usage    : $ctfss = OPOSSUM::ConservedTFBSSet->new();
 Function : Create a new OPOSSUM::ConservedTFBSSet object.
 Returns  : An OPOSSUM::ConservedTFBSSet object.
 Args     : level - optionally specify the conservation level of this
 		    set

=cut

sub new
{
    my ($class, %args) = @_;

    my $level = $args{-level};

    my $self = bless {
        -level	    => $level,
        _site_array	=> []
    }, ref $class || $class;

    return $self;
}

=head2 level

 Title    : level
 Usage    : $level = $ctfss->level() or $ctfss->level($level);
 Function : Get/set the conservation level of the set
 Returns  : An integer conservation level
 Args     : [1] optionally specify a new conservation level for this
            set, can only be set for a non-empty set

=cut

sub level
{
    my ($self, $level) = @_;

    if ($level) {
        if (!$self->{_site_array} || !$self->{_site_array}->[0]) {
            $self->{-level} = $level;
        } else {
            carp "can't change conservation level of non-empty set";
        }
    }
    return $self->{-level};
}

=head2 size

 Title    : size
 Usage    : $size = $ctfss->size();
 Function : Get the size of the set (number of ConservedTFBS objects
	    contained)
 Returns  : An integer size
 Args     : None

=cut

sub size
{
    my $self = shift;

    return scalar @{$self->{_site_array}} || 0;
}

=head2 tf_sites

 Title    : tf_sites
 Usage    : $tf_sites = $ctfss->tf_sites();
 Function : Return the list of conserved TF sites in the set
 Returns  : A reference to an array of ConservedTFBS objects
 Args     : None

=cut

sub tf_sites
{
    my ($self, $sort_by) = @_;

	my @sorted_array;
    if ($sort_by) {
        $sort_by = lc $sort_by;

        if ($sort_by eq 'name' || $sort_by eq 'id' || $sort_by eq 'tf_id') {
            @sorted_array = sort {$a->tf_id cmp $b->tf_id
                        || $a->start <=> $b->start
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'score' || $sort_by eq 'score') {
            @sorted_array = sort {$a->score <=> $b->score
                        || $a->start <=> $b->start
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'start' || $sort_by eq 'start') {
            @sorted_array = sort {$a->start <=> $b->start
                        || $a->tf_id cmp $b->tf_id
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'end' || $sort_by eq 'end') {
            @sorted_array = sort {$a->end <=> $b->end
                        || $a->tf_id cmp $b->tf_id
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'conservation') {
            @sorted_array = sort {$a->conservation <=> $b->conservation
                        || $a->start <=> $b->start
                        || $a->tf_id cmp $b->tf_id
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } else {
            carp "unrecognized sort field";
			return $self->{_site_array};
        }
		return \@sorted_array;
    }
	
    return $self->{_site_array};
}

=head2 add_tf_site

 Title    : add_tf_site
 Usage    : $ctfss->add_tf_site($tf_site);
 Function : Add a new conserved TF site to the set
 Returns  : Nothing
 Args     : A ConservedTFBS object

=cut

sub add_tf_site
{
    my ($self, $tf) = @_;

    return if !$tf;
    if (!ref $tf || !$tf->isa("OPOSSUM::ConservedTFBS")) {
	carp "conserved TF site is not an OPOSSUM::ConservedTFBS object";
	return;
    }
    my $level = $self->level;
    if ($level && $tf->level > $level) {
	carp "region conservation level does not match set";
	return;
    }

    push @{$self->{_site_array}}, $tf;
}

=head2 add_tf_site_list

 Title    : add_tf_site_list
 Usage    : $tfs->add_tf_site_list();
 Function : Add a list of TFBS::Matrix objects to the set
 Returns  : Nothing
 Args     : A listref of TFBS::Matrix objects

=cut

sub add_tf_site_list
{
    my ($self, $tf_site_list) = @_;

    return unless $tf_site_list && $tf_site_list->[0];

    unless ($tf_site_list->[0]->isa('OPOSSUM::ConservedTFBS')) {
        carp("Not a OPOSSUM::ConservedTFBS listref");
        return;
    }

    foreach my $site (@$tf_site_list) {
        $self->add_tf_site($site);
    }
}

=head2 add_tf_site_set

 Title    : add_tf_site_set
 Usage    : $tfs->add_tf_site_set();
 Function : Add a set of OPOSSUM::ConservedTFBS objects to the set
 Returns  : Nothing
 Args     : A OPOSSUM::ConservedTFBSSet object

=cut

sub add_tf_site_set
{
    my ($self, $tf_site_set) = @_;

    return unless $tf_site_set;;

    unless ($tf_site_set->isa('OPOSSUM::ConservedTFBSSet')) {
        carp("Not a OPOSSUM::ConservedTFBSSet");
        return;
    }

    return unless $tf_site_set->size > 0;

    my $tf_sites = $tf_site_set->tf_sites();
    foreach my $tf_site (@$tf_sites) {
        $self->add_tf_site($tf_site);
    }
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
