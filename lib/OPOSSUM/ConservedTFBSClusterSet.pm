=head1 NAME

ConservedTFBSSet.pm - module to hold a set of ConservedTFBS objects

=head1 DESCRIPTION

Implements a set of ConservedTFBS objects

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

    if ($sort_by) {
        if ($sort_by eq 'name' || $sort_by eq 'tfbs_name') {
            return sort {$a->tfbs_name cmp $b->tfbs_name
                        || $a->start <=> $b->start
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'score' || $sort_by eq 'score') {
            return sort {$a->score <=> $b->score
                        || $a->start <=> $b->start
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'start' || $sort_by eq 'start') {
            return sort {$a->start <=> $b->start
                        || $a->tfbs_name cmp $b->tfbs_name
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'end' || $sort_by eq 'end') {
            return sort {$a->end <=> $b->end
                        || $a->tfbs_name cmp $b->tfbs_name
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } elsif ($sort_by eq 'conservation') {
            return sort {$a->conservation <=> $b->conservation
                        || $a->start <=> $b->start
                        || $a->tfbs_name cmp $b->tfbs_name
                        || $a->strand cmp $b->strand
            } @{$self->{_site_array}};
        } else {
            carp "unrecognized sort field";
        }
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

1;
