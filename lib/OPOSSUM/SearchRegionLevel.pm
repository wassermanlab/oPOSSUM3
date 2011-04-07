=head1 NAME

OPOSSUM::SearchRegionLevel - SearchRegionLevel object (search_region_levels
DB record)

=head1 DESCRIPTION

A SearchRegionLevel object models a record retrieved from the
search_region_levels table of the oPOSSUM DB. This table stores the various
amounts of upstream and downstream nucleotides around the TSS that
were used to search for TFBSs.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut
package OPOSSUM::SearchRegionLevel;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $srl = OPOSSUM::SearchRegionLevel->new(
			    -level		    => 1,
			    -upstream_bp	=> 5000,
			    -downstream_bp	=> 5000);

 Function: Construct a new SearchRegionLevel object
 Returns : a new OPOSSUM::SearchRegionLevel object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
		    %args
		}, ref $class || $class;

    return $self;
}

=head2 level

 Title   : level
 Usage   : $level = $srl->level() or $srl->level($level);

 Function: Get/set the search region level.
 Returns : An integer.
 Args    : None or an integer.

=cut

sub level
{
    my ($self, $level) = @_;

    if (defined $level) {
        $self->{-level} = $level;
    }
    return $self->{-level};
}

sub search_region_level
{
    my ($self, $level) = @_;

    return $self->level($level);
}

=head2 upstream_bp

 Title   : upstream_bp
 Usage   : $upbp = $srl->upstream_bp() or $srl->upstream_bp($upbp);

 Function: Get/set the amount of upstream nucleotides in the search region
           corresponding to this level.
 Returns : An integer.
 Args    : None or an integer.

=cut

sub upstream_bp
{
    my ($self, $upstream_bp) = @_;

    if (defined $upstream_bp) {
        $self->{-upstream_bp} = $upstream_bp;
    }
    return $self->{-upstream_bp};
}

=head2 downstream_bp

 Title   : downstream_bp
 Usage   : $upbp = $srl->downstream_bp() or $srl->downstream_bp($upbp);

 Function: Get/set the amount of downstream nucleotides in the search region
           corresponding to this level.
 Returns : An integer.
 Args    : None or an integer.

=cut

sub downstream_bp
{
    my ($self, $downstream_bp) = @_;

    if (defined $downstream_bp) {
        $self->{-downstream_bp} = $downstream_bp;
    }
    return $self->{-downstream_bp};
}

1;
