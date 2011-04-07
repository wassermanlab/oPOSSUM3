=head1 NAME

OPOSSUM::TFBSCount - TFBSCount object (tfbs_counts DB record)

=head1 DESCRIPTION

A TFBSCount object models a record retrieved from the tfbs_counts table of
the oPOSSUM DB. It contains the number of times a given TFBS was detected
for a given gene at a given level of conservation, PWM score threshold
and search region.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

 Modified by Shannan Ho Sui on Dec 21, 2006 to account for db schema changes

=head1 METHODS

=cut

package OPOSSUM::TFBSCount;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $tc = OPOSSUM::TFBSCount->new(
				-gene_id                => $gid,
				-tf_id                  => $tfid,
				-conservation_level     => $clevel,
				-threshold_level        => $tlevel,
				-search_region_level    => $srlevel,
				-count                  => $count);

 Function: Construct a new TFBSCount object
 Returns : a new OPOSSUM::TFBSCount object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        %args
    }, ref $class || $class;

    return $self;
}

=head2 gene_id

 Title   : gene_id
 Usage   : $id = $tc->gene_id() or $tc->gene_id($gid);
 Function: Get/set the gene ID
 Returns : A gene ID
 Args    : Optionally a gene ID

=cut

sub gene_id
{
    my ($self, $gid) = @_;

    if (defined $gid) {
        $self->{-gene_id} = $gid;
    }
    return $self->{-gene_id};
}

=head2 tf_id

 Title   : tf_id
 Usage   : $id = $tc->tf_id() or $tc->tf_id($tfid);
 Function: Get/set the TF ID
 Returns : A TF ID
 Args    : Optionally a TF ID

=cut

sub tf_id
{
    my ($self, $tf_id) = @_;

    if (defined $tf_id) {
        $self->{-tf_id} = $tf_id;
    }
    return $self->{-tf_id};
}

=head2 conservation_level

 Title   : conservation_level
 Usage   : $clevel = $tc->conservation_level()
           or $tc->conservation_level($level);
 Function: Get/set the conservation level
 Returns : An integer conservation level
 Args    : Optionally an integer conservation level

=cut

sub conservation_level
{
    my ($self, $conservation_level) = @_;

    if (defined $conservation_level) {
        $self->{-conservation_level} = $conservation_level;
    }
    return $self->{-conservation_level};
}

=head2 threshold_level

 Title   : threshold_level
 Usage   : $tlevel = $tc->threshold_level()
           or $tc->threshold_level($level);
 Function: Get/set the threshold level
 Returns : An integer threshold level
 Args    : Optionally an integer threshold level

=cut

sub threshold_level
{
    my ($self, $threshold_level) = @_;

    if (defined $threshold_level) {
        $self->{-threshold_level} = $threshold_level;
    }
    return $self->{-threshold_level};
}

=head2 search_region_level

 Title   : search_region_level
 Usage   : $srlevel = $tc->search_region_level()
           or $tc->search_region_level($level);
 Function: Get/set the search region level
 Returns : An integer search region level
 Args    : Optionally an integer search region level

=cut

sub search_region_level
{
    my ($self, $search_region_level) = @_;

    if (defined $search_region_level) {
        $self->{-search_region_level} = $search_region_level;
    }
    return $self->{-search_region_level};
}

=head2 count

 Title   : count
 Usage   : $count = $tc->count()
           or $tc->search_count($count);
 Function: Get/set the TFBS count
 Returns : An integer TFBS count
 Args    : Optionally an integer TFBS count

=cut

sub count
{
    my ($self, $count) = @_;

    if (defined $count) {
        $self->{-count} = $count;
    }
    return $self->{-count};
}

1;
