=head1 NAME

OPOSSUM::TFBSClusterCount - TFBSClusterCount object (tfbs_cluster_counts DB record)

=head1 DESCRIPTION

A TFBSClusterCount object models a record retrieved from the tfbs_cluster_counts
table of the oPOSSUM DB. It contains the number of times a given TFBS cluster was
detected for a given gene at a given level of conservation, PWM score threshold
and search region.

=head1 AUTHOR

 Andrew Kwon, modifying code by Dave Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::TFBSClusterCount;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $tc = OPOSSUM::TFBSClusterCount->new(
				-gene_id                => $gid,
				-cluster_id             => $cid,
				-conservation_level     => $clevel,
				-threshold_level        => $tlevel,
				-search_region_level    => $srlevel,
				-count                  => $count,
				-sum_length				=> $length);

 Function: Construct a new TFBSClusterCount object
 Returns : a new OPOSSUM::TFBSClusterCount object

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

=head2 cluster_id

 Title   : cluster_id
 Usage   : $id = $tc->cluster_id() or $tc->cluster_id($tfid);
 Function: Get/set the TFCluster ID
 Returns : A TFCluster ID
 Args    : Optionally a TFCluster ID

=cut

sub cluster_id
{
    my ($self, $cluster_id) = @_;

    if (defined $cluster_id) {
        $self->{-cluster_id} = $cluster_id;
    }
    return $self->{-cluster_id};
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
 Function: Get/set the TFBS cluster count
 Returns : An integer TFBS cluster count
 Args    : Optionally an integer TFBS cluster count

=cut

sub count
{
    my ($self, $count) = @_;

    if (defined $count) {
        $self->{-count} = $count;
    }
    return $self->{-count};
}

=head2 sum_length

 Title   : sum_length
 Usage   : $sum_length = $tc->sum_length()
           or $tc->search_sum_length($sum_length);
 Function: Get/set the sum of the bases covered by TFBS cluster 
 Returns : An integer TFBS cluster sum_length
 Args    : Optionally an integer TFBS cluster sum_length

=cut

sub sum_length
{
    my ($self, $sum_length) = @_;

    if (defined $sum_length) {
        $self->{-sum_length} = $sum_length;
    }
    return $self->{-sum_length};
}

1;
