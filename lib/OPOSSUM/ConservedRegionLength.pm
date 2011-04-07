=head1 NAME

OPOSSUM::ConservedRegionLength - ConservedRegionLength object
(conserved_region_lengths DB record)

=head1 DESCRIPTION

A ConservedRegionLength object models a record retrieved from the
conserved_region_lengths table of the oPOSSUM DB. The ConservedRegionLength
object contains the sum of the lengths of all the conserved regions for a
particular gene at a given level of conservation and amount of search region.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ConservedRegionLength;

use strict;
use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);

=head2 new

 Title   : new
 Usage   : $crl = OPOSSUM::ConservedRegionLength->new(
                -gene_id                => 28479,
                -conservation_level     => 1,
                -search_region_level    => 1,
                -length                 => 481
           );

 Function: Construct a new ConservedRegionLength object
 Returns : a new OPOSSUM::ConservedRegionLength object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 gene_id

 Title   : gene_id
 Usage   : $gid = $crl->gene_id()
           or $crl->gene_id($gid)
 Function: Get/set the gene ID of the conserved region
 Returns : A gene ID.
 Args    : Optionally a gene ID.

=cut

sub gene_id
{
    my ($self, $gid) = @_;

    if (defined $gid) {
        $self->{-gene_id} = $gid;
    }
    return $self->{-gene_id};
}

=head2 conservation_level

 Title   : conservation_level
 Usage   : $clevel = $crl->conservation_level()
           or $crl->conservation_level($level);
 Function: Get/set the conservation level of the conserved region
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

=head2 search_region_level

 Title   : search_region_level
 Usage   : $srlevel = $crl->search_region_level()
           or $crl->search_region_level($level);
 Function: Get/set the search region level of the conserved region
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

=head2 length

 Title   : length
 Usage   : $length = $crl->length() or $crl->length($length);
 Function: Get/set the length of the conserved region
 Returns : Integer length of this conserved region
 Args    : Optionally the integer length of this conserved region

=cut

sub length
{
    my ($self, $length) = @_;

    if (defined $length) {
        $self->{-length} = $length;
    }
    return $self->{-length};
}

1;
