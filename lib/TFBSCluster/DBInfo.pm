=head1 NAME

TFBSCluster::DBInfo - DBInfo object (db_info DB record)

=head1 DESCRIPTION

A DBInfo object models the (single) record contained in the db_info
table of the oPOSSUM DB. The DBInfo object contains information about
how the oPOSSUM database was built, including the databases and software
versions which were used.

=head1 AUTHOR

 David Arenillas, modified by Andrew Kwon
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut

package TFBSCluster::DBInfo;

use strict;
use Carp;
use TFBSCluster::DBObject;

use vars qw(@ISA);

@ISA = qw(TFBSCluster::DBObject);

=head2 new

 Title   : new
 Usage   : $db_info = TFBSCluster::DBInfo->new(
                -build_date	        => '2010/08/31 17:23:06',
				-jaspar_db			=> 'JASPAR_2010',
			    -collections	    => 'CORE,PBM,PENDING'
			    -tax_groups		    => 'vertebrates,insects,nematodes,chordates',
                -min_ic             => 7,
                -cluster_threshold  => 1.8,
			    -radius_margin		=> 0.05,
				-comparison_method	=> 'mat_align',
			    -cluster_method		=> 'neighbour_search'
           );

 Function: Construct a new DBInfo object
 Returns : a new TFBSCluster::DBInfo object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 build_date

 Title   : build_date
 Usage   : $date = $dbbi->build_date() or $dbbi->build_date($date);

 Function: Get/set the DB date.
 Returns : A string.
 Args    : None or a new date.

=cut

sub build_date
{
    my ($self, $build_date) = @_;

    if ($build_date) {
        $self->{-build_date} = $build_date;
    }

    return $self->{-build_date};
}

=head2 jaspar_db

 Title   : jaspar_db
 Usage   : $jaspar_db = $dbbi->jaspar_db()
           or $dbbi->jaspar_db($jaspar_db);

 Function: Get/set the JASPAR DB which this cluster DB is based on.
 Returns : An jaspar_db string.
 Args    : None or a new jaspar_db.

=cut

sub jaspar_db
{
    my ($self, $jaspar_db) = @_;

    if (defined $jaspar_db) {
        $self->{-jaspar_db} = $jaspar_db;
    }

    return $self->{-jaspar_db};
}

=head2 collections

 Title   : collections
 Usage   : $collections = $dbbi->collections()
           or $dbbi->collections($collections);

 Function: Get/set the JASPAR collections that were used to build the clusters
            for this DB. The list should be made of those collections in the
            JASPAR DB used: i.e. 'CORE','PBM','PENDING' etc.
 Returns : A collections string.
 Args    : None or a new collections list string.

=cut

sub collections
{
    my ($self, $collections) = @_;

    if (defined $collections) {
        $self->{-collections} = $collections;
    }

    return $self->{-collections};
}

=head2 tax_groups

 Title   : tax_groups
 Usage   : $tax_groups = $dbbi->tax_groups()
           or $dbbi->tax_groups($tax_groups);

 Function: Get/set the taxonomic supergroups that were used for clustering.
           This should be the tax groups from the JASPAR CORE
           profile collection, e.g.: 'vertebrates', 'plants', 'nematodes',
           'fungi' etc. Note the pluralization of the name.
 Returns : A tax group string.
 Args    : None or a new tax group list string.

=cut

sub tax_groups
{
    my ($self, $tax_groups) = @_;

    if (defined $tax_groups) {
        $self->{-tax_groups} = $tax_groups;
    }

    return $self->{-tax_groups};
}

=head2 min_ic

 Title   : min_ic
 Usage   : $min_ic = $dbbi->min_ic() or $dbbi->min_ic($min_ic);

 Function: Get/set the name of the min_ic
 Returns : A float
 Args    : None or a new min_ic

=cut

sub min_ic
{
    my ($self, $min_ic) = @_;

    if ($min_ic) {
        $self->{-min_ic} = $min_ic;
    }

    return $self->{-min_ic};
}

=head2 cluster_threshold

 Title   : cluster_threshold
 Usage   : $th = $dbbi->cluster_threshold() or $dbbi->cluster_threshold($th);

 Function: Get/set the name of the cluster_threshold
 Returns : A float
 Args    : None or a new cluster_threshold

=cut

sub cluster_threshold
{
    my ($self, $th) = @_;

    if ($th) {
        $self->{-cluster_threshold} = $th;
    }

    return $self->{-cluster_threshold};
}

=head2 radius_margin

 Title   : radius_margin
 Usage   : $margin = $dbbi->radius_margin() or $dbbi->radius_margin($margin);

 Function: Get/set the name of the radius_margin
 Returns : A float
 Args    : None or a new radius_margin

=cut

sub radius_margin
{
    my ($self, $margin) = @_;

    if ($margin) {
        $self->{-radius_margin} = $margin;
    }

    return $self->{-radius_margin};
}

=head2 comparison_method

 Title   : comparison_method
 Usage   : $cm = $dbbi->comparison_method() or $dbbi->comparison_method($cm);

 Function: Get/set the name of the comparison_method
 Returns : A string
 Args    : None or a new comparison_method

=cut

sub comparison_method
{
    my ($self, $cm) = @_;

    if ($cm) {
        $self->{-comparison_method} = $cm;
    }

    return $self->{-comparison_method};
}

=head2 cluster_method

 Title   : cluster_method
 Usage   : $cm = $dbbi->cluster_method() or $dbbi->cluster_method($cm);

 Function: Get/set the name of the cluster_method
 Returns : A string
 Args    : None or a new cluster_method

=cut

sub cluster_method
{
    my ($self, $cm) = @_;

    if ($cm) {
        $self->{-cluster_method} = $cm;
    }

    return $self->{-cluster_method};
}


1;
