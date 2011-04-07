=head1 NAME

OPOSSUM::TFInfo - TFInfo object (tf_info DB record)

NOTE: This module is DEPRECATED! We will not longer store TF information
in the oPOSSUM database. Instead TF (matrix) information should be retrieved
directly from the JASPAR/PAZAR database.

=head1 DESCRIPTION

A TFInfo object models a record retrieved from the tf_info table of the
oPOSSUM DB. It stores information about transcription factor binding site
profiles.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::TFInfo;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $ti = OPOSSUM::TFInfo->new(
			    -id			    => '1',
			    -name		    => 'AGL3',
			    -source	        => 'JASPAR',
			    -external_id	=> 'MA0001',
			    -collection	    => 'CORE',
			    -class		    => 'Other Alpha-Helix',
			    -family		    => 'MADS',
			    -tax_group		=> 'plants',
			    -width		    => 10,
			    -ic			    => 10.5882);

 Function: Construct a new TFInfo object
 Returns : a new OPOSSUM::TFInfo object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        %args
    }, ref $class || $class;

    if ($args{-external_db}) {
        carp "deprecated tag '-external_db'; please use '-source' instead\n";

        if (!$args{-source}) {
            $self->{-source} = $args{-external_db};
        }
    }

    return $self;
}

=head2 id

 Title   : id
 Usage   : $id = $ti->id() or $ti->id($id);

 Function: Get/set the ID of the TFBS profile. This should be a unique
 	   identifier for this object within the implementation. If
	   the GenePair object was read from the oPOSSUM database,
	   this should be set to the value in the tf_id column.
 Returns : A string
 Args    : None or an id string

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
        $self->{-id} = $id;
    }
    return $self->{-id};
}

=head2 name

 Title   : name
 Usage   : $name = $ti->name() or $ti->name($name);

 Function: Get/set the name of the TFBS profile.
 Returns : A string
 Args    : None or an id string

=cut

sub name
{
    my ($self, $name) = @_;

    if (defined $name) {
        $self->{-name} = $name;
    }
    return $self->{-name};
}

=head2 external_id

 Title   : external_id
 Usage   : $id = $ti->external_id() or $ti->external_id($id);

 Function: Get/set the external ID of the TFBS profile. Within the context
 	   of the oPOSSUM database this is the unique ID of this profile
	   in the originating DB (i.e. JASPAR2).
 Returns : An external ID string
 Args    : None or an ID string

=cut

sub external_id
{
    my ($self, $external_id) = @_;

    if (defined $external_id) {
        $self->{-external_id} = $external_id;
    }
    return $self->{-external_id};
}

=head2 source

 Title   : source
 Usage   : $source = $ti->source() or $ti->source($db);

 Function: Get/set the source (formerly called external DB) of the TFBS
           profile (i.e. 'JASPAR').
 Returns : The TFBS profile source
 Args    : None or a new source

=cut

sub source
{
    my ($self, $source) = @_;

    if (defined $source) {
        $self->{-source} = $source;
    }

    return $self->{-source};
}

sub external_db
{
    my ($self, $db_name) = @_;

    carp "deprecated method 'external_db()'; please use 'source()' instead\n";

    return $self->source($db_name);
}

=head2 class

 Title   : class
 Usage   : $class = $ti->class() or $ti->class($class);

 Function: Get/set the class of the TFBS profile, e.g. 'Helix-Turn-Helix'
 Returns : A string
 Args    : None or a string

=cut

sub class
{
    my ($self, $class) = @_;

    if (defined $class) {
        $self->{-class} = $class;
    }
    return $self->{-class};
}

=head2 family

 Title   : family
 Usage   : $family = $ti->family() or $ti->family($family);

 Function: Get/set the family of the TFBS profile, e.g. 'Leucine Zipper'
 Returns : A string
 Args    : None or a string

=cut

sub family
{
    my ($self, $family) = @_;

    if (defined $family) {
        $self->{-family} = $family;
    }
    return $self->{-family};
}

=head2 collection

 Title   : collection
 Usage   : $collection = $ti->collection() or $ti->collection($collection);

 Function: Get/set the collection of the TFBS profile, e.g. 'CORE'
 Returns : A string
 Args    : None or a string

=cut

sub collection
{
    my ($self, $collection) = @_;

    if (defined $collection) {
        $self->{-collection} = $collection;
    }
    return $self->{-collection};
}

=head2 tax_group

 Title   : tax_group
 Usage   : $tax_group = $ti->tax_group() or $ti->tax_group($tax_group);

 Function: Get/set the taxonomic supergroup of this TFBS profile,
           e.g. 'vertebrates'
 Returns : A string
 Args    : None or a string

=cut

sub tax_group
{
    my ($self, $tax_group) = @_;

    if (defined $tax_group) {
        $self->{-tax_group} = $tax_group;
    }
    return $self->{-tax_group};
}

#
# Old method name
#
sub phylum
{
    my ($self, $tax_group) = @_;

    carp "Method 'phylum' is deprecated - please use 'tax_group' instead\n";

    return $self->tax_group($tax_group);
}

=head2 width

 Title   : width
 Usage   : $width = $ti->width() or $ti->width($width);

 Function: Get/set the width of the TFBS profile in nucleotides.
 Returns : An integer
 Args    : None or an integer

=cut

sub width
{
    my ($self, $width) = @_;

    if (defined $width) {
        $self->{-width} = $width;
    }
    return $self->{-width};
}

=head2 ic

 Title   : ic
 Usage   : $ic = $ti->ic() or $ti->ic($ic);

 Function: Get/set the information content of the TFBS profile. This is also
 	   known as specificity.
 Returns : A float
 Args    : None or a float

=cut

sub ic
{
    my ($self, $ic) = @_;

    if (defined $ic) {
        $self->{-ic} = $ic;
    }
    return $self->{-ic};
}

1;
