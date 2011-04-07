=head1 NAME

OPOSSUM::ExternalGeneIDType - ExternalGeneIDType object (external_gene_id_types
DB record)

=head1 DESCRIPTION

An ExternalGeneIDType models a record retrieved from the external_gene_id_types
table of the oPOSSUM DB.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut
package OPOSSUM::ExternalGeneIDType;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $xid_type = OPOSSUM::ExternalGeneIDType->new(
			    -id_type	    => 5,
			    -name		    => 'UniProt/TrEMBL',
                -dblink_name    => 'UniProt/SPTREMBL'
           );

 Function: Construct a new ExternalGeneIDType object
 Returns : a new OPOSSUM::ExternalGeneIDType object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
		    %args
		}, ref $class || $class;

    return $self;
}

=head2 id_type

 Title   : id_type
 Usage   : $id_type = $xid_type->id_type() or $xid_type->id_type($id_type);

 Function: Get/set the ID type.
 Returns : A string.
 Args    : None or an id type string.

=cut

sub id_type
{
    my ($self, $id_type) = @_;

    if (defined $id_type) {
        $self->{-id_type} = $id_type;
    }
    return $self->{-id_type};
}

=head2 name

 Title   : name
 Usage   : $name = $xid_type->name() or $xid_type->name($name);

 Function: Get/set the name of this ID type.
 Returns : A string.
 Args    : None or a new name.

=cut

sub name
{
    my ($self, $name) = @_;

    if (defined $name) {
        $self->{-name} = $name;
    }
    return $self->{-name};
}

=head2 dblink_name

 Title   : dblink_name
 Usage   : $dblink_name = $xid_type->dblink_name()
           or $xid_type->dblink_name($dblink_name);

 Function: Get/set the dblink_name of this ID type. This is the external
           DB link name used by Ensembl when using e.g.,
           $ensembl_gene->get_all_DBLinks().
 Returns : A string.
 Args    : None or a new dblink_name.

=cut

sub dblink_name
{
    my ($self, $dblink_name) = @_;

    if (defined $dblink_name) {
        $self->{-dblink_name} = $dblink_name;
    }
    return $self->{-dblink_name};
}

1;
