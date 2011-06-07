
=head1 NAME

TFBSCluster::DBSQL::BaseAdaptor - Base Adaptor for DBSQL adaptors

=head1 DESCRITION

Base class for Adaptors in the oPOSSUM DBSQL. This is adapted from the
Bio::EnsEMBL::DBSQL::BaseAdaptor module.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package TFBSCluster::DBSQL::BaseAdaptor;

use strict;

use Carp;
use DBI qw(:sql_types);
use Data::Dumper;

require Exporter;
use vars qw(@ISA @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 new

 Title    : new
 Usage    : $adaptor = AdaptorInheritedFromBaseAdaptor->new($dba);
 Function : Construct a new BaseAdaptor object. This should not be called
            directly but from an class inherited from this one.
 Returns  : a new TFBSCluster::DBSQL::BaseAdaptor object
 Args	  : an TFBSCluster::DBSQL::DBAdaptor object

=cut

sub new
{
    my ($class, $dba) = @_;

    if (!defined $dba || !ref $dba) {
        carp "No DBAdaptor object for new adaptor";
        return;
    }

    my $self = bless {}, ref $class || $class;

    $self->db($dba);
    $self->dbc($dba->dbc);

    return $self;
}

=head2 prepare

 Title    : prepare
 Usage    : $sth = $adaptor->prepare($sql_statement);
 Function : Return a DBI statement handle from the adaptor.
 Returns  : A DBI statement handle.
 Args     : An SQL statement to be prepared by this adaptor's database.

=cut

sub prepare
{
    my ($self, $string) = @_;

    return $self->dbc->prepare($string);
}

=head2 errstr

 Title    : errstr
 Usage    : $err = $adaptor->errstr();
 Function : Return the last DBI error string.
 Returns  : A string.
 Args     : None.

=cut

sub errstr
{
    $_[0]->dbc->db_handle->errstr || "";
}

=head2 db

 Title    : db
 Usage    : $db = $adaptor->db() or $adaptor->db($dbobj);
 Function : Get/set the database object this adaptor is using.
 Returns  : An TFBSCluster::DBSQL::DBConnection object.
 Args     : None or an TFBSCluster::DBSQL::DBConnection object.

=cut

sub db
{
    my ($self, $db) = @_;

    if ($db) {
        $self->{-db} = $db;
    }

    return $self->{-db};
}

=head2 dbc

  Title     : dbc
  Usage     : $dbc = $adaptor->dbc();
  Function  : Get/Set the DatabaseConnection that this adaptor is using.
  Returns   : An TFBSCluster::DBSQL::DBConnection
  Args      : An optional TFBSCluster::DBSQL::DBConnection object.

=cut

sub dbc
{
    my ($self, $dbc) = @_;

    if (defined($dbc)) {
        $self->{'-dbc'} = $dbc;
    }

    return $self->{'-dbc'};
}

1;
