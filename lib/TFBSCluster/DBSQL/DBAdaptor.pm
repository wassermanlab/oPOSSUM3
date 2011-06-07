=head1 NAME

0POSSUM_CLUSTER::DBSQL::DBAdaptor

=head1 DESCRIPTION

This object represents a database. Once created you can retrieve database
adaptors specific to various database objects that allow the retrieval and
creation of objects from the database.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=cut

package TFBSCluster::DBSQL::DBAdaptor;

use strict;

use Carp;
use TFBSCluster::DBSQL::DBConnection;

=head2 new

 Title    : new
 Usage    : $db_adaptor = TFBSCluster::DBSQL::DBAdaptor->new(
                -dbconn => $dbc
            );

            OR

            $db_adaptor = TFBSCluster::DBSQL::DBAdaptor->new(
                -user	=> 'opossum_r',
                -host	=> 'localhost',
                -dbname	=> 'oPOSSUM'
            );

 Function : Construct a new DBAdaptor object.
 Returns  : A new TFBSCluster::DBSQL::DBAdaptor object
 Args	  : Either a DBConnection or args for a DBConnection (passed through)

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {}, ref $class || $class;

    my $con = $args{-dbconn};
    if (defined $con) {
        $self->dbc($con);
    } else {
        $self->dbc(TFBSCluster::DBSQL::DBConnection->new(%args));
    }

    return $self;
}

=head2 dbc

  Title     : dbc
  Usage     : $dbc = $dba->dbc();
  Function  : Get/set DBConnection.
  Returns   : An TFBSCluster::DBSQL::DBConnection
  Args      : Optional Bio::EnsEMBL::DBSQL::DBConnection

=cut

sub dbc
{
  my $self = shift;

    if (@_) {
        my $arg = shift;

        if (defined($arg)) {
            if (!$arg->isa('TFBSCluster::DBSQL::DBConnection')) {
                carp "not a DBConnection\n";
            }
        }

        $self->{_dbc} = $arg;
    }

    return $self->{_dbc};
}


=head2 get_TFClusterAdaptor

 Title    : get_TFClusterAdaptor
 Usage    : $dbia = $db_adaptor->get_TFClusterAdaptor();
 Function : Construct a new TFBSCluster::DBSQL::TFClusterAdaptor object.
 Returns  : A new TFBSCluster::DBSQL::TFClusterAdaptor object
 Args	  : None.

=cut

sub get_TFClusterAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("TFBSCluster::DBSQL::TFClusterAdaptor");
}


=head2 get_TFInfoAdaptor

 Title    : get_TFInfoAdaptor
 Usage    : $dbia = $db_adaptor->get_TFInfoAdaptor();
 Function : Construct a new TFBSCluster::DBSQL::TFInfoAdaptor object.
 Returns  : A new TFBSCluster::DBSQL::TFInfoAdaptor object
 Args	  : None.

=cut

sub get_TFInfoAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("TFBSCluster::DBSQL::TFInfoAdaptor");
}


=head2 get_DBInfoAdaptor

 Title    : get_DBInfoAdaptor
 Usage    : $dbia = $db_adaptor->get_DBInfoAdaptor();
 Function : Construct a new TFBSCluster::DBSQL::DBInfoAdaptor object.
 Returns  : A new TFBSCluster::DBSQL::DInfoAdaptor object
 Args	  : None.

=cut

sub get_DBInfoAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("TFBSCluster::DBSQL::DBInfoAdaptor");
}

=head2 _get_adaptor

 Title    : _get_adaptor
 Usage    : $adpator = $self->_get_adaptor("full::adaptor::name");
 Function : Used by subclasses to obtain adaptor objects from this
            database adaptor using the fully qualified module name
            of the adaptor. If the adaptor has not been retrieved before
            it is created, otherwise it is retrieved from the adaptor
            cache.
 Returns  : Adaptor object.
 Args	  : Fully qualified adaptor module name,
            optional arguments to be passed to the adaptor constructor.

=cut

sub _get_adaptor
{
    my ($self, $module) = @_;

    my ($adaptor, $internal_name);
  
    #Create a private member variable name for the adaptor by replacing
    #:: with _
  
    $internal_name = $module;

    $internal_name =~ s/::/_/g;

    unless (defined $self->{'_adaptors'}{$internal_name}) {
        eval "require $module";
        
        if ($@) {
            carp "$module cannot be found.\nException $@\n";
            return undef;
        }
          
        $adaptor = "$module"->new($self);

        $self->{'_adaptors'}{$internal_name} = $adaptor;
    }

    return $self->{'_adaptors'}{$internal_name};
}

1;
