=head1 NAME

0POSSUM::DBSQL::DBAdaptor

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

package OPOSSUM::DBSQL::DBAdaptor;

use strict;

use Carp;
use OPOSSUM::DBSQL::DBConnection;

=head2 new

 Title    : new
 Usage    : $db_adaptor = OPOSSUM::DBSQL::DBAdaptor->new(
                -dbconn => $dbc
            );

            OR

            $db_adaptor = OPOSSUM::DBSQL::DBAdaptor->new(
                -user	=> 'opossum_r',
                -host	=> 'localhost',
                -dbname	=> 'oPOSSUM'
            );

 Function : Construct a new DBAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::DBAdaptor object
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
        $self->dbc(OPOSSUM::DBSQL::DBConnection->new(%args));
    }

    return $self;
}

=head2 dbc

  Title     : dbc
  Usage     : $dbc = $dba->dbc();
  Function  : Get/set DBConnection.
  Returns   : An OPOSSUM::DBSQL::DBConnection
  Args      : Optional Bio::EnsEMBL::DBSQL::DBConnection

=cut

sub dbc
{
  my $self = shift;

    if (@_) {
        my $arg = shift;

        if (defined($arg)) {
            if (!$arg->isa('OPOSSUM::DBSQL::DBConnection')) {
                carp "not a DBConnection\n";
            }
        }

        $self->{_dbc} = $arg;
    }

    return $self->{_dbc};
}

=head2 get_GeneAdaptor

 Title    : get_GeneAdaptor
 Usage    : $gpa = $db_adaptor->get_GeneAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::GeneAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::GeneAdaptor object
 Args	  : None.

=cut

sub get_GeneAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::GeneAdaptor");
}

=head2 get_OperonAdaptor

 Title    : get_OperonAdaptor
 Usage    : $opa = $db_adaptor->get_OperonAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::OperonAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::OperonAdaptor object
 Args	  : None.

=cut

sub get_OperonAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::OperonAdaptor");
}

=head2 get_PromoterAdaptor

 Title    : get_PromoterAdaptor
 Usage    : $ppa = $db_adaptor->get_PromoterAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::PromoterAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::PromoterAdaptor object
 Args	  : None.

=cut

sub get_PromoterAdaptor
{
    my ($self) = @_;
 
    return $self->_get_adaptor("OPOSSUM::DBSQL::PromoterAdaptor");
}

=head2 get_SequenceAdaptor

 Title    : get_SequenceAdaptor
 Usage    : $sa = $db_adaptor->get_SequenceAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::SequenceAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::SequenceAdaptor object
 Args	  : None.

=cut

sub get_SequenceAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::SequenceAdaptor");
}

=head2 get_ConservedRegionAdaptor

 Title    : get_ConservedRegionAdaptor
 Usage    : $cra = $db_adaptor->get_ConservedRegionAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservedRegionAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ConservedRegionAdaptor object
 Args	  : None.

=cut

sub get_ConservedRegionAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservedRegionAdaptor");
}

=head2 get_ConservedTFBSAdaptor

 Title    : get_ConservedTFBSAdaptor
 Usage    : $ctfsa = $db_adaptor->get_ConservedTFBSAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservedTFBSAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ConservedTFBSAdaptor object
 Args	  : None.

=cut

sub get_ConservedTFBSAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservedTFBSAdaptor");
}

=head2 get_ConservedTFBSDimerAdaptor

 Title    : get_ConservedTFBSDimerAdaptor
 Usage    : $ctfsa = $db_adaptor->get_ConservedTFBSDimerAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor object
 Args	  : None.

=cut

sub get_ConservedTFBSDimerAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor");
}

=head2 get_ConservedNHRSiteAdaptor

 Title    : get_ConservedNHRSiteAdaptor
 Usage    : $ctfsa = $db_adaptor->get_ConservedNHRSiteAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservedNHRSiteAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ConservedNHRSiteAdaptor object
 Args	  : None.

=cut

sub get_ConservedNHRSiteAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservedNHRSiteAdaptor");
}

=head2 get_TFInfoAdaptor

 Title    : get_TFInfoAdaptor
 Usage    : $ctia = $db_adaptor->get_TFInfoAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::TFInfoAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::TFInfoAdaptor object
 Args	  : None.

=cut

sub get_TFInfoAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::TFInfoAdaptor");
}

=head2 get_DBInfoAdaptor

 Title    : get_DBInfoAdaptor
 Usage    : $dbia = $db_adaptor->get_DBInfoAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::DBInfoAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::DBInfoAdaptor object
 Args	  : None.

=cut

sub get_DBInfoAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::DBInfoAdaptor");
}

=head2 get_TFBSCountAdaptor

 Title    : get_TFBSCountAdaptor
 Usage    : $tca = $db_adaptor->get_TFBSCountAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::TFBSCountAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::TFBSCountAdaptor object
 Args	  : None.

=cut

sub get_TFBSCountAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::TFBSCountAdaptor");
}

=head2 get_TFBSClusterCountAdaptor

 Title    : get_TFBSClusterCountAdaptor
 Usage    : $tca = $db_adaptor->get_TFBSClusterCountAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::TFBSClusterCountAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::TFBSClusterCountAdaptor object
 Args	  : None.

=cut

sub get_TFBSClusterCountAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::TFBSClusterCountAdaptor");
}

=head2 get_AnalysisCountsAdaptor

 Title    : get_AnalysisCountsAdaptor
 Usage    : $aca = $db_adaptor->get_AnalysisCountsAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::Analysis::CountsAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::Analysis::CountsAdaptor object
 Args	  : None.

=cut

sub get_AnalysisCountsAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::Analysis::CountsAdaptor");
}

=head2 get_AnalysisClusterCountsAdaptor

 Title    : get_AnalysisClusterCountsAdaptor
 Usage    : $aca = $db_adaptor->get_AnalysisClusterCountsAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor object
 Args	  : None.

=cut

sub get_AnalysisClusterCountsAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor");
}

=head2 get_ConservedRegionLengthAdaptor

 Title    : get_ConservedRegionLengthAdaptor
 Usage    : $crla = $db_adaptor->get_ConservedRegionLengthAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservedRegionLengthAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ConservedRegionLengthAdaptor object
 Args	  : None.

=cut

sub get_ConservedRegionLengthAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservedRegionLengthAdaptor");
}

=head2 get_ConservationLevelAdaptor

 Title    : get_ConservationLevelAdaptor
 Usage    : $cla = $db_adaptor->get_ConservationLevelAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservationLevelAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ConservationLevelAdaptor object
 Args	  : None.

=cut

sub get_ConservationLevelAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservationLevelAdaptor");
}

=head2 get_SearchRegionLevelAdaptor

 Title    : get_SearchRegionLevelAdaptor
 Usage    : $srla = $db_adaptor->get_SearchRegionLevelAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::SearchRegionLevelAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::SearchRegionLevelAdaptor object
 Args	  : None.

=cut

sub get_SearchRegionLevelAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::SearchRegionLevelAdaptor");
}

=head2 get_ThresholdLevelAdaptor

 Title    : get_ThresholdLevelAdaptor
 Usage    : $tla = $db_adaptor->get_ThresholdLevelAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ThresholdLevelAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ThresholdLevelAdaptor object
 Args	  : None.

=cut

sub get_ThresholdLevelAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ThresholdLevelAdaptor");
}

=head2 get_ExternalGeneIDAdaptor

 Title    : get_ExternalGeneIDAdaptor
 Usage    : $xgia = $db_adaptor->get_ExternalGeneIDAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ExternalGeneIDAdaptor object.
 Returns  : A new OPOSSUM::DBSQL::ExternalGeneIDAdaptor object
 Args	  : None.

=cut

sub get_ExternalGeneIDAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ExternalGeneIDAdaptor");
}

=head2 get_ExternalGeneIDTypeAdaptor

 Title    : get_ExternalGeneIDTypeAdaptor
 Usage    : $xgita = $db_adaptor->get_ExternalGeneIDTypeAdaptor
 Function : Construct a new OPOSSUM::DBSQL::ExternalGeneIDTypeAdaptor
 Returns  : A new OPOSSUM::DBSQL::ExternalGeneIDTypeAdaptor
 Args	  : None.

=cut

sub get_ExternalGeneIDTypeAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ExternalGeneIDTypeAdaptor");
}

=head2 get_ExonAdaptor

 Title    : get_ExonAdaptor
 Usage    : $ea = $db_adaptor->get_ExonAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ExonAdaptor
            object.
 Returns  : A new OPOSSUM::DBSQL::ExonAdaptor object
 Args	  : None.

=cut

sub get_ExonAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ExonAdaptor");
}

=head2 get_ConservationAnalysisAdaptor

 Title    : get_ConservationAnalysisAdaptor
 Usage    : $caa = $db_adaptor->get_ConservationAnalysisAdaptor();
 Function : Construct a new OPOSSUM::DBSQL::ConservationAnalysisAdaptor
            object.
 Returns  : A new OPOSSUM::DBSQL::ConservationAnalysisAdaptor object
 Args	  : None.

=cut

sub get_ConservationAnalysisAdaptor
{
    my ($self) = @_;
    
    return $self->_get_adaptor("OPOSSUM::DBSQL::ConservationAnalysisAdaptor");
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
