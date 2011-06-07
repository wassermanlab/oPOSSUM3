=head1 NAME

TFBSCluster::TFCluster - TFCluster object (tf_clusters DB record)

=head1 DESCRIPTION


=head1 AUTHOR

 Andrew Kwon
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut
package TFBSCluster::TFCluster;


use strict;

use Carp;
use TFBSCluster::DBObject;use TFBS::Matrix;

use vars qw(@ISA);
@ISA = qw(TFBSCluster::DBObject);

=head2 new

 Title   : new
 Usage   : $tfc = TFBSCluster::TFCluster->new(
				-cluster_id		=> 3,
				-tf_class => 'Other',
				-tf_family => 'Other',
				-avg_width => 10,
				-tf_ids => $tfids);

 Function: Construct a new TFCluster object
 Returns : a new TFBSCluster::TFCluster object

=cut

sub new
{
	my ($class, %args) = @_;

	my $self = bless {
			%args
		}, ref $class || $class;
	
	return $self;
}

=head2 id

 Title	: id
 Usage	: $id = $tc->id() or $tc->id($id);
 Function: Get/set the ID of the TFCluster.  This should be a unique 
 	identifier for this object within the implementation. If the TFCluster 
	object was read from the oPOSSUM DB, this should be set to the value in 
	the cluster_id column.
 Returns: An integer
 Args	: None or an id integer
=cut

sub id
{
	my ($self, $id) = @_;

	if (defined $id) {
		$self->{-cluster_id} = $id;
	}
	return $self->{-cluster_id};
}

=head2 name

 Title	: name
 Usage	: $name = $tc->name()
 Function: Get the name of the TFCluster. This is name is produced on the fly to
 provide the user with an informative description of the cluster, based on the
 cluster ID and TF classes.
 Returns: A string
 Args	: None
=cut

sub name
{
	my $self = shift;
	
	my $name = "C" . $self->id;
	$name .= " - " . $self->class;
	$name .= "::" . $self->family;

	return $name;
}

=head2 tf_ids

  Title   : tf_ids
  Usage   : $tfids = $tfc->tf_ids();
  Function: Gets the TF IDs associated with this cluster.
  Returns : A listref of TF IDs
  Args    : None

=cut

sub tf_ids
{
	my $self = shift;
	
	if (defined $self->{-tf_ids}) {
		return $self->{-tf_ids};
	}
	
	return;
}

=head2 add_tf_ids

  Title   : add_tf_ids
  Usage   : $tfc->add_tf_ids($tf_ids);
  Function: Adds the given listref of TF IDs to the internal TF ID list.
  Returns : None
  Args    : listref of TF IDs

=cut

sub add_tf_ids
{
	my ($self, $tf_ids) = @_;
	
	if (ref $tf_ids eq 'ARRAY') {
		push @{$self->{-tf_ids}}, @$tf_ids;
	} elsif (ref $tf_ids eq 'SCALAR') {
		push @{$self->{-tf_ids}}, $tf_ids;
	}
	
	return;
}

=head2 tf_class

 Title   : tf_class
 Usage   : $tf_class = $tc->tf_class();

 Function: Get the TF class for this cluster.
 Returns : A TF class string
 Args	: None or a TF class string

=cut

sub tf_class
{
	my ($self, $tf_class) = @_;
	
	if (defined $tf_class) {
		$self->{-class} = $tf_class;
	}
	
	return $self->{-class};
}

=head2 class

 Title   : class
 Usage   : $tf_class = $tc->tf_class();

 Function: Synonym of tf_class method.
 Returns : A TF class string
 Args	: None or a TF class string

=cut

sub class
{
	my ($self, $class) = @_;
	
	return $self->tf_class($class);
}

=head2 family

 Title   : family
 Usage   : $family = $tc->family();

 Function: Get the TF family for this cluster.
 Returns : A TF family string
 Args	: None

=cut

sub family
{
	my ($self, $family) = @_;
	
	if (defined $family) {
		$self->{-family} = $family;
	}
	
	return $self->{-family};
}

sub tf_family
{
	my ($self, $family) = @_;
	
	return $self->family($family);
}

=head2 rep_tf_id

 Title   : rep_tf_id
 Usage   : $rep_tf_id = $tc->rep_tf_id();

 Function: Returns a representative TF ID for the cluster (first one in the list for now)
 Returns : TF ID string
 Args	: None

=cut

sub rep_tf_id
{
	my $self = shift;
	
	return $self->{-tf_ids}->[0];
}

=head2 contains_tf_id

 Title   : contains_tf_id
 Usage   : $tf_found = $tc->contains_tf_id($tfid);

 Function: Checks whether the cluster contains this TF ID
 Returns : 0 or 1
 Args	: TF ID integer

=cut

sub contains_tf_id
{
	my ($self, $tf_id) = @_;
	
	my $id_string = join (',', @{$self->{-tf_ids}});
	
	return 0 if $tf_id !~ /$id_string/;
	
	return 1;
}



=head2 match_tf_class

 Title   : match_tf_class
 Usage   : $match = $tc->match_tf_class($tfclass);

 Function: Checks whether the cluster has TFs belonging to the specified class
 Returns : 0 or 1
 Args	: TF class string

=cut

sub match_tf_class
{
	my ($self, $tf_class) = @_;
	
	return 0 if !$self->class;
	
	return 1 if $tf_class eq $self->class;
	
	return 0;
}

=head2 match_tf_family

 Title   : match_tf_family
 Usage   : $match = $tc->match_tf_family($tffamily);

 Function: Checks whether the cluster has TFs belonging to the specified family
 Returns : 0 or 1
 Args	: TF family string

=cut

sub match_tf_family
{
	my ($self, $tf_family) = @_;
	
	return 0 if !$self->family;
	
	return 1 if $tf_family eq $self->family;
	
	return 0;
}

=head2 avg_width

 Title	: avg_width
 Usage	: $avg_width = $tfc->avg_width();
 Function : Return the avg. TFBs profile width of the matrices belonging to the cluster
 Returns  : A float
 Args	 : None

=cut

sub avg_width
{
	return $_[0]->{-avg_width} || 0;
}

=head2 size

 Title	: size
 Usage	: $size = $tfc->size();
 Function : Return the size of the set (number of TFBS matrices in the cluster)
 Returns  : An integer
 Args	 : None

=cut

sub size
{
	return scalar(@{$_[0]->{-tf_ids}}) || 0;
}

1;
