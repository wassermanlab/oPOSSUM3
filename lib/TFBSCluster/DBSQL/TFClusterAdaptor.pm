=head1 NAME

OPOSSUM_::DBSQL::TFClusterAdaptor - Adaptor for MySQL queries to retrieve
and store TF cluster infomation.

=head1 SYNOPSIS

 $tfca = $db_adaptor->get_TFClusterAdaptor();

=head1 DESCRIPTION

The tf_cluster table contains records which store information about the TFBS
position weight matrices.

=head1 AUTHOR

 Andrew Kwon
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut

package TFBSCluster::DBSQL::TFClusterAdaptor;

use strict;

use Carp;

use TFBSCluster::DBSQL::BaseAdaptor;
use TFBSCluster::TFClusterSet;
use TFBSCluster::TFCluster;

use vars '@ISA';
@ISA = qw(TFBSCluster::DBSQL::BaseAdaptor);

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

sub fetch_where
{
    my ($self, $where) = @_;
	
    if ($where && $where !~ /^\s*where /) {
        $where = "where $where";
    }

    my $sql = "select * from clusters";
    $sql .= " $where" if $where;

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch clusters:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch clusters:\n$sql\n" . $self->errstr;
        return;
    }
	
    my @tf_clusters;
    while (my @row = $sth->fetchrow_array) {
        my $tf_cluster = TFBSCluster::TFCluster->new(
            -adaptor	=> $self,
            -cluster_id => $row[0],
            -class	=> $row[1],
            -family	=> $row[2],
            -avg_width	=> $row[3]
        );
            
        my $tf_sql = "select tf_id from tfs where cluster_id = " . $row[0];
        my $tf_sth = $self->prepare($tf_sql);
        if (!$tf_sth) {
            carp "Error preparing fetch tfs:\n$tf_sql\n" . $self->errstr;
            return;
        }

        if (!$tf_sth->execute) {
            carp "Error executing fetch tfs:\n$tf_sql\n" . $self->errstr;
            return;
        }

        my @tf_ids;
        while (my ($tf_id) = $tf_sth->fetchrow_array) {
            push @tf_ids, $tf_id;
        }

        $tf_sth->finish;

        $tf_cluster->add_tf_ids(\@tf_ids);

        push @tf_clusters, $tf_cluster;
    }

    $sth->finish;

    #carp "Array size = " . scalar @tf_clusters . "\n";
	
    return @tf_clusters if wantarray();

	return @tf_clusters ? \@tf_clusters : undef;
}

=head2 fetch_all

 Title    : fetch_all
 Usage    : $cl = $tfca->fetch_all;
 Function : Fetch all TFCluster objects
 Returns  : A listref of TFCluster objects.
 Args	  : None

=cut

sub fetch_all
{
    my $self = shift;

    return $self->fetch_where();
}

=head2 fetch_by_cluster_id

 Title    : fetch_by_cluster_id
 Usage    : $cl = $tfca->fetch_by_cluster_id($id);
 Function : Fetch the TFCluster object with the specified ID.
 Returns  : An TFBSCluster::TFCluster object.
 Args	  : An integer TF cluster id

=cut

sub fetch_by_cluster_id
{
    my ($self, $id) = @_;

    my $where = "cluster_id = $id";

    my $clusters = $self->fetch_where($where);

    if ($clusters && $clusters->[0]) {
        return $clusters->[0];
    }

    return undef;
}

=head2 fetch_by_tf_id

 Title    : fetch_by_tf_id
 Usage    : $cl = $tfca->fetch_by_tf_id($tfid);
 Function : Fetch the TFCluster object with the specified TF ID.
 Returns  : An TFBSCluster::TFCluster object.
 Args	  : A string TF id

=cut

sub fetch_by_tf_id
{
    my ($self, $tfid) = @_;

	my $sql = "select cluster_id from tfs where tf_id = '$tfid'";
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch tfs:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch tfs:\n$sql\n" . $self->errstr;
        return;
    }

	my ($cluster_id) = $sth->fetchrow_array;
	if (!$cluster_id) {
		carp "No cluster associated with the given tf id $tfid\n";
		return;
	}
	
	return $self->fetch_by_cluster_id($cluster_id);
}

=head2 fetch_by_tf_id_list

 Title   : fetch_by_tf_id_list
 Usage   : $clusters = $cla->fetch_by_tf_id_list($id_list);
 Function: Fetch a list of cluster objects from the DB according to a list
           of TF IDs.
 Args    : A reference to a list of unique TF IDs.
 Returns : A reference to a list of TFBSCluster::TFCluster objects.

=cut

sub fetch_by_tf_id_list
{
    my ($self, $tf_ids) = @_;

    my $sql = qq{select distinct(cluster_id) from tfs where tf_id = ?};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching tfs\n" . $self->errstr;
        return;
    }

    my @clusters;
    foreach my $tfid (@$tf_ids) {
        if (!$sth->execute($tfid)) {
            carp "error fetching cluster_id for TF ID $tfid\n"
                . $self->errstr;
            return;
        }
        
		while (my ($cluster_id) = $sth->fetchrow_array) {
			my $cluster = $self->fetch_by_cluster_id($cluster_id);
			push @clusters, $cluster;
        }
    }
    $sth->finish;

    return @clusters ? \@clusters : undef;
}

=head2 fetch_by_tf_classes

 Title    : fetch_by_tf_classes
 Usage    : $clset = $tfca->fetch_by_tf_classes($classes);
 Function : Fetch the TFCluster objects of the specified TF classes.
 Returns  : A listref of TFBSCluster::TFCluster objects.
 Args	  : A string or listref of strings of tf classes

=cut

sub fetch_by_tf_classes
{
    my ($self, $classes) = @_;

    my $class_string;
    if (ref $classes eq 'ARRAY') {
        $class_string = join("','", @$classes);
        $class_string = "'" . $class_string . "'";
    } elsif (ref $classes eq 'SCALAR') {
        $class_string = "'$classes'";
    } else {
        carp "Provided argument is not a string or listref of strings\n";
        return;
    }
    
    my $where = "class in ($class_string)";

    return $self->fetch_where($where);
}

=head2 fetch_by_tf_families

 Title    : fetch_by_tf_families
 Usage    : $clset = $tfca->fetch_by_tf_families($families);
 Function : Fetch the TFCluster objects of the specified TF families.
 Returns  : A listref of TFBSCluster::TFCluster objects.
 Args	  : A string or listref of strings of tf families

=cut

sub fetch_by_tf_families
{
    my ($self, $families) = @_;

    my $family_string;
    if (ref $families eq 'ARRAY') {
        $family_string = join("','", @$families);
        $family_string = "'" . $family_string . "'";
    } elsif (ref $families eq 'SCALAR') {
        $family_string = "'$families'";
    } else {
        carp "Provided argument is not a string or listref of strings\n";
        return;
    }
    
    my $where = "family in ($family_string)";

    return $self->fetch_where($where);
}

=head2 fetch_cluster_ids

 Title   : fetch_cluster_ids
 Usage   : $cids = $ga->fetch_cluster_ids($where);
 Function: Fetch list of cluster IDs from the DB.
 Args    : Optionally a where clause.
 Returns : Reference to a list of internal cluster IDs. If no where clause
           is provided, returns all cluster IDs in the database.

=cut

sub fetch_cluster_ids
{
    my ($self, $where_clause) = @_;

    my $sql = "select cluster_id from clusters";
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching cluster IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching cluster IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_by_cluster_id_list

 Title   : fetch_by_cluster_id_list
 Usage   : $clusters = $cla->fetch_by_cluster_id_list($id_list);
 Function: Fetch a list of cluster objects from the DB according to a list
           of cluster IDs.
 Args    : A reference to a list of unique internal cluster IDs.
 Returns : A reference to a list of TFBSCluster::TFCluster objects.

=cut

sub fetch_by_cluster_id_list
{
    my ($self, $cluster_ids) = @_;

    my @clusters;
    foreach my $cid (@$cluster_ids) {
		my $cluster = $self->fetch_by_cluster_id($cid);
		push @clusters, $cluster;
    }

    return @clusters ? \@clusters : undef;
}

=head2 store

Forget this one for now.

=cut

1;
