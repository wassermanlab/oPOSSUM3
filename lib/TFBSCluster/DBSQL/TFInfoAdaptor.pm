=head1 NAME

TFBSCluster::DBSQL::TFInfoAdaptor - Adaptor for MySQL queries to retrieve
and store TF infomation.

=head1 SYNOPSIS

 $tfia = $db_adaptor->get_TFInfoAdaptor();

=head1 DESCRIPTION

The tf_info table contains records which store information about the TFBS
position weight matrices.

=head1 AUTHOR

 Andrew Kwon
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package TFBSCluster::DBSQL::TFInfoAdaptor;

use strict;

use Carp;

use TFBSCluster::DBSQL::BaseAdaptor;
use TFBSCluster::TFInfo;

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

    my $sql =
        qq{select tf_id, source, collection, external_id, name, class, family,
		tax_group, width, ic
        from tf_info};

    $sql .= " $where" if $where;

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch tf_info:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch tf_info:\n$sql\n" . $self->errstr;
        return;
    }

	my @tf_info_list;
    my $row_count = 0;
    while (my @row = $sth->fetchrow_array)
    {
        $row_count++;
        
        my $tfi = TFBSCluster::TFInfo->new(
            -tf_id			=> $row[0],
            -source			=> $row[1],
			-collection		=> $row[2],
            -external_id	=> $row[3],
			-name			=> $row[4],
			-class			=> $row[5],
			-family			=> $row[6],
			-tax_group		=> $row[7],
			-width			=> $row[8],
			-ic				=> $row[9]
        );

        push @tf_info_list, $tfi;
    }
    
    $sth->finish;

    return @tf_info_list if wantarray();

    #if ($row_count) {
    #    if ($row_count > 1) {
    #        return \@tf_info_list;
    #    }
    #    return $tf_info_list[0];
    #}

    return @tf_info_list ? \@tf_info_list : undef;
}

=head2 fetch_by_tf_id

 Title    : fetch_by_tf_id
 Usage    : $tfi = $tfia->fetch_by_tf_id($id);
 Function : Fetch the TFInfo object with the given tf_id.
 Returns  : An TFBSCluster::TFInfo object.
 Args	  : An Integer.

=cut

sub fetch_by_tf_id
{
    my ($self, $tf_id) = @_;
    
    my $where = "tf_id = '$tf_id'";
    
    my $tf_info_list = $self->fetch_where($where);

    if ($tf_info_list && $tf_info_list->[0]) {
        return $tf_info_list->[0];
    }

    return undef;
}

=head2 fetch_by_external_id

 Title    : fetch_by_external_id
 Usage    : $tfi = $tfia->fetch_by_external_id($jid);
 Function : Fetch TFInfo objects with the specified external DB ID.
 Returns  : An TFBSCluster::TFInfo object.
 Args	  : An external TF DB ID string.

=cut

sub fetch_by_external_id
{
    my ($self, $id) = @_;

    my $where = "external_id = '$id'";

    my $tf_info_list = $self->fetch_where($where);

    if ($tf_info_list && $tf_info_list->[0]) {
        return $tf_info_list->[0];
    }

    return undef;
}

=head2 fetch_by_name

 Title    : fetch_by_name
 Usage    : $tfis = $tfia->fetch_by_name($name);
 Function : Fetch the TFInfo object with the specified TF name.
 Returns  : A listref of TFBSCluster::TFInfo objects.
 Args	  : A TF name string

=cut

sub fetch_by_name
{
    my ($self, $name) = @_;

    my $where = "name = $name";

    return $self->fetch_where($where);
}

=head2 fetch_by_source

 Title    : fetch_by_source
 Usage    : $tfi = $tfia->fetch_by_($);
 Function : Fetch the TFInfo objects with the given source.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A source string.

=cut

sub fetch_by_source
{
    my ($self, $tf_source) = @_;
    
    my $where = "source = '$tf_source'";
    
    return $self->fetch_where($where);
}

=head2 fetch_by_collection

 Title    : fetch_by_collection
 Usage    : $tfi = $tfia->fetch_by_collection($collection);
 Function : Fetch the TFInfo objects with the given collection.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A collection string.

=cut

sub fetch_by_collection
{
    my ($self, $tf_collection) = @_;
    
    my $where = "collection = '$tf_collection'";
    
    return $self->fetch_where($where);
}

=head2 fetch_by_class

 Title    : fetch_by_class
 Usage    : $tfi = $tfia->fetch_by_class($cl);
 Function : Fetch the TFInfo objects with the given class.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A class string.

=cut

sub fetch_by_class
{
    my ($self, $tf_class) = @_;
    
    my $where = "class = '$tf_class'";
    
    return $self->fetch_where($where);
}

=head2 fetch_by_family

 Title    : fetch_by_family
 Usage    : $tfi = $tfia->fetch_by_family($fl);
 Function : Fetch the TFInfo objects with the given family.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A family string.

=cut

sub fetch_by_family
{
    my ($self, $tf_family) = @_;
    
    my $where = "family = '$tf_family'";
    
    return $self->fetch_where($where);
}

=head2 fetch_by_tax_group

 Title    : fetch_by_tax_group
 Usage    : $tfi = $tfia->fetch_by_tax_group($tax_group);
 Function : Fetch the TFInfo objects with the given tax_group.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A tax_group string.

=cut

sub fetch_by_tax_group
{
    my ($self, $tf_tax_group) = @_;
    
    my $where = "tax_group = '$tf_tax_group'";
    
    return $self->fetch_where($where);
}

=head2 fetch_by_min_width

 Title    : fetch_by_min_width
 Usage    : $tfi = $tfia->fetch_by_min_width($min_width);
 Function : Fetch the TFInfo objects with the given min_width.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A min_width string.

=cut

sub fetch_by_min_width
{
    my ($self, $tf_min_width) = @_;
    
    my $where = "width >= '$tf_min_width'";
    
    return $self->fetch_where($where);
}

=head2 fetch_by_min_ic

 Title    : fetch_by_min_ic
 Usage    : $tfi = $tfia->fetch_by_min_ic($min_ic);
 Function : Fetch the TFInfo objects with the given min_ic.
 Returns  : An array ref of TFBSCluster::TFInfo object.
 Args	  : A min_ic string.

=cut

sub fetch_by_min_ic
{
    my ($self, $tf_min_ic) = @_;
    
    my $where = "ic >= '$tf_min_ic'";
    
    return $self->fetch_where($where);
}

=head2 store

 Title   : store
 Usage   : $id = $tfia->store($tfi);
 Function: Store the given TFInfo object in the database.
 Args    : The TFBSCluster::TFInfo object to store.
 Returns : A database ID of the newly stored tf_cluster record.

=cut

sub store
{
    my ($self, $tfi) = @_;

    if (!ref $tfi || !$tfi->isa('TFBSCluster::TFInfo')) {
        carp "Not an TFBSCluster::TFInfo object";
        return;
    }

    my $sql = qq{
        insert into tf_info
            (tf_id, source, collection, external_id, name, class, family,
			tax_group, width, ic)
        values (?,?,?,?,?,?,?,?,?,?)
    };
    
	my $sth = $self->prepare($sql);
	
    if (!$sth) {
        carp "Error preparing insert tf_info statement\n" . $self->errstr;
        return;
    }
    
    if (!$sth->execute($tfi->id, $tfi->source, $tfi->collection, $tfi->external_id,
					   $tfi->name, $tfi->class, $tfi->family, $tfi->tax_group,
					   $tfi->width, $tfi->ic))
    {
        carp "Error inserting tf_info record\n" . $self->errstr;
        return;
    }
    $sth->finish;
	
	return $sth->{'mysql_insertid'};
}

1;
