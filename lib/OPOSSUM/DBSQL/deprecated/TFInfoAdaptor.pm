=head1 NAME

OPOSSUM::DBSQL::TFInfoAdaptor - Adaptor for MySQL queries to retrieve
and store TF infomation.

NOTE: This module is DEPRECATED! We will not longer store TF information
in the oPOSSUM databse. Instead TF (matrix) information should be retrieved
directly from the appropriate JASPAR/PAZAR database.

=head1 SYNOPSIS

 $tfia = $db_adaptor->get_TFInfoAdaptor();

 #
 # Fetch a TF info list (listref of OPOSSUM::TFInfo objects) by source(s),
 # collection(s) and tax group(s), ordered by tax group and name.
 #
 my $tf_info_list = $tfia->fetch_list(
                        -sources        => 'JASPAR',
                        -collections    => 'CORE',
                        -tax_groups     => ['vertebrates', 'insects'],
                        -order_by       => ['tax_group', 'name']
                    );

 #
 # Fetch a TF info set (OPOSSUM::TFInfoSet) by names
 #
 my $tf_info_set  = $tfia->fetch_set(
                        -names  => ['RUNX1',
                                    'TFAP2A',
                                    'Arnt',
                                    'Arnt::Ahr',
                                    'Ar']
                    );

=head1 DESCRIPTION

The tf_info table contains records which store information about the TFBS
position weight matrices.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::TFInfoAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::TFInfoSet;
use OPOSSUM::TFInfo;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_list

 Title    : fetch_list
 Usage    : $tf_info_list = $tfia->fetch_list(
                -sources        => $sources,
                -collections    => $collections,
                -tax_groups     => $tax_groups,
                -external_ids   => $external_ids,
                -names          => $names,
                -classes        => $classes,
                -families       => $families,
                -min_ic         => $min_ic,
                -order_by       => $order_column
            );
 Function : Fetch a list of TFInfo objects.
 Returns  : An listref of OPOSSUM::TFInfo objects.
 Args	  : Optionally key value pairs indicating selection criteria and
            table column to order on. All values may be a single values or
            a listref of values except min_ic which is a single value.

=cut

sub fetch_list
{
    my ($self, %args) = @_;

    my $where = '';

    my $sources = $args{-sources};
    if ($sources) {
        $where .= " and " if $where;
        if (ref $sources eq 'ARRAY') {
            $where .= "source in ('";
            $where .= join "','", @$sources;
            $where .= "')";
        } else {
            $where .= " source = '$sources'";
        }
    }

    my $collections = $args{-collections};
    if ($collections) {
        $where .= " and " if $where;
        if (ref $collections eq 'ARRAY') {
            $where .= "collection in ('";
            $where .= join "','", @$collections;
            $where .= "')";
        } else {
            $where .= " collection = '$collections'";
        }
    }

    my $tax_groups = $args{-tax_groups};
    if ($tax_groups) {
        $where .= " and " if $where;
        if (ref $tax_groups eq 'ARRAY') {
            $where .= "tax_group in ('";
            $where .= join "','", @$tax_groups;
            $where .= "')";
        } else {
            $where .= " tax_group = '$tax_groups'";
        }
    }

    my $names = $args{-names};
    if ($names) {
        $where .= " and " if $where;
        if (ref $names eq 'ARRAY') {
            $where .= "name in ('";
            $where .= join "','", @$names;
            $where .= "')";
        } else {
            $where .= " name = '$names'";
        }
    }

    my $tf_ids = $args{-tf_ids};
    if ($tf_ids) {
        $where .= " and " if $where;
        if (ref $tf_ids eq 'ARRAY') {
            $where .= "tf_id in (";
            $where .= join ",", @$tf_ids;
            $where .= ")";
        } else {
            $where .= " tf_id = $tf_ids";
        }
    }

    my $external_ids = $args{-external_ids};
    if ($external_ids) {
        $where .= " and " if $where;
        if (ref $external_ids eq 'ARRAY') {
            $where .= "external_id in (";
            $where .= join ",", @$external_ids;
            $where .= ")";
        } else {
            $where .= " external_id = $external_ids";
        }
    }

    my $min_ic = $args{-min_ic};
    if (defined $min_ic) {
        $where .= " and " if $where;
        $where .= "ic >= $min_ic";
    }

    my $order_cols = $args{-order_by} || undef;

    return $self->fetch_list_where($where, $order_cols);
}

=head2 fetch_list_where

 Title    : fetch_list_where
 Usage    : $tfi_list = $tfia->fetch_list_where($where);
 Function : Fetch TF info object(s)
 Returns  : An OPOSSUM::TFInfo object or listref of OPOSSUM::TFInfo objects.
 Args	  : Optional where clause,
            Optional column(s) to order on. May be single value or listref
            of values.

=cut

sub fetch_list_where
{
    my ($self, $where, $order_cols) = @_;

    if ($where && $where !~ /^\s*where\s+/) {
        $where = "where $where";
    }

    my $sql = qq{
        select tf_id, name, external_id, source, collection, class, family,
               tax_group, width, ic
        from tf_info
    };

    $sql .= " $where" if $where;

    my $order_clause = '';
    if ($order_cols) {
        if (ref $order_cols eq 'ARRAY') {
            $order_clause .= join ", ", @$order_cols;
        } else {
            $order_clause .= " $order_cols";
        }
    }

    $sql .= " order by $order_clause" if $order_clause;

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch TF info\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch TF info\n$sql\n" . $self->errstr;
        return;
    }

    my $row_count = 0;
    my @tf_info;
    while (my @row = $sth->fetchrow_array) {
        $row_count++;
        push @tf_info, OPOSSUM::TFInfo->new(
            -id          => $row[0],
            -name        => $row[1],
            -external_id => $row[2],
            -source      => $row[3],
            -collection  => $row[4],
            -class       => $row[5],
            -family      => $row[6],
            -tax_group   => $row[7],
            -width       => $row[8],
            -ic          => $row[9]
        );
    }
    $sth->finish;

    return @tf_info if wantarray();

    if ($row_count) {
        if ($row_count > 1) {
            return \@tf_info;
        }
        return $tf_info[0];
    }

    return undef;
}

=head2 fetch_set

 Title    : fetch_set
 Usage    : $tf_info_set = $tfia->fetch_set(
                -sources        => $sources,
                -collections    => $collections,
                -tax_groups     => $tax_groups,
                -external_ids   => $external_ids,
                -names          => $names,
                -classes        => $classes,
                -families       => $families,
                -min_ic         => $min_ic
            );
 Function : Fetch a set of TFInfo objects.
 Returns  : An OPOSSUM::TFInfoSet object.
 Args	  : Optionally key value pairs indicating selection criteria.
            All values may be a single values or a listref of values
            except min_ic which is a single value.

=cut

sub fetch_set
{
    my ($self, %args) = @_;

    my $tfi_list = $self->fetch_list(%args);

    return undef if !$tfi_list;

    my $tfi_set = OPOSSUM::TFInfoSet->new();
    if (ref $tfi_list eq 'ARRAY') {
        foreach my $tfi (@$tfi_list) {
            $tfi_set->add_tf_info($tfi);
        }
    } else {
        $tfi_set->add_tf_info($tfi_list);
    }

    foreach my $arg (keys %args) {
        $arg =~ s/^-//;
        if ($arg eq 'sources' || $arg eq 'collections' || $arg eq 'tax_groups'
            || $arg eq 'min_ic')
        {
            $tfi_set->param($arg, $args{"-$arg"});
        }
    }

    return $tfi_set;
}

=head2 fetch_set_where

 Title    : fetch_set_where
 Usage    : $tfi = $tfia->fetch_set_where($where);
 Function : Fetch a set of TFInfo objects.
 Returns  : An OPOSSUM::TFInfoSet object.
 Args	  : An optional where clause.

=cut

sub fetch_set_where
{
    my ($self, $where) = @_;

    my $tfi_list = $self->fetch_list_where($where);

    return undef if !$tfi_list;

    my $tfi_set = OPOSSUM::TFInfoSet->new();
    if (ref $tfi_list eq 'ARRAY') {
        foreach my $tfi (@$tfi_list) {
            $tfi_set->add_tf_info($tfi);
        }
    } else {
        $tfi_set->add_tf_info($tfi_list);
    }

    return $tfi_set;
}

=head2 fetch_tf_ids

 Title    : fetch_tf_ids
 Usage    : $ids = $tfia->fetch_tf_ids($where_clause);
 Function : Fetch a list of all the TF IDs, optionally using
 	    a where clause.
 Returns  : A list ref of integer TF IDs.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_tf_ids
{
    my ($self, $where_clause) = @_;

    my $sql = qq{select tf_id from tf_info};
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    $sql .= " order by tf_id";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TF IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TF IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($ids) = $sth->fetchrow_array) {
        push @ids, $ids;
    }

    return @ids ? \@ids : undef;
}

=head2 fetch_external_ids

 Title    : fetch_external_ids
 Usage    : $ids = $tfia->fetch_external_ids($where_clause);
 Function : Fetch a list of all the TFBS external IDs, optionally using
            a where clause.
 Returns  : A list ref of TFBS external ID strings.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_external_ids
{
    my ($self, $where_clause) = @_;

    my $sql = qq{select external_id from tf_info};
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    $sql .= " order by external_id";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TFBS external IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TFBS external IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($ids) = $sth->fetchrow_array) {
        push @ids, $ids;
    }

    return @ids ? \@ids : undef;
}

=head2 fetch_names

 Title    : fetch_names
 Usage    : $names = $tfia->fetch_names($where_clause);
 Function : Fetch a list of all the TFBS profile names, optionally using
            a where clause.
 Returns  : A list ref of TFBS profile name strings.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_names
{
    my ($self, $where_clause) = @_;

    my $sql = qq{select name from tf_info};
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    $sql .= " order by name";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TFBS names\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TFBS names\n" . $self->errstr;
        return;
    }

    my @names;
    while (my ($name) = $sth->fetchrow_array) {
        push @names, $name;
    }

    return @names ? \@names : undef;
}

=head2 fetch_sources

 Title    : fetch_sources
 Usage    : $sources = $tfia->fetch_sources($where_clause);
 Function : Fetch a list of all the TFBS profile sources, optionally using
            a where clause.
 Returns  : A list ref of TFBS profile sources.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_sources
{
    my ($self, $where_clause) = @_;

    my $sql = qq{
        select distinct source from tf_info where (source is not NULL
        and source != '')
    };

    if ($where_clause) {
        $sql .= " and ($where_clause)";
    }

    $sql .= " order by source";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TFBS sources\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TFBS sources\n" . $self->errstr;
        return;
    }

    my @sources;
    while (my ($source) = $sth->fetchrow_array) {
        push @sources, $source;
    }

    return @sources ? \@sources : undef;
}

=head2 fetch_collections

 Title    : fetch_collections
 Usage    : $collections = $tfia->fetch_collections($where_clause);
 Function : Fetch a list of all the TFBS profile collections, optionally using
            a where clause.
 Returns  : A list ref of TFBS profile collections.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_collections
{
    my ($self, $where_clause) = @_;

    my $sql = qq{
        select distinct collection from tf_info where (collection is not NULL
        and collection != '')
    };

    if ($where_clause) {
        $sql .= " and ($where_clause)";
    }

    $sql .= " order by collection";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TFBS collections\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TFBS collections\n" . $self->errstr;
        return;
    }

    my @collections;
    while (my ($collection) = $sth->fetchrow_array) {
        push @collections, $collection;
    }

    return @collections ? \@collections : undef;
}

=head2 fetch_tax_groups

 Title    : fetch_tax_groups
 Usage    : $tax_groups = $tfia->fetch_tax_groups($where_clause);
 Function : Fetch a list of all the TFBS profile tax_groups, optionally using
            a where clause.
 Returns  : A list ref of TFBS profile tax_groups.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_tax_groups
{
    my ($self, $where_clause) = @_;

    my $sql = qq{
        select distinct tax_group from tf_info where (tax_group is not NULL
        and tax_group != '')
    };

    if ($where_clause) {
        $sql .= " and ($where_clause)";
    }

    $sql .= " order by tax_group";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TFBS tax_groups\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TFBS tax_groups\n" . $self->errstr;
        return;
    }

    my @tax_groups;
    while (my ($tax_group) = $sth->fetchrow_array) {
        push @tax_groups, $tax_group;
    }

    return @tax_groups ? \@tax_groups : undef;
}

=head2 fetch_by_id

 Title    : fetch_by_id
 Usage    : $tfi = $tfia->fetch_by_id($id);
 Function : Fetch a TF info object by it's ID.
 Returns  : An OPOSSUM::TFInfo object.
 Args	  : An integer TF ID.

=cut

sub fetch_by_id
{
    my ($self, $id) = @_;

    my $where = "where tf_id = $id";

    return $self->fetch_where($where);
}

=head2 fetch_by_tf_id

 Title    : fetch_by_tf_id
 Usage    : $tfi = $tfia->fetch_by_tf_id($id);
 Function : Alternate name for fetch_by_id. Fetch a TF info object by it's
            ID.
 Returns  : An OPOSSUM::TFInfo object.
 Args	  : An integer TF ID.

=cut

sub fetch_by_tf_id
{
    my ($self, $id) = @_;

    return $self->fetch_by_id($id);
}

=head2 fetch_by_name

 Title    : fetch_by_name
 Usage    : $tfi = $tfia->fetch_by_name($name);
 Function : Fetch a TF info object by it's name.
 Returns  : An OPOSSUM::TFInfo object.
 Args	  : Optionally a TFBS name string.

=cut

sub fetch_by_name
{
    my ($self, $name) = @_;

    my $where = "where name = '$name'";

    return $self->fetch_where($where);
}

=head2 fetch_by_external_id

 Title    : fetch_by_external_id
 Usage    : $tfi = $tfia->fetch_by_external_id($xid);
 Function : Fetch a TF info object by it's external ID.
 Returns  : An OPOSSUM::TFInfo object.
 Args	  : Optionally a TFBS external ID string.

=cut

sub fetch_by_external_id
{
    my ($self, $external_id) = @_;

    my $where = "where external_id = '$external_id'";

    return $self->fetch_where($where);
}

=head2 store

 Title   : store
 Usage   : $id = $tfia->store($tf_info);
 Function: Store tf_info record in the database.
 Args    : The OPOSSUM::TFInfo object to store.
 Returns : A database ID of the newly stored tf_info record.

=cut

sub store
{
    my ($self, $tf_info) = @_;

    if (!ref $tf_info || !$tf_info->isa('OPOSSUM::TFInfo')) {
        carp "Not an OPOSSUM::TFInfo object";
        return;
    }

    my $sql = qq{
        insert into tf_info
            (source, collection, external_id, name, class, family, tax_group,
             width, ic)
        values (?,?,?,?,?,?,?,?,?)
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert tf_info statement\n" . $self->errstr;
        return;
    }

    if (
        !$sth->execute(
            $tf_info->source,      $tf_info->collection,
            $tf_info->external_id, $tf_info->name,
            $tf_info->class,       $tf_info->family,
            $tf_info->tax_group,   $tf_info->width,
            $tf_info->ic
        )
        )
    {
        carp "Error inserting tf_info record\n" . $self->errstr;
        return;
    }
    $sth->finish;

    return $sth->{'mysql_insertid'};
}

1;
