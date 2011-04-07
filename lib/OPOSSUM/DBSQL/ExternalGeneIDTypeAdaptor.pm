
=head1 NAME

OPOSSUM::DBSQL::ExternalGeneIDTypeAdaptor - Adaptor for MySQL queries to
retrieve and store external gene ID types.

=head1 SYNOPSIS

$xgita = $db_adaptor->get_ExternalGeneIDTypeAdaptor();

=head1 DESCRIPTION

The external_gene_id_types table contains records which map the external gene
ID type numbers to their corresponding names.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ExternalGeneIDTypeAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::ExternalGeneIDType;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_id_types

 Title    : fetch_id_types
 Usage    : $id_types = $xgita->fetch_id_types();
 Function : Fetch a list of all the external ID type numbers from the
            database.
 Returns  : A list ref of integer ID type numbers.
 Args	  : None.

=cut

sub fetch_id_types
{
    my ($self) = @_;

    my $sql = qq{select id_type from external_gene_id_types order by id_type};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching ID types\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching ID types\n" . $self->errstr;
        return;
    }

    my @id_types;
    while (my ($id_type) = $sth->fetchrow_array) {
        push @id_types, $id_type;
    }

    return @id_types ? \@id_types : undef;
}

=head2 fetch_names

 Title    : fetch_names
 Usage    : $names = $xgita->fetch_names();
 Function : Fetch a list of all the external ID type names from the
            database.
 Returns  : A list ref of external ID type name strings.
 Args	  : None.

=cut

sub fetch_names
{
    my ($self) = @_;

    my $sql = qq{select name from external_gene_id_types order by name};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch external gene ID type names\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch external gene ID type names\n"
            . $self->errstr;
        return;
    }

    my @names;
    while (my ($name) = $sth->fetchrow_array) {
        push @names, $name;
    }

    return @names ? \@names : undef;
}

=head2 fetch_dblink_names

 Title    : fetch_dblink_names
 Usage    : $dblink_names = $xgita->fetch_dblink_names();
 Function : Fetch a list of all the external ID type dblink_names from the
            database.
 Returns  : A list ref of external ID type dblink_name strings.
 Args	  : None.

=cut

sub fetch_dblink_names
{
    my ($self) = @_;

    my $sql = qq{select dblink_name from external_gene_id_types
        order by dblink_name};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch external gene ID type dblink_names\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch external gene ID type dblink_names\n"
            . $self->errstr;
        return;
    }

    my @names;
    while (my ($name) = $sth->fetchrow_array) {
        push @names, $name;
    }

    return @names ? \@names : undef;
}

=head2 fetch_by_id_type

 Title    : fetch_by_id_type
 Usage    : $xgit = $xgita->fetch_by_id_type($id_type);
 Function : Fetch an external gene ID type record by it's external ID type
            number.
 Returns  : An ExternalGeneIDType object.
 Args	  : An external ID type number.

=cut

sub fetch_by_id_type
{
    my ($self, $id) = @_;

    return if !$id;

    my $sql = qq{select id_type, name, dblink_name from external_gene_id_types
    		where id_type = $id};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch external gene ID type\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch external gene ID type\n" . $self->errstr;
        return;
    }

    my $ext_id_type;
    if (my @row = $sth->fetchrow_array) {
        $ext_id_type = OPOSSUM::ExternalGeneIDType->new(
            -id_type        => $row[0],
            -name           => $row[1],
            -dblink_name    => $row[2]
        );
    }
    $sth->finish;

    return $ext_id_type;
}

=head2 fetch_by_name

 Title    : fetch_by_name
 Usage    : $xgit = $xgita->fetch_by_name($name);
 Function : Fetch an external gene ID type record by it's external ID type
 	        name.
 Returns  : An ExternalGeneIDType object.
 Args	  : An external ID type name string.

=cut

sub fetch_by_name
{
    my ($self, $name) = @_;

    return if !$name;

    my $sql = qq{select id_type, name, dblink_name from external_gene_id_types
        where name like '%$name%'};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch external gene ID type\n" . $self->errstr;
        return;
    }

    if (!$sth->execute()) {
        carp "error excuting fetch external gene ID type\n" . $self->errstr;
        return;
    }

    my $ext_id_type;
    if (my @row = $sth->fetchrow_array) {
        $ext_id_type = OPOSSUM::ExternalGeneIDType->new(
            -id_type        => $row[0],
            -name           => $row[1],
            -dblink_name    => $row[2]
        );
    }
    $sth->finish;

    return $ext_id_type;
}

=head2 fetch_external_gene_id_type_hash

 Title    : fetch_external_gene_id_type_hash
 Usage    : $xgit = $xgita->fetch_external_gene_id_hash();
 Function : Fetch all external gene ID type records and return as a
            hash keyed on the id_type.
 Returns  : A hash ref of ExternalGeneIDType objects.
 Args	  : Optional where clause.

=cut

sub fetch_external_gene_id_type_hash
{
    my ($self, $where_clause) = @_;

    my $xgit_list = $self->fetch_where($where_clause);

    my %xgit_hash;
    foreach my $xgit (@$xgit_list) {
        $xgit_hash{$xgit->id_type()} = $xgit;
    }

    return %xgit_hash ? \%xgit_hash : undef;
}

=head2 fetch_where

 Title    : fetch_where
 Usage    : $xgit = $xgita->fetch_where($where_clause);
 Function : Fetch one or more external gene ID type records with an
            optional where clause.
 Returns  : A list ref of ExternalGeneIDType objects.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_where
{
    my ($self, $where_clause) = @_;

    my $sql = qq{select id_type, name, dblink_name from external_gene_id_types};

    if ($where_clause) {
        unless ($where_clause =~ /^\s*where /) {
            $sql .= " where";
        }
        $sql .= " $where_clause";
    }

    $sql .= " order by id_type";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch external gene ID types\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch external gene ID types\n" . $self->errstr;
        return;
    }

    my @ext_id_types;
    while (my @row = $sth->fetchrow_array()) {
        push @ext_id_types, OPOSSUM::ExternalGeneIDType->new(
            -id_type        => $row[0],
            -name           => $row[1],
            -dblink_name    => $row[2]
        );
    }
    $sth->finish;

    return @ext_id_types ? \@ext_id_types : undef;
}

1;
