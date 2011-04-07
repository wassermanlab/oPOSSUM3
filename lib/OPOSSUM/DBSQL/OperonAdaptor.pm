=head1 NAME

OPOSSUM::DBSQL::OperonAdaptor - Adaptor for MySQL queries to retrieve and
store Operon objects.

=head1 SYNOPSIS

$ga = $db_adaptor->get_OperonAdaptor();

=head1 DESCRIPTION

The operons table of the oPOSSUM database stores operonral operon information.

=head1 AUTHOR

 Andrew Kwon
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::OperonAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Operon;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title   : new
 Usage   : $ga = OPOSSUM::DBSQL::OperonAdaptor->new($db_adaptor);
 Function: Construct a new OperonAdaptor object
 Args    : An OPOSSUM::DBSQL::DBAdaptor object
 Returns : a new OPOSSUM::DBSQL::OperonAdaptor object

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_where

 Title   : fetch_where
 Usage   : $operon = $opa->fetch_where($where);
 Function: Generic fetch method. Fetch operon object(s) from the DB
           with the given where clause.
 Args    : Optionally, a where clause.
 Returns : Either an OPOSSUM::Operon object or a reference to an array of
           OPOSSUM::Operon objects depending on whether query returns 1 or
           more rows.

=cut

sub fetch_where
{
    my ($self, $where) = @_;

    if ($where && $where !~ /^\s*where /) {
        $where = "where $where";
    }

    my $sql =
        qq{select operon_id,
            symbol,
            gene_id
        from operons};

    $sql .= " $where" if $where;

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch operons:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch operons:\n$sql\n" . $self->errstr;
        return;
    }

    my %operons;
	my @operons;
    while (my @row = $sth->fetchrow_array) {
		
		my $operon = $operons{$row[0]};
		if (!$operon) {
			$operon = OPOSSUM::Operon->new(
				-adaptor	=> $self,
				-id			=> $row[0],
				-symbol		=> $row[1]
			);
			push @operons, $operon;
			$operons{$row[0]} = $operon;
		}
        $operon->add_gene_by_id($row[2]);
    }
    $sth->finish;
	
    return @operons if wantarray();

    if (scalar @operons == 1) {
		return $operons[0];
    } else {
		return \@operons;
	}

    return undef;
}

=head2 fetch_operon_ids

 Title   : fetch_operon_ids
 Usage   : $gids = $ga->fetch_operon_ids($where);
 Function: Fetch list of operon IDs from the DB.
 Args    : Optionally a where clause.
 Returns : Reference to a list of internal operon IDs. If no where clause
           is provided, returns all operon IDs in the database.

=cut

sub fetch_operon_ids
{
    my ($self, $where_clause) = @_;

    my $sql = "select operon_id from operons";
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching operon IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching operon IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_gene_ids_to_operon_ids_hash

 Title   : fetch_gene_ids_to_operon_ids_hash
 Usage   : $op_gid_hash = $oa->fetch_gene_ids_to_operon_ids_hash($gids);
 Function: Fetch the hash of gene ids to operon ids.
 Args    : listref of opossum gene ids. If no argument is passed, all gene ids
           are matched to operons.
 Returns : Hash ref

=cut

sub fetch_gene_ids_to_operon_ids_hash
{
	my ($self, $gids) = @_;
	
	my $sql = "select gene_id, operon_id from operons";
	
	if ($gids and ref $gids eq 'ARRAY') {
		$sql .= " where gene_id in ('";
		$sql .= join("','", @$gids);
		$sql .= "')";
	}
	
	my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching operon IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching operon IDs\n" . $self->errstr;
        return;
    }

    my %gid_oid;
    while (my ($gid, $oid) = $sth->fetchrow_array) {
        $gid_oid{$gid} = $oid;
    }
    $sth->finish;

    return %gid_oid ? \%gid_oid : undef;
}

=head2 fetch_operon_id_fields

 Title   : fetch_operon_id_fields
 Usage   : $ids = $ga->fetch_operon_id_fields($field, $where);
 Function: Fetch list of IDs of type specified by $field from the DB.
 Args    : Optionally a field type and where clause.
 Returns : Reference to a list of IDs. If no field type is specified,
           returns internal operon IDs, otherwise the IDs stored in the
           given field. If no where clause is provided, returns all IDs in
           the database.

=cut

sub fetch_operon_id_fields
{
    my ($self, $id_field, $where_clause) = @_;

    if (!$id_field) {
        return $self->fetch_operon_ids;
    }

    my $sql = "select distinct $id_field from operons";
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching operon $id_field fields\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching operon $id_field fields\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_by_id

 Title   : fetch_by_id
 Usage   : $operon = $ga->fetch_by_id($id);
 Function: Fetch a operon object from the DB using it's ID.
 Args    : The unique internal operon ID.
 Returns : An OPOSSUM::Operon object.

=cut

sub fetch_by_id
{
    my ($self, $id) = @_;

    my $where = "where operon_id = $id";

    return $self->fetch_where($where);
}

#
# Synonym for fetch_by_id
#
sub fetch_by_operon_id
{
    my ($self, $id) = @_;

    return $self->fetch_by_id($id);
}

=head2 fetch_by_operon_id_list

 Title   : fetch_by_operon_id_list
 Usage   : $operons = $ga->fetch_by_operon_id_list($id_list);
 Function: Fetch a list of operon objects from the DB according to a list
           of operon IDs.
 Args    : A reference to a list of unique internal operon IDs.
 Returns : A reference to a list of OPOSSUM::Operon objects.

=cut

sub fetch_by_operon_id_list
{
    my ($self, $operon_ids) = @_;

    my $sql =
        qq{select operon_id,
			symbol,
			gene_id
		from operons
		where operon_id = ?};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching operons\n" . $self->errstr;
        return;
    }

    my @operons;
	my %operons;
    foreach my $operon_id (@$operon_ids) {
        if (!$sth->execute($operon_id)) {
            carp "error fetching operon for operon ID $operon_id\n"
                . $self->errstr;
            return;
        }

        if (my @row = $sth->fetchrow_array) {
			my $operon = $operons{$row[0]};
			if (!$operon) {
				my $operon = OPOSSUM::Operon->new(
				-adaptor		=> $self,
                -id             => $row[0],
                -symbol         => $row[1]
				);
				push @operons, $operon;
			}
			$operon->add_gene_by_id($row[2]);
        }
    }
    $sth->finish;

    return @operons ? \@operons : undef;
}

sub fetch_by_gene_id
{
    my ($self, $gene_id) = @_;

    return if !$gene_id;
	
    my $sql = "select operon_id from operons where gene_id = $gene_id";
	
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch operons:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch operons:\n$sql\n" . $self->errstr;
        return;
    }

    my ($operon_id) = $sth->fetchrow_array;

    $sth->finish;
	
	return if !$operon_id;
	
	return $self->fetch_by_id($operon_id);
}

sub fetch_by_symbol
{
    my ($self, $symbol) = @_;

    return if !$symbol;

    my $where = "where symbol = '$symbol'";

    return $self->fetch_where($where);
}

=head2 store

 Title   : store
 Usage   : $id = $ga->store($operon);
 Function: Store operon in the database.
 Args    : The operon (OPOSSUM::Operon) to store.
 Returns : A database ID of the newly stored operon.

=cut

sub store
{
    my ($self, $operon) = @_;

    if (!ref $operon || !$operon->isa('OPOSSUM::Operon')) {
        carp "Not an OPOSSUM::Operon object";
        return;
    }

    my $sql
        = qq{insert into operons
		    (symbol, gene_id)
		    values (?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert operon statement\n" . $self->errstr;
        return;
    }

    if (!$sth->execute($operon->id, $operon->symbol, $operon->gene_id))
    {
        carp "Error inserting operon\n" . $self->errstr;
        return;
    }
    $sth->finish;

    return $sth->{'mysql_insertid'};
}

1;
