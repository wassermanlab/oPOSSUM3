
=head1 NAME

OPOSSUM::DBSQL::PromoterAdaptor - Adaptor for MySQL queries to retrieve
and store Promoter objects.

=head1 SYNOPSIS

$pa = $db_adaptor->get_PromoterAdaptor();

=head1 DESCRIPTION

The promoters table contains records which store promoter information
for the associated gene. Currently there is a one-to-one mapping of
genes to promoters although in theory there could be multiple
promoters for a given gene. This table stores the various chromosomal
coordinate information which was used to retrieve sequences, align them
and scan for TFBSs.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::PromoterAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Promoter;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_gene_ids

 Title    : fetch_gene_ids
 Usage    : $ids = $pa->fetch_gene_ids($where_clause);
 Function : Fetch a list of all the gene IDs in the database which
            have associated promoters, optionally using a where clause.
 Returns  : A list ref of integer gene IDs.
 Args	  : Optionally an SQL where clause.

=cut

sub fetch_gene_ids
{
    my ($self, $where_clause) = @_;

    my $sql = "select distinct gene_id from promoters";
    if ($where_clause) {
        $sql .= " where $where_clause";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching gene IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching gene IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_gene_id_fields

 Title    : fetch_gene_id_fields
 Usage    : $ids = $pa->fetch_gene_id_fields($id_field);
 Function : Fetch a list of all $id_field column values of the genes
            which have associated promoters.
 Returns  : A list ref of gene $id_field column values.
 Args	  : None.

=cut

sub fetch_gene_id_fields
{
    my ($self, $id_field) = @_;

    if (!$id_field) {
        return $self->fetch_gene_ids;
    }

    my $sql = qq{select distinct g.$id_field
                from genes g, promoters p
                where g.gene_id = p.gene_id};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching gene $id_field fields\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching gene $id_field fields\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_by_gene_id

 Title    : fetch_by_gene_id
 Usage    : $pps = $pa->fetch_by_gene_id($id);
 Function : Fetch promoters by their associated gene ID.
 Returns  : A list ref of OPOSSUM::Promoter objects.
 Args	  : A gene ID.

=cut

sub fetch_by_gene_id
{
    my ($self, $gene_id) = @_;

    my $sql = qq{select gene_id, tss, ensembl_transcript_id from promoters
                 where gene_id = $gene_id order by tss};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching promoters with gene_id ="
            . " $gene_id\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching promoters with gene_id ="
            . " $gene_id\n"
            . $self->errstr;
        return;
    }

    my @promoters;
    while (my @row = $sth->fetchrow_array) {
        push @promoters,
            OPOSSUM::Promoter->new(
            -adaptor               => $self,
            -gene_id               => $row[0],
            -tss                   => $row[1],
            -ensembl_transcript_id => $row[2]
            );
    }

    return @promoters ? \@promoters : undef;
}

=head2 store

 Title   : store
 Usage   : $id = $pa->store($promoter);
 Function: Store promoter in the database.
 Args    : The promoter (OPOSSUM::Promoter) to store.
 Returns : A database ID of the newly stored promoter.

=cut

sub store
{
    my ($self, $promoter) = @_;

    if (!ref $promoter || !$promoter->isa('OPOSSUM::Promoter')) {
        carp "Not an OPOSSUM::Promoter object";
        return;
    }

    my $sql = qq{insert into promoters (gene_id, tss, ensembl_transcript_id)
           values (?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert promoter statement\n" . $self->errstr;
        return;
    }

    if (
        !$sth->execute(
            $promoter->gene_id, $promoter->tss,
            $promoter->ensembl_transcript_id
        )
        )
    {
        carp "Error inserting promoter\n" . $self->errstr;
        return;
    }

    return $sth->{'mysql_insertid'};
}

1;
