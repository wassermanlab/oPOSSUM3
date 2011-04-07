
=head1 NAME

OPOSSUM::DBSQL::ExonAdaptor - Adaptor for MySQL queries to retrieve and
store exon information.

=head1 SYNOPSIS

$ea = $db_adaptor->get_ExonAdaptor();

=head1 DESCRIPTION

The exons table of the oPOSSUM database stores the exons associated with
the gene records.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ExonAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Exon;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title   : new
 Usage   : $ea = OPOSSUM::DBSQL::ExonAdaptor->new($db_adaptor);
 Function: Construct a new ExonAdaptor object
 Args    : An OPOSSUM::DBSQL::DBAdaptor object
 Returns : a new OPOSSUM::DBSQL::ExonAdaptor object

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_by_gene_id

 Title   : fetch_by_gene_id
 Usage   : $exons = $gpa->fetch_by_gene_id($gene_id);
 Function: Fetch all the exons from the DB by gene_id. Exons are ordered
           by start coordinate.
 Args    : The gene ID.
 Returns : A listref of OPOSSUM::Exon objects.

=cut

sub fetch_by_gene_id
{
    my ($self, $gid) = @_;

    my $sql =
        qq{select gene_id, start, end from exons where gene_id = $gid
           order by start};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching exons for gene_id $gid\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching exons for gene_id $gid\n" . $self->errstr;
        return;
    }

    my @exons;
    while (my @row = $sth->fetchrow_array) {
        push @exons, OPOSSUM::Exon->new(
            -adaptor    => $self,
            -gene_id    => $row[0],
            -start      => $row[1],
            -end        => $row[2]
        );
    }

    return @exons ? \@exons : undef;
}

=head2 fetch_by_region

 Title   : fetch_by_region
 Usage   : $exons = $gpa->fetch_by_region($chr, $start, $end);
 Function: Fetch all the exons from the DB in the region defined by
           chromosom, start and end. Exons are ordered by start coordinate.
 Args    : Chromosome name,
           Start coordinate,
           End coordinate
 Returns : A listref of OPOSSUM::Exon objects.

=cut

sub fetch_by_region
{
    my ($self, $chr, $start, $end) = @_;

    my $sql =
        qq{select x.gene_id, x.start, x.end from exons x, genes g
           where x.gene_id = g.gene_id and g.chr = '$chr' and
           x.start <= $end and x.end >= $start order by x.start};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching exons by start/end $start/$end\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching exons by start/end $start/$end\n" . $self->errstr;
        return;
    }

    my @exons;
    while (my @row = $sth->fetchrow_array) {
        push @exons, OPOSSUM::Exon->new(
            -adaptor    => $self,
            -gene_id    => $row[0],
            -start      => $row[1],
            -end        => $row[2]
        );
    }

    return @exons ? \@exons : undef;
}

=head2 store

 Title   : store
 Usage   : $ea->store($exon);
 Function: Store exon in the database.
 Args    : An OPOSSUM::Exon object
 Returns : True on success, false otherwise.

=cut

sub store
{
    my ($self, $exon) = @_;

    my $sql = qq{insert into exons (gene_id, start, end) values (?, ?, ?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert exon statement\n" . $self->errstr;
        return;
    }

    if (!$sth->execute($exon->gene_id, $exon->start, $exon->end)) {
        carp sprintf(
            "Error inserting exon for gene_id %d" . $self->errstr,
            $exon->gene_id
        );
        return 0;
    }

    return 1;
}

=head2 store_list

 Title   : store_list
 Usage   : $ea->store_list($exons, $gid);
 Function: Store exons in the database.
 Args    : A listref of OPOSSUM::Exon objects,
           A Gene ID,
 Returns : True on success, false otherwise.

=cut

sub store_list
{
    my ($self, $exons, $gid) = @_;

    my $sql = qq{insert into exons (gene_id, start, end) values (?, ?, ?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert exon statement\n" . $self->errstr;
        return;
    }

    my $ok = 1;
    foreach my $exon (@$exons) {
        if (!$sth->execute($gid, $exon->start, $exon->end)) {
            carp sprintf("Error inserting exon for gene_id %s\n%s",
                $gid, $self->errstr);
            # keep trying to store exons but return error status...
            $ok = 0;
        }
    }

    return $ok;
}

1;
