=head1 NAME

OPOSSUM::DBSQL::SequenceAdaptor - Adaptor for MySQL queries to retrieve and
store Sequence objects.

=head1 SYNOPSIS

$aa = $db_adaptor->get_SequenceAdaptor();

=head1 DESCRIPTION

The sequences table of the oPOSSUM database stores the raw and masked 
sequences.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::SequenceAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Sequence;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title   : new
 Usage   : $aa = OPOSSUM::DBSQL::SequenceAdaptor->new($db_adaptor);
 Function: Construct a new SequenceAdaptor object
 Args    : An OPOSSUM::DBSQL::DBAdaptor object
 Returns : a new OPOSSUM::DBSQL::SequenceAdaptor object

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
 Usage   : $sequence = $aa->fetch_by_gene_id($gid);
 Function: Fetch an Sequence object from the DB using it's gene ID.
 Args    : The Gene ID.
 Returns : An OPOSSUM::Sequence object.

=cut

sub fetch_by_gene_id
{
    my ($self, $gid) = @_;

    my $sql = qq{
        select gene_id, seq, masked_seq from sequences
		where gene_id = $gid
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error fetching sequence with gene_id = $gid\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error fetching sequence with gene_id = $gid\n"
            . $self->errstr;
        return;
    }

    my $sequence;
    if (my @row = $sth->fetchrow_array) {
        $sequence = OPOSSUM::Sequence->new(
            -adaptor    => $self,
            -gene_id	=> $row[0],
            -seq 		=> $row[1],
            -masked_seq => $row[2])
    }
    $sth->finish;

    return $sequence;
}

=head2 store

 Title   : store
 Usage   : $id = $aa->store($sequence);
 Function: Store sequence in the database.
 Args    : The sequence (OPOSSUM::Sequence) to store
 Returns : True on success, false otherwise.

=cut

sub store
{
    my ($self, $sequence) = @_;

    if (!ref $sequence || !$sequence->isa('OPOSSUM::Sequence')) {
    	carp "Not an OPOSSUM::Sequence object";
        return;
    }

    my $sql = qq{insert into sequences (gene_id, seq, masked_seq)
        values (?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
    	carp "Error preparing insert sequence statement - "
		. $self->errstr . "\n";
        return;
    }

    if (!$sth->execute(
        $sequence->gene_id,
        $sequence->seq,
        $sequence->masked_seq))
    {
    	carp "Error inserting sequence - " . $self->errstr . "\n";
        return 0;
    }

    return 1;
}

1;
