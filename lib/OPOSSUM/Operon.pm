
=head1 NAME

OPOSSUM::Operon - Operon object (operon DB record)

=head1 DESCRIPTION

A Operon object models a record retrieved from the operons table of the
oPOSSUM DB. It contains the gene id's of those genes that are associated
with this operon.

=head1 AUTHOR

 Andrew Kwon
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Operon;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);


=head2 new

 Title   : new
 Usage   : $operon = OPOSSUM::Operon->new(
		    -id				    => 1,
		    -symbol		        => CEOP1906,
		    );

 Function: Construct a new Operon object
 Returns : a new OPOSSUM::Operon object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 id

 Title   : id
 Usage   : $id = $operon->id() or $operon->id('33992');
 
 Function: Get/set the ID of the Operon. This should be a unique
           identifier for this object within the implementation. If
           the Operon object was read from the oPOSSUM database,
           this should be set to the value in the operon_id column.
 Returns : A numeric ID
 Args    : None or a numeric ID

=cut

sub id
{
    my ($self, $id) = @_;

    if ($id) {
        $self->{-id} = $id;
    }

    return $self->{-id};
}

=head2 symbol

 Title   : symbol
 Usage   : $op_symbol = $operon->symbol() or $operon->symbol($op_symbol);

 Function: Get/set the symbol(name) of the Operon object associated with this
           Operon.
 Returns : An operon symbol
 Args    : None or an operon symbol

=cut

sub symbol
{
    my ($self, $symbol) = @_;

    if ($symbol) {
        $self->{-symbol} = $symbol;
    }

    return $self->{-symbol};
}


=head2 genes

 Title   : genes
 Usage   : $genes = $op->genes()
           or $op->genes($genes);

 Function: Get/set the Gene objects associated with this operon. 
 Returns : A listref of Gene objects. The order of genes should reflect the order
 in the operon, 5' to 3'
 Args    : None or a listref of new Gene objects

=cut

sub genes
{
    my ($self, $genes) = @_;

    if ($genes) {
        if (ref $genes eq 'ARRAY'
            && $genes->[0]->isa("OPOSSUM::Gene"))
        {
            $self->{-genes} = $self->sort_genes($genes);
        } else {
            carp "not a listref of OPOSSUM::Gene objects";
            return undef;
        }
    }

    return $self->{-genes};
}

=head2 contains_gene

 Title   : contains_gene
 Usage   : $gene_id = $op->contains_gene($gid);

 Function: Checks whether the Operon object contains the Gene object with the given Gene ID. 
 Returns : null or Gene object with the specified Gene ID
 Args    : Gene ID Integer

=cut

sub contains_gene
{
    my ($self, $gid) = @_;

	my $gene_list = $self->{-gene_id_list};
    if (!$gene_list) {
		return;
	}
	
	return $$gene_list{$gid};
}

=head2 add_gene_by_id

 Title   : add_gene_by_id
 Usage   : $gene = $operon->add_gene_by_id($gene_id);

 Function: Adds a gene object to the operon's gene list, specified by the gene id
 Returns : Gene object with the specified gene_id
 Args    : Gene ID integer

=cut

sub add_gene_by_id
{
	my ($self, $gid) = @_;
	
	if (!$gid) {
		carp "must provide a valid gene_id";
		return undef;
	}
	
	my $gid_list = $self->{-gene_id_list};
	if ($gid_list and $$gid_list{$gid}) {
		return $$gid_list{$gid};
	}
	
    my $adaptor = $self->adaptor()->db()->get_GeneAdaptor();
    if (!$adaptor) {
        carp "could not get GeneAdaptor";
        return;
    }

	my $gene = $adaptor->fetch_by_gene_id($gid);
	if (!$gene) {
		carp "not a valid gene_id: $gid";
		return;
	}
	
	if (!$self->{-genes}) {
		push @{$self->{-genes}}, $gene;
	} elsif (scalar @{$self->{-genes}} > 0) {
		push @{$self->{-genes}}, $gene;
		$self->{-genes} = $self->sort_genes($self->{-genes});
	} else {
		$self->{-genes} = ($gene);
	}
	
	$self->{-gene_id_list}->{$gid} = $gene;
	
	return $gene;
}

=head2 fetch_first_gene

 Title   : fetch_first_gene
 Usage   : $op_first_gene = $operon->fetch_first_gene();

 Function: Get the first Gene object associated with this Operon.
 Returns : A Gene object
 Args    : None

=cut

sub fetch_first_gene
{
    my ($self) = @_;

    if (!$self->{-genes}) {
		carp "No genes stored in this operon";
		return;
	}
	
    return $self->{-genes}->[0];
}

=head2 fetch_last_gene

 Title   : fetch_last_gene
 Usage   : $op_last_gene = $operon->fetch_last_gene();

 Function: Get the last Gene object associated with this Operon.
 Returns : A Gene object
 Args    : None

=cut

sub fetch_last_gene
{
    my ($self) = @_;

    if (!$self->{-genes}) {
		carp "No genes stored in this operon";
		return;
	}
	
	my $last_index = scalar(@{$self->{-genes}}) - 1;
    
	return $self->{-genes}->[$last_index];
}

=head2 chr

 Title   : chr
 Usage   : $chr = $operon->chr();

 Function: Get the chromosome name of the operon
 Returns : A string
 Args    : None

=cut

sub chr
{
    my $self = shift;

	if (!defined $self->fetch_first_gene()) {
		carp "No genes stored in this operon";
		return;
	}
	
    return $self->fetch_first_gene()->chr();
}

=head2 strand

 Title   : strand
 Usage   : $strand = $operon->strand();

 Function: Get the strand of this operon
 Returns : 1 or -1
 Args    : None

=cut

sub strand
{
    my $self = shift;
	if (!defined $self->fetch_first_gene()) {
		carp "No genes stored in this operon";
		return;
	}
    return $self->fetch_first_gene()->strand();
}

=head2 start

 Title   : start
 Usage   : $start = $operon->start();

 Function: Get the chromosomal start position of this operon
 Returns : An integer
 Args    : None

=cut

sub start
{
    my $self = shift;
	if (!defined $self->fetch_first_gene()) {
		carp "No genes stored in this operon";
		return;
	}

    # added code for -ve strand operons DJA 2010/11/30
    if ($self->strand() == 1) {
        return $self->fetch_first_gene()->start();
    } else {
        return $self->fetch_last_gene()->start();
    }
}

=head2 end

 Title   : end
 Usage   : $end = $operon->end();

 Function: Get the chromosomal end position of this operon
 Returns : An integer
 Args    : None

=cut

sub end
{
    my $self = shift;
	if (!defined $self->fetch_last_gene()) {
		carp "No genes stored in this operon";
		return;
	}

    # added code for -ve strand operons DJA 2010/11/30
    if ($self->strand() == 1) {
        return $self->fetch_last_gene()->end();
    } else {
        return $self->fetch_first_gene()->end();
    }
}

sub sort_genes
{
	my ($self, $genes) = @_;
	# make sure that these are in order
	my $strand = $$genes[0]->strand;
	my @sorted_genes;
	if ($strand == 1) {
		@sorted_genes = sort {$a->start <=> $b->start} @$genes;
	} else {
		@sorted_genes = sort {$b->start <=> $a->start} @$genes;
	}
	
	return \@sorted_genes;
}

1;
