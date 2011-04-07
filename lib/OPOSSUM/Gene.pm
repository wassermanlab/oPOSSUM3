=head1 NAME

OPOSSUM::Gene - Gene object (genes DB record)

=head1 DESCRIPTION

A Gene object models a record retrieved from the genes table of the
oPOSSUM DB. Note that the start/end also should encompass the upstream
part of the gene sequence.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Gene;

use strict;

use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);

=head2 new

 Title   : new
 Usage   : $gene = OPOSSUM::Gene->new(
			    -id			    => '123',
			    -ensembl_id	    => 'ENSG00000165029',
			    -symbol 		=> 'ABCA1',
			    -biotype 		=> 'protein_coding',
			    -chr		    => '9',
			    -start		    => 106583104,
			    -end		    => 106730339,
			    -tss		    => 106720339,
			    -strand 		=> -1
            );

 Function: Construct a new Gene object
 Returns : a new OPOSSUM::Gene object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 id

 Title   : id
 Usage   : $id = $gene->id() or $gene->id('123');

 Function: Get/set the ID of the Gene. This should be a unique
           identifier for this object within the implementation. If
           the Gene object was read from the oPOSSUM database,
           this should be set to the value in the gene_id column.
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

#
# Synonym for id method
#
sub gene_id
{
    my ($self, $id) = @_;

    return $self->id($id);
}

=head2 ensembl_id

 Title   : ensembl_id
 Usage   : $ensembl_id = $gene->ensembl_id()
           or $gene->ensembl_id('ENSG00000165029');

 Function: Get/set the sequence ID of the Gene (same as Ensembl ID).
           This should be a unique identifier for this object within
           the implementation. 
 Returns : A Gene sequence ID
 Args    : None or a Gene sequence ID

=cut

sub ensembl_id
{
	my ($self, $id) = @_;

	if ($id) {
		$self->{-ensembl_id} = $id;
	}

	return $self->{-ensembl_id};
}

=head2 symbol

 Title   : symbol
 Usage   : $label = $gene->symbol() or $gene->symbol('ABCA1');

 Function: Get/set the symbol of the Gene. This should
           be a generally accepted name for this gene, i.e. a HUGO.
 Returns : A string
 Args    : None or a string

=cut

sub symbol
{
    my ($self, $symbol) = @_;

    if ($symbol) {
        $self->{-symbol} = $symbol;
    }

    return $self->{-symbol};
}

=head2 description

 Title   : description
 Usage   : $desc = $gene->description() or $gene->description($desc);

 Function: Get/set the gene description of the Gene.
 Returns : A string
 Args    : None or a string

=cut

sub description
{
    my ($self, $description) = @_;

    if ($description) {
        $self->{-description} = $description;
    }

    return $self->{-description};
}

=head2 biotype

 Title   : biotype
 Usage   : $biotype = $gene->biotype() or $gene->biotype('protein_coding');

 Function: Get/set the biotype of the Gene. This should be the same biotype as
			given by Ensembl
 Returns : A string
 Args	: None or a string

=cut

sub biotype
{
	my ($self, $biotype) = @_;
	
	if ($biotype) {
		$self->{-biotype} = $biotype;
	}
	
	return $self->{-biotype}
}

=head2 chr

 Title   : chr
 Usage   : $chr = $gene->chr() or $gene->chr($chr);

 Function: Get/set the chromosome name of the gene.
 Returns : A string
 Args    : None or a chromosome name

=cut

sub chr
{
    my ($self, $chr) = @_;

    if ($chr) {
        $self->{-chr} = $chr;
    }

    return $self->{-chr};
}

=head2 strand

 Title   : strand
 Usage   : $strand = $gene->strand() or $gene->strand(1);

 Function: Get/set the strand of the PromoterPair.
 Returns : 1 or -1
 Args    : None or a new strand value

=cut

sub strand
{
    my ($self, $strand) = @_;

    if ($strand) {
        $self->{-strand} = $strand;
    }

    return $self->{-strand};
}

=head2 start

 Title   : start
 Usage   : $start = $gene->start() or $gene->start($start);

 Function: Get/set the gene start position
 Returns : An integer
 Args    : None or a new gene start site 

=cut

sub start
{
    my ($self, $start) = @_;

    if ($start) {
        $self->{-start} = $start;
    }

    return $self->{-start};
}

=head2 end

 Title   : end
 Usage   : $end = $gene->end() or $gene->end($end);

 Function: Get/set the gene end position
 Returns : An integer
 Args    : None or a new gene end site 

=cut

sub end
{
    my ($self, $end) = @_;

    if ($end) {
        $self->{-end} = $end;
    }

    return $self->{-end};
}

=head2 tss

 Title   : tss
 Usage   : $tss = $gene->tss() or $gene->tss($tss);

 Function: Get/set the gene tss position
 Returns : An integer
 Args    : None or a new gene tss

=cut

sub tss
{
    my ($self, $tss) = @_;

    if ($tss) {
        $self->{-tss} = $tss;
    }

    return $self->{-tss};
}

=head2 operon

 Title   : operon
 Usage   : $label = $gene->operon() or $gene->operon($operon);

 Function: Get/set the Operon object associated with this gene
 Returns : An Operon object. If this is not an operon gene, Operon object
           would contain empty values (maybe NA as operon name?)
 Args	 : None or a new Operon object

=cut

sub operon
{
	my ($self, $operon) = @_;

	if ($operon) {
		if ($operon->isa("OPOSSUM::Operon")) {
            $self->{-operon} = $operon;
        } else {
            carp "not an OPOSSUM::Operon";
            return undef;
        }
    } elsif (!$self->{-operon}) {
        $self->fetch_operon();
    }

	return $self->{-operon};
}

=head2 fetch_operon

 Title   : fetch_operon
 Usage   : $member_genes = $gene->fetch_other_genes_in_operon()

 Function: Fetch the Operon object associated with this gene
 Returns : An Operon object
 Args    : None

=cut

sub fetch_operon
{
    my $self = shift;

    if (!$self->adaptor()) {
        carp "no adaptor defined trying to get Operon";
        return;
    }

    my $adaptor = $self->adaptor()->db()->get_OperonAdaptor();
    if (!$adaptor) {
        carp "could not get OperonAdaptor";
        return;
    }

    my $operon = $adaptor->fetch_by_gene_id($self->id());

    $self->{-operon} = $operon;
}

=head2 sequence

 Title   : sequence
 Usage   : $sequence = $gene->sequence()
           or $gene->sequence($sequence);

 Function: Get/set the Sequence object associated with this gene
 Returns : A Sequence object
 Args    : None or a new Sequence object 

=cut

sub sequence
{
    my ($self, $sequence) = @_;

    if ($sequence) {
        if ($sequence->isa("OPOSSUM::Sequence")) {
            $self->{-sequence} = $sequence;
        } else {
            carp "not an OPOSSUM::Sequence";
            return undef;
        }
    } elsif (!$self->{-sequence}) {
        $self->fetch_sequence();
    }

    return $self->{-sequence};
}

=head2 fetch_sequence

 Title   : fetch_sequence
 Usage   : $sequence = $gene->fetch_sequence()

 Function: Get the Sequence object associated with this gene.
 Returns : A Sequence object
 Args    : None

=cut

sub fetch_sequence
{
    my $self = shift;

    if (!$self->adaptor()) {
        carp "no adaptor defined trying to get Sequence";
        return;
    }

    my $adaptor = $self->adaptor()->db()->get_SequenceAdaptor();
    if (!$adaptor) {
        carp "could not get SequenceAdaptor";
        return;
    }

    my $sequence = $adaptor->fetch_by_gene_id($self->id());

    $self->{-sequence} = $sequence;
}

=head2 promoters

 Title   : promoters
 Usage   : $promoters = $gene->promoters()
           or $gene->promoters($promoters);

 Function: Get/set the Promoter objects associated with this gene
 Returns : A listref of Promoter objects
 Args    : None or a listref of new Promoter objects

=cut

sub promoters
{
    my ($self, $promoters) = @_;

    if ($promoters) {
        if (ref $promoters eq 'ARRAY'
            && $promoters->[0]->isa("OPOSSUM::Promoter"))
        {
            $self->{-promoters} = $promoters;
        } else {
            carp "not a listref of OPOSSUM::Promoter objects";
            return undef;
        }
    } elsif (!$self->{-promoters}) {
        $self->fetch_promoters();
    }

    return $self->{-promoters};
}

=head2 fetch_promoters

 Title   : fetch_promoters
 Usage   : $promoters = $gene->fetch_promoters()

 Function: Fetch the Promoter objects associated with this gene.
 Returns : A listref of Promoter objects
 Args    : None

=cut

sub fetch_promoters
{
    my $self = shift;

    if (!$self->adaptor()) {
        carp "no adaptor defined trying to get Promoter";
        return;
    }

    my $adaptor = $self->adaptor()->db()->get_PromoterAdaptor();
    if (!$adaptor) {
        carp "could not get PromoterAdaptor";
        return;
    }

    my $promoters = $adaptor->fetch_by_gene_id($self->id());

    $self->{-promoters} = $promoters;
}



=head2 conserved_regions

 Title   : conserved_regions
 Usage   : $conserved_regions = $gene->conserved_regions()
           or
           $gene->conserved_regions($conservation_level, $conserved_regions);

 Function: Get/set the list ConservedRegion objects associated with this
           gene
 Returns : A listref of ConservedRegion objects
 Args    : An integer conservation level
           A listref of ConservedRegion objects

=cut

sub conserved_regions
{
    my ($self, $level, $conserved_regions) = @_;

    $level = 1 if !$level;

    if ($conserved_regions) {
        if (ref $conserved_regions eq 'ARRAY'
            && $conserved_regions->[0]->isa("OPOSSUM::ConservedRegion"))
        {
            $self->{-conserved_regions} = $conserved_regions;
        } else {
            carp "not a listref of OPOSSUM::ConservedRegion objects";
            return undef;
        }
    } elsif (!$self->{-conserved_regions}) {
        $self->fetch_conserved_regions($level);
    }

    return $self->{-conserved_regions}{$level};
}

=head2 fetch_conserved_regions

 Title   : fetch_conserved_regions
 Usage   : $conserved_regions = $gene->fetch_conserved_regions()

 Function: Fetch the ConservedRegion objects associated with this gene.
 Returns : A listref of ConservedRegion objects
 Args    : None

=cut

sub fetch_conserved_regions
{
    my ($self, $cons_level) = @_;

    $cons_level = 1 if !$cons_level;

    if (!$self->adaptor()) {
        carp "no adaptor defined trying to get ConservedRegion";
        return;
    }

    my $adaptor = $self->adaptor()->db()->get_ConservedRegionAdaptor();
    if (!$adaptor) {
        carp "could not get ConservedRegionAdaptor";
        return;
    }

    my $conserved_regions = $adaptor->fetch_by_gene_id(
        $self->id(), $cons_level
    );

    $self->{-conserved_regions}{$cons_level} = $conserved_regions;
}

=head2 exons

 Title   : exons
 Usage   : $exons = $gene->exons()
           or $gene->exons($exons);

 Function: Get/set the Exon objects associated with this gene
 Returns : A listref of Exon objects
 Args    : None or a listref of new Exon objects

=cut

sub exons
{
    my ($self, $exons) = @_;

    if ($exons) {
        if (ref $exons eq 'ARRAY'
            && $exons->[0]->isa("OPOSSUM::Exon"))
        {
            $self->{-exons} = $exons;
        } else {
            carp "not a listref of OPOSSUM::Exon objects";
            return undef;
        }
    } elsif (!$self->{-exons}) {
        $self->fetch_exons();
    }

    return $self->{-exons};
}

=head2 fetch_exons

 Title   : fetch_exons
 Usage   : $exons = $gene->fetch_exons()

 Function: Fetch the Exon objects associated with this gene.
 Returns : A listref of Exon objects
 Args    : None

=cut

sub fetch_exons
{
    my $self = shift;

    if (!$self->adaptor()) {
        carp "no adaptor defined trying to get Exon";
        return;
    }

    my $adaptor = $self->adaptor()->db()->get_ExonAdaptor();
    if (!$adaptor) {
        carp "could not get ExonAdaptor";
        return;
    }

    my $exons = $adaptor->fetch_by_gene_id($self->id());

    $self->{-exons} = $exons;
}

=head2 promoter_search_regions

 Title   : promoter_search_regions
 Usage   : $regions = $gene->promoter_search_regions(
                $upstream,
                $downstream
           );
 Function: Define promoter search regions based on the amount of upstream
           and downstream sequence specified around the TSSs (promoters)
           of the gene. Truncate regions at 3' end of gene. Combine regions
           into unique non-overlapping regions.
           NOTE: Return search regions in sequence (1-based) coordinates.
 Returns : A listref of Bio::SeqFeature::Generic objects defining the
           promoter search regions.
 Args    : Optional amount of upstream sequence,
           Optional amount of downstream sequence

=cut

sub promoter_search_regions
{
    my ($self, $upstream, $downstream) = @_;

    my $start     = $self->start();
    my $end       = $self->end();
    my $strand    = $self->strand();
    my $promoters = $self->promoters();

    my @search_regions;
    foreach my $p (@$promoters) {
        my $tss = $p->tss;

        my $sr_start;
        my $sr_end;
        if ($strand == 1) {
            if (defined $upstream) {
                $sr_start = $tss - $upstream;
                $sr_start = $start if $start > $sr_start;
            } else {
                $sr_start = $start;
            }

            if (defined $downstream) {
                $sr_end = $tss + $downstream - 1;
                $sr_end = $end if $end < $sr_end;
            } else {
                $sr_end = $end;
            }
        } elsif ($strand == -1) {
            if (defined $upstream) {
                $sr_end = $tss + $upstream;
                $sr_end = $end if $end < $sr_end;
            } else {
                $sr_end = $end;
            }

            if (defined $downstream) {
                $sr_start = $tss - $downstream + 1;
                $sr_start = $start if $start > $sr_start;
            } else {
                $sr_start = $start;
            }
        } else {
            carp "no strand defined";
            return;
        }

        # Convert to sequence (1-based) coords
        $sr_start = $sr_start - $start + 1;
        $sr_end   = $sr_end - $start + 1;

        push @search_regions, Bio::SeqFeature::Generic->new(
            -start	=> $sr_start,
            -end	=> $sr_end
        );
    }

    if (@search_regions) {
        return _combine_search_regions(\@search_regions);
    } else {
        return undef;
    }
}

=head2 conserved_region_length

 Title   : conserved_region_length
 Usage   : $length = $gene->conserved_region_length(
                $consevation_level,
                $upstream,
                $downstream
           );
 Function: Compute the total length of all the conserved regions for this
           gene at the given level of conservation. If upstream or
           downstream bp are specified, include only the portions of the
           conserved regions which fall strictly within the promoter
           regions.
 Returns : An integer value representing the total length of the conserved
           regions at the specified level of conservation.
 Args    : Conservation level
           Optional amount of upstream sequence,
           Optional amount of downstream sequence

=cut

sub conserved_region_length
{
    my ($self, $cons_level, $upstream, $downstream) = @_;

    my $cons_regions = $self->conserved_regions($cons_level);
    return 0 if !$cons_regions;

    my $ps_regions;
    if (defined $upstream || defined $downstream) {
        $ps_regions = $self->promoter_search_regions($upstream, $downstream);
    }

    my $length = 0;
    if ($ps_regions) {
    } else {
        foreach my $cr (@$cons_regions) {
            $length += $cr->length();
        }
    }

    return $length;
}

#
# Combine any promoter search regions which may overlap 
#
sub _combine_search_regions
{
    my ($regs) = @_;

    return if !$regs || !$regs->[0];

    @$regs = sort {$a->start <=> $b->start} @$regs;

    my $num_regs = scalar @$regs;
    for (my $i = 0; $i < $num_regs; $i++) {
        my $reg1 = $regs->[$i] if exists($regs->[$i]);
        if ($reg1) {
            for (my $j = $i+1; $j < $num_regs; $j++) {
                my $reg2 = $regs->[$j] if exists ($regs->[$j]);
                if ($reg2) {
                    if (_do_features_combine($reg1, $reg2)) {
                        if ($reg2->start < $reg1->start) {
                            $reg1->start($reg2->start);
                        }

                        if ($reg2->end > $reg1->end) {
                            $reg1->end($reg2->end);
                        }
                        delete $regs->[$j];
                    } else {
                        last;
                    }
                }
            }
        }
    }

    my @unique_regs;
    foreach my $reg (@$regs) {
        if (defined $reg) {
            push @unique_regs, $reg;
        }
    }

    return @unique_regs ? \@unique_regs : undef;
}

#
# Check if features should be combined. They should if they overlap or are
# adjacent (no gap between them).
#
sub _do_features_combine
{
    my ($feat1, $feat2) = @_;

    my $combine = 1;
    $combine = 0 if $feat1->start > $feat2->end + 1
        || $feat1->end < $feat2->start - 1;

    return $combine;
}

1;
