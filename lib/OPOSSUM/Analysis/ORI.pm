=head1 NAME

OPOSSUM::ORI.pm - module to compute the over-representation index (ORI)
for over-represented binding sites

=head1 SYNOPSIS

 $ori = OPOSSUM::Analysis::ORI->new();

 $results = $ori->calculate_ori(
				$background_counts,
				$target_counts);

=head1 DESCRIPTION

Given TFBS position weight matrix counts for each gene in a target and
background set, calculate the over-representation index that each 
binding site is overrepresented in the set of target genes.

The ORI is defined as

 ORI = (Dc/Db) * (Pc/Pb)

where

 Dc = density of binding sites in the promoters of co-expressed genes
 Db = density of binding sites in the promoters of background genes
 Pc = proportion of co-expressed genes containing the binding site
 Pb = proportion of background genes containing the binding site

Density is defined as

 Di = total number of binding sites / total promoter length analysed

Proportion is defined as

 Pi = number of genes containing the binding site / total number of genes

=head1 AUTHORS

 Original coding: Shannan Ho Sui

 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: shosui@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::ORI;

use strict;

use Carp;

use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::CountsIO;
use OPOSSUM::Analysis::ORIResult;
use OPOSSUM::Analysis::ORIResultSet;


=head2 new

 Title    : new
 Usage    : $ori = OPOSSUM::Analysis::ORI->new();
 Function : Create a new OPOSSUM::Analysis::ORI object.
 Returns  : An OPOSSUM::Analysis::ORI object.
 Args     : None

=cut

sub new
{
    my ($class) = @_;

    my $self = bless {}, ref $class || $class;

    return $self;
}

=head2 calculate_ORI

 Title    : calculate_ORI
 Usage    : $results = $ori->calculate_ORI(
					   $background_counts,
					   $target_counts);
 Function : Perform the ORI analysis to determine over-representation
 	    of TFBS profiles and return the results.
 Returns  : An OPOSSUM::Analysis::ORIResultSet object.
 Args     : background_counts	- An OPOSSUM::Analysis::Counts object
 				  containing the TFBS counts in the
				  background set of genes
	    target_counts	- An OPOSSUM::Analysis::Counts object
 				  containing the TFBS counts in the
				  target set of genes

=cut

sub calculate_ORI
{
    my ($self, $bg_counts, $t_counts) = @_;

    return if !$bg_counts || !$t_counts;

    if (!$bg_counts->isa("OPOSSUM::Analysis::Counts")) {
    	carp "background counts is not an OPOSSUM::Analysis::Counts object";
	return;
    }

    if (!$t_counts->isa("OPOSSUM::Analysis::Counts")) {
    	carp "test counts is not an OPOSSUM::Analysis::Counts object";
	return;
    }

    my $t_cr_len = $t_counts->total_cr_length;
    if (!$t_cr_len) {
	carp "error computing target set conserved region length";
	return;
    }

    my $bg_cr_len = $bg_counts->total_cr_length;
    if (!$bg_cr_len) {
	carp "error computing background set conserved region length";
	return;
    }

    my $results = OPOSSUM::Analysis::ORIResultSet->new();
    if (!$results) {
	carp "error creating ORI result set";
	return;
    }
    my $t_tfbs_ids = $t_counts->get_all_tfbs_ids;
    foreach my $tfbs_id (@$t_tfbs_ids) {
    	if (!$bg_counts->tfbs_exists($tfbs_id)) {
	    carp "TFBS ID $tfbs_id does not exist in background set"
		    . " - excluding from analysis";
	    next;
	}

	# target and background TFBS hits
    	my $Ht = $t_counts->tfbs_count($tfbs_id);
    	my $Hb = $bg_counts->tfbs_count($tfbs_id);

	# target and background TFBS densities
	my $Dt = $Ht / $t_cr_len;
	my $Db = $Hb / $bg_cr_len;

	# target and background gene counts
	my $t_gene_hits = $t_counts->tfbs_gene_count($tfbs_id);
	my $bg_gene_hits = $bg_counts->tfbs_gene_count($tfbs_id);
	my $total_t_genes = $t_counts->num_genes;
	my $total_bg_genes = $bg_counts->num_genes;

	# target and background gene proportions
	my $Pt = $t_gene_hits / $total_t_genes;
	my $Pb = $bg_gene_hits / $total_bg_genes;
	
	my $ORI;

	if($Db != 0 && $Pb != 0) {
	    $ORI = ($Dt/$Db) * ($Pt/$Pb);
	}
	else {
	    $ORI = 1e999;    # set to infinite
	}

	if (defined $ORI) {
	    $results->add_result(
			OPOSSUM::Analysis::ORIResult->new(
					-id		=> $tfbs_id,
					-t_rate		=> $Dt,
					-bg_rate	=> $Db,
					-t_gene_prop	=> $Pt,
					-bg_gene_prop	=> $Pb,
					-ori      	=> $ORI));
	}
    }
    $results->param('bg_cr_length', $bg_cr_len);
    $results->param('t_cr_length', $t_cr_len);

    return $results;
}

1;
