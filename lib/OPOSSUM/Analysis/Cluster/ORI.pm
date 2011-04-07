=head1 NAME

OPOSSUM::Analysis::Cluster::ORI.pm - module to compute the over-representation
index (ORI) for over-represented binding sites

=head1 SYNOPSIS

 $ori = OPOSSUM::Analysis::Cluster::ORI->new();

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

package OPOSSUM::Analysis::Cluster::ORI;

use strict;

use Carp;

use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::CountsIO;
use OPOSSUM::Analysis::Cluster::ORIResult;
use OPOSSUM::Analysis::Cluster::ORIResultSet;


=head2 new

 Title    : new
 Usage    : $ori = OPOSSUM::Analysis::Cluster::ORI->new();
 Function : Create a new OPOSSUM::Analysis::Cluster::ORI object.
 Returns  : An OPOSSUM::Analysis::Cluster::ORI object.
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
 	    of TFBS profile clusters and return the results.
 Returns  : An OPOSSUM::Analysis::Cluster::ORIResultSet object.
 Args     : background_counts	- An OPOSSUM::Analysis::Cluster::Counts object
 				  containing the TFBS cluster counts in the
				  background set of genes
	    target_counts	- An OPOSSUM::Analysis::Cluster::Counts object
 				  containing the TFBS cluster counts in the
				  target set of genes

=cut

sub calculate_ORI
{
    my ($self, $bg_counts, $t_counts) = @_;

    return if !$bg_counts || !$t_counts;

    if (!$bg_counts->isa("OPOSSUM::Analysis::Cluster::Counts")) {
    	carp "background counts is not an OPOSSUM::Analysis::Cluster::Counts object";
	return;
    }

    if (!$t_counts->isa("OPOSSUM::Analysis::Cluster::Counts")) {
    	carp "test counts is not an OPOSSUM::Analysis::Cluster::Counts object";
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

    my $results = OPOSSUM::Analysis::Cluster::ORIResultSet->new();
    if (!$results) {
	carp "error creating ORI result set";
	return;
    }
    my $t_cluster_ids = $t_counts->get_all_cluster_ids;
    foreach my $cluster_id (@$t_cluster_ids) {
    	if (!$bg_counts->cluster_exists($cluster_id)) {
	    carp "TFCluster ID $cluster_id does not exist in background set"
		    . " - excluding from analysis";
	    next;
	}

	# target and background TFBS hits
    	my $Ht = $t_counts->cluster_count($cluster_id);
    	my $Hb = $bg_counts->cluster_count($cluster_id);

	# target and background TFBS densities
	my $Dt = $Ht / $t_cr_len;
	my $Db = $Hb / $bg_cr_len;

	# target and background gene counts
	my $t_gene_hits = $t_counts->cluster_gene_count($cluster_id);
	my $bg_gene_hits = $bg_counts->cluster_gene_count($cluster_id);
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
			OPOSSUM::Analysis::Cluster::ORIResult->new(
					-id		=> $cluster_id,
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
