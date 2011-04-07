=head1 NAME

OPOSSUM::Zscore.pm - module to compute z-scores of over-represented binding
sites

=head1 SYNOPSIS

 $zscore = OPOSSUM::Analysis::Zscore->new();

 $results = $zscore->calculate_zscore($background_counts, $target_counts);

=head1 DESCRIPTION

Given TFBS position weight matrix counts for each gene in a target and
background set, calculate the probability that each binding site is
overrepresented in the set of target genes as a z-score.

Z-score is computed such that if the standard deviation S, is given by:

 S = sqrt(N * Rb * (1 - Rb)) 

where

 N  = total number of nucleotides searched for binding sites in the target
      set of genes
 Rb = rate of binding site nucleotides in the background set and is equal
      to the number of observed binding site nucleotides in the background
      set divided by the total number of nucleotides scanned in the
      background set

then the z-score is given by:

 Z = (Nt - Ne - 0.5) / S

where,

 Nt = number of observed binding site nucleotides in the target set of
      genes
 Ne = number of expected binding site nucleotides in the target set of
      genes and is equal to the number of observed binding site
      nucleotides in the backgtound set of genes multiplied by the ratio
      of total number of target nucleotides to total number of background
      nucleotides scanned

=head1 AUTHORS

 Original coding: Chris Walsh

 Modified by: Shannan Ho Sui

 Final module: David Arenillas

 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Zscore;

use strict;

use Carp;

use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::ZscoreResult;
use OPOSSUM::Analysis::ZscoreResultSet;
use Statistics::Distributions;

=head2 new

 Title    : new
 Usage    : $zscore = OPOSSUM::Analysis::Zscore->new();
 Function : Create a new OPOSSUM::Analysis::Zscore object.
 Returns  : An OPOSSUM::Analysis::Zscore object.
 Args     : None

=cut

sub new
{
    my ($class) = @_;

    my $self = bless {}, ref $class || $class;

    return $self;
}

=head2 calculate_Zscore

 Title    : calculate_Zscore
 Usage    : $results = $fisher->calculate_Zscore(
                $bg_counts,
                $t_counts,
                $bg_seq_len,
                $t_seq_len,
                $tf_set
            );

 Function : Perform the z-score analysis to determine over-representation
            of TFBS profiles and return the results.
 Returns  : An OPOSSUM::Analysis::ZscoreResultSet object.
 Args     : bg_counts   - An OPOSSUM::Analysis::Counts object containing
                          the TFBS counts in the background set of genes
            t_counts    - An OPOSSUM::Analysis::Counts object containing
                          the TFBS counts in the target set of genes
            bg_seq_len  - An integer defining the total search space
                          (sequence nucleotides) used to search TFBSs in
                          the background set of genes
            t_seq_len   - An integer defining the total search space
                          (sequence nucleotides) used to search TFBSs in
                          the target set of genes
            tf_set      - An OPOSSUM::TFSet containing the TF information
                          used in the analysis (background and target)

=cut

sub calculate_Zscore
{
    my ($self, $bg_counts, $t_counts, $bg_seq_len, $t_seq_len, $tf_set) = @_;

    if (!$bg_counts) {
        carp "background counts not specified\n";
        return;
    }
    
    if (!$t_counts) {
        carp "target counts not specified\n";
        return;
    }

    if (!$bg_seq_len) {
        carp  "background sequence length (search space) not"
            . " specified or 0 length\n";
        return;
    }

    if (!$t_seq_len) {
        carp  "target sequence length (search space) not"
            . " specified or 0 length\n";
        return;
    }

    if (!$tf_set) {
        carp  "TF set not specified\n";
        return;
    }

    if (!$bg_counts->isa("OPOSSUM::Analysis::Counts")) {
        carp "background counts is not an OPOSSUM::Analysis::Counts object\n";
        return;
    }

    if (!$t_counts->isa("OPOSSUM::Analysis::Counts")) {
        carp "target counts is not an OPOSSUM::Analysis::Counts object\n";
        return;
    }

    if (!$tf_set->isa("OPOSSUM::TFSet")) {
        carp "TF set is not an OPOSSUM::TFSet object\n";
        return;
    }

    my $C;
    if ($bg_seq_len) {
        $C = $t_seq_len / $bg_seq_len;
    } else {
        carp "Ratio of target to background sequence lengths is undefined";
        return;
    }

    #
    # Actually $C should never be > 1 but all it means is that the target
    # set is larger than the background set which theoretically doesn't matter
    # for the z-score. Also if $C == 0 this is not an error either and could
    # legitimately happen in cases where a small upstream/downstream region
    # was selected. It just means the target rate would be 0.
    #
    #if (!$C || $C > 1) {
    #	carp "error computing conserved length ratio coefficient";
    #	return;
    #}

    my $results = OPOSSUM::Analysis::ZscoreResultSet->new();
    if (!$results) {
        carp "error creating Zscore result set";
        return;
    }
    my $t_tf_ids = $t_counts->tf_ids();
    foreach my $tf_id (@$t_tf_ids) {
        if (!$bg_counts->tf_exists($tf_id)) {
            carp "TF ID $tf_id does not exist in background set"
                . " - excluding from analysis";
            next;
        }

        # TFBS profile width
        my $width = $tf_set->get_matrix($tf_id)->length();

        # target and background TFBS hits
        my $Ht = $t_counts->tfbs_count($tf_id);
        my $Hb = $bg_counts->tfbs_count($tf_id);

        # target and background TFBS nucleotides
        my $Nt = $width * $Ht;
        my $Nb = $width * $Hb;

        # target and background TFBS nucleotide rates
        my $Rt = $Nt / $t_seq_len;
        my $Rb = $Nb / $bg_seq_len;

        # number of trials
        #my $N = $t_seq_len / $width;
        #my $N = $bg_seq_len / $width;
        #my $N = $bg_seq_len;
        my $N = $t_seq_len;

        #
        # standard deviation
        # sqrt(n * p (1 - p))
        # where n = number of trials,
        #	p = probability of success
        #
        #my $S = sqrt($Nb * $C * (1 - $C));
        my $S = sqrt($N * $Rb * (1 - $Rb));

        # expected rate of success (expected number of hits)
        my $Ne = $Nb * $C;
        my $He = $Hb * $C;

        my $Z;
        my $p_value;
        if ($S != 0) {
            $Z = ($Nt - $Ne - 0.5) / $S;
            #$Z = ($Ht - $He - 0.5) / $S;
            $p_value = Statistics::Distributions::uprob($Z);
        } else {
            #
            # Leave as undefined DJA 2010/03/11
            #
            #$Z       = 1e999;    # set to infinite
            #$p_value = 0;
        }

        #if (defined $Z) {
            my $t_gene_hits  = $t_counts->tfbs_gene_count($tf_id);
            my $bg_gene_hits = $bg_counts->tfbs_gene_count($tf_id);
            $results->add_result(
                OPOSSUM::Analysis::ZscoreResult->new(
                    -id           => $tf_id,
                    -t_hits       => $Ht,
                    -bg_hits      => $Hb,
                    -t_rate       => $Rt,
                    -bg_rate      => $Rb,
                    -t_gene_hits  => $t_gene_hits,
                    -bg_gene_hits => $bg_gene_hits,
                    -z_score      => $Z,
                    -p_value      => $p_value
                )
            );
        #}
    }
    $results->param('bg_seq_length', $bg_seq_len);
    $results->param('t_seq_length',  $t_seq_len);

    return $results;
}

1;
