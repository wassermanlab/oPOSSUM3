=head1 NAME

OPOSSUM::CombinedResult.pm - module to hold the combined result of a
Fisher, Z-score and KS analysis for a single TFBS

=head1 MODIFICATIONS

 Andrew Kwon, Nov. 16, 2011
 - added KS parts
 
=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::CombinedResult;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $result = OPOSSUM::Analysis::CombinedResult->new(
                -id                 => $id,
                -t_gene_hits        => $t_gene_hits,
                -bg_gene_hits       => $bg_gene_hits,
                -t_gene_no_hits     => $t_gene_no_hits,
                -bg_gene_no_hits    => $bg_gene_no_hits,
                -t_tfbs_hits        => $t_tfbs_hits,
                -bg_tfbs_hits       => $bg_tfbs_hits,
                -t_tfbs_rate        => $t_tfbs_rate,
                -bg_tfbs_rate       => $bg_tfbs_rate,
                -zscore             => $zscore,
                -zscore_p_value     => $zscore_pvalue,
                -fisher_p_value     => $fisher_pvalue,
                -ks_p_value         => $ks_pvalue
            );
 Function : Create a new OPOSSUM::Analysis::CombinedResult object.
 Returns  : An OPOSSUM::Analysis::CombinedResult object.
 Args     : id              - Result ID, i.e. the TF ID
            t_gene_hits     - Number of genes in the target set which had
                              at least one binding site for this TF
            t_gene_no_hits  - Number of genes in the target set which did
                              not have at least one binding site for this
                              TF
            bg_gene_hits    - Number of genes in the background set which
                              had at least one binding site for this TF
            bg_gene_no_hits - Number of genes in the background set of
                              which did not have at least one binding site
                              for this TF
            t_tfbs_hits     - Total number of binding sites for this TF in
                              the target set of genes
            bg_tfbs_hits    - Total number of binding sites for this TF in
                              the background set of genes
            t_tfbs_rate     - Rate of binding site nucleotides vs.
                              non-binding site nucleotides for this TF in
                              the target set of genes
            bg_tfbs_rate    - Rate of binding site nucleotides vs.
                              non-binding site nucleotides for this TF in
                              the background set of genes
            zscore          - Z-score measure of significance of number
                              of binding site nucleotides in target vs.
                              background set of genes for this TF.
            zscore_p_value  - Z-score analysis probability measure of
                              significance of number of binding site
                              nucleotides in target set versus the
                              background set of genes for this TF
            fisher_p_value  - Fisher probability measure of number of
                              target genes which had at least one site for
                              this TF vs. the background set
            ks_p_value      - K-S test p-value for the test and background
                              distributions of peak-to-site distances

=cut

sub new
{
    my ($class, %args) = @_;

    unless (%args && $args{-id}) {
        carp "Must provide an analysis (TF) ID\n";
        return;
    }

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 id

 Title    : id
 Usage    : $id = $result->id();
            $result->id($id);
 Function : Get/set the (TF) ID of this result.
 Returns  : Result (TF) ID
 Args     : Optional result (TF) ID

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
        $self->{-id} = $id;
    }

    return $self->{-id};
}

=head2 t_gene_hits

 Title    : t_gene_hits
 Usage    : $gene_hits = $result->t_gene_hits();
            $result->t_gene_hits($gene_hits);
 Function : Get/set the number of target genes for which at least one
            binding site of this TF was detected
 Returns  : Target gene hits
 Args     : Optional target gene hits

=cut

sub t_gene_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-t_gene_hits} = $hits;
    }

    return $self->{-t_gene_hits};
}

=head2 t_gene_no_hits

 Title    : t_gene_no_hits
 Usage    : $gene_no_hits = $result->t_gene_no_hits();
            $result->t_gene_no_hits($gene_no_hits);
 Function : Get/set the number of target genes for which no binding
            sites of this TF were detected
 Returns  : Target gene non-hits.
 Args     : Optional target gene non-hits.

=cut

sub t_gene_no_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-t_gene_no_hits} = $hits;
    }

    return $self->{-t_gene_no_hits};
}

=head2 bg_gene_hits

 Title    : bg_gene_hits
 Usage    : $gene_hits = $result->bg_gene_hits();
            $result->bg_gene_hits($gene_hits);
 Function : Get/set the number of target genes for which at least one
            binding site of this TF was detected
 Returns  : Background gene hits.
 Args     : Optional background gene hits.

=cut

sub bg_gene_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_gene_hits} = $hits;
    }

    return $self->{-bg_gene_hits};
}

=head2 bg_gene_no_hits

 Title    : bg_gene_no_hits
 Usage    : $gene_no_hits = $result->bg_gene_no_hits();
            $result->bg_gene_no_hits($gene_no_hits);
 Function : Get/set the number of background genes for which no binding
            sites for this TF were detected
 Returns  : Background gene non-hits.
 Args     : Optional background gene non-hits.

=cut

sub bg_gene_no_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_gene_no_hits} = $hits;
    }

    return $self->{-bg_gene_no_hits};
}

=head2 t_tfbs_hits

 Title    : t_tfbs_hits
 Usage    : $tfbs_hits = $result->t_tfbs_hits();
            $result->t_tfbs_hits($tfbs_hits);
 Function : Get/set the number of binding sites for this TF which were
            detected in the target set of genes
 Returns  : Target gene TF binding sites
 Args     : Optional target gene TF binding sites

=cut

sub t_tfbs_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-t_tfbs_hits} = $hits;
    }

    return $self->{-t_tfbs_hits};
}

=head2 bg_tfbs_hits

 Title    : bg_tfbs_hits
 Usage    : $tfbs_hits = $result->bg_tfbs_hits();
            $result->bg_tfbs_hits($tfbs_hits);
 Function : Get/set the number of binding sites for this TF which were
            detected in the background set of genes
 Returns  : Background gene TF binding sites
 Args     : Optional background gene TF binding sites

=cut

sub bg_tfbs_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_tfbs_hits} = $hits;
    }

    return $self->{-bg_tfbs_hits};
}

=head2 t_tfbs_rate

 Title    : t_tfbs_rate
 Usage    : $rate = $result->t_tfbs_rate();
            $result->t_tfbs_rate($rate);
 Function : Get/set the rate of binding site nucleotides for this TF
            which were detected in the target set of genes
 Returns  : Target gene TF binding site nucleotide rate
 Args     : Optional target gene TF binding site nucleotide rate

=cut

sub t_tfbs_rate
{
    my ($self, $rate) = @_;

    if (defined $rate) {
        $self->{-t_tfbs_rate} = $rate;
    }

    return $self->{-t_tfbs_rate};
}

=head2 bg_tfbs_rate

 Title    : bg_tfbs_rate
 Usage    : $rate = $result->bg_tfbs_rate();
            $result->bg_tfbs_rate($rate);
 Function : Get/set the rate of binding site nucleotides for this TF
            which were detected in the background set of genes
 Returns  : Background gene TF binding site nucleotide rate
 Args     : Optional background gene TF binding site nucleotide rate

=cut

sub bg_tfbs_rate
{
    my ($self, $rate) = @_;

    if (defined $rate) {
        $self->{-bg_tfbs_rate} = $rate;
    }

    return $self->{-bg_tfbs_rate};
}

=head2 zscore

 Title    : zscore
 Usage    : $zscore = $result->zscore();
            $result->zscore($zscore);
 Function : Get/set the Z-score of this result.
 Returns  : Real z-score.
 Args     : Optional real z-score.

=cut

sub zscore
{
    my ($self, $zscore) = @_;

    if (defined $zscore) {
        $self->{-zscore} = $zscore;
    }

    return $self->{-zscore};
}

=head2 zscore_p_value

 Title    : zscore_p_value
 Usage    : $p_value = $result->zscore_p_value();
            $result->zscore_p_value($p_value);
 Function : Get/set the Z-score p_value of this result.
 Returns  : Real p_value.
 Args     : Optional real p_value.

=cut

sub zscore_p_value
{
    my ($self, $p_value) = @_;

    if (defined $p_value) {
        $self->{-zscore_p_value} = $p_value;
    }

    return $self->{-zscore_p_value};
}

=head2 fisher_p_value

 Title    : fisher_p_value
 Usage    : $p_value = $result->fisher_p_value();
            $result->fisher_p_value($p_value);
 Function : Get/set the Fisher's exact test p_value of this result.
 Returns  : Real p_value.
 Args     : Optional real p_value.

=cut

sub fisher_p_value
{
    my ($self, $p_value) = @_;

    if (defined $p_value) {
        $self->{-fisher_p_value} = $p_value;
    }

    return $self->{-fisher_p_value};
}

#
# Synonym for fisher_p_value()
#
sub fisher_score
{
    my ($self, $score) = @_;

    return $self->fisher_p_value($score);
}
   
=head2 ks_p_value

 Title    : ks_p_value
 Usage    : $p_value = $result->ks_p_value();
            $result->ks_p_value($p_value);
 Function : Get/set the KS test p_value of this result.
 Returns  : Real p_value.
 Args     : Optional real p_value.

=cut

sub ks_p_value
{
    my ($self, $p_value) = @_;

    if (defined $p_value) {
        $self->{-ks_p_value} = $p_value;
    }

    return $self->{-ks_p_value};
}
 
1;
