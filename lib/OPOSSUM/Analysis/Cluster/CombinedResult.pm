=head1 NAME

OPOSSUM::Analysis::Cluster::CombinedResult.pm - module to hold the combined
result of a Fisher and Z-score analysis for a single TFBS

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::CombinedResult;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $result = OPOSSUM::Analysis::Cluster::CombinedResult->new(
                -id                 => $id,
                -t_gene_hits        => $t_gene_hits,
                -bg_gene_hits       => $bg_gene_hits,
                -t_gene_no_hits     => $t_gene_no_hits,
                -bg_gene_no_hits    => $bg_gene_no_hits,
                -t_cluster_hits        => $t_cluster_hits,
                -bg_cluster_hits       => $bg_cluster_hits,
                -t_cluster_rate        => $t_cluster_rate,
                -bg_cluster_rate       => $bg_cluster_rate,
                -zscore             => $zscore,
                -zscore_p_value     => $zscore_pvalue,
                -fisher_p_value     => $fisher_pvalue
            );
 Function : Create a new OPOSSUM::Analysis::Cluster::CombinedResult object.
 Returns  : An OPOSSUM::Analysis::Cluster::CombinedResult object.
 Args     : id              - Result ID, i.e. the TFCluster ID
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
            t_cluster_hits     - Total number of binding sites for this TF cluster
                              in the target set of genes
            bg_cluster_hits    - Total number of binding sites for this TF cluster
                              in the background set of genes
            t_cluster_rate     - Rate of binding site nucleotides vs.
                              non-binding site nucleotides for this TF cluster in
                              the target set of genes
            bg_cluster_rate    - Rate of binding site nucleotides vs.
                              non-binding site nucleotides for this TF cluster in
                              the background set of genes
            zscore          - Z-score measure of significance of number
                              of binding site nucleotides in target vs.
                              background set of genes for this TF cluster.
            zscore_p_value  - Z-score analysis probability measure of
                              significance of number of binding site
                              nucleotides in target set versus the
                              background set of genes for this TF cluster
            fisher_p_value  - Fisher probability measure of number of
                              target genes which had at least one site for
                              this TF cluster vs. the background set

=cut

sub new
{
    my ($class, %args) = @_;

    unless (%args && $args{-id}) {
        carp "Must provide an analysis (TFCluster) ID\n";
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
            binding site of this TF cluster was detected
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
            sites of this TF cluster were detected
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
            binding site of this TF cluster was detected
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
            sites for this TF cluster were detected
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

=head2 t_cluster_hits

 Title    : t_cluster_hits
 Usage    : $cluster_hits = $result->t_cluster_hits();
            $result->t_cluster_hits($cluster_hits);
 Function : Get/set the number of binding sites for this TF cluster which were
            detected in the target set of genes
 Returns  : Target gene TF binding site clusters
 Args     : Optional target gene TF binding site clusters

=cut

sub t_cluster_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-t_cluster_hits} = $hits;
    }

    return $self->{-t_cluster_hits};
}

=head2 bg_cluster_hits

 Title    : bg_cluster_hits
 Usage    : $cluster_hits = $result->bg_cluster_hits();
            $result->bg_cluster_hits($cluster_hits);
 Function : Get/set the number of binding sites for this TF cluster which were
            detected in the background set of genes
 Returns  : Background gene TF binding site clusters
 Args     : Optional background gene TF binding site clusters

=cut

sub bg_cluster_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_cluster_hits} = $hits;
    }

    return $self->{-bg_cluster_hits};
}

=head2 t_cluster_rate

 Title    : t_cluster_rate
 Usage    : $rate = $result->t_cluster_rate();
            $result->t_cluster_rate($rate);
 Function : Get/set the rate of binding site nucleotides for this TF cluster
            which were detected in the target set of genes
 Returns  : Target gene TF binding site cluster nucleotide rate
 Args     : Optional target gene TF binding site cluster nucleotide rate

=cut

sub t_cluster_rate
{
    my ($self, $rate) = @_;

    if (defined $rate) {
        $self->{-t_cluster_rate} = $rate;
    }

    return $self->{-t_cluster_rate};
}

=head2 bg_cluster_rate

 Title    : bg_cluster_rate
 Usage    : $rate = $result->bg_cluster_rate();
            $result->bg_cluster_rate($rate);
 Function : Get/set the rate of binding site nucleotides for this TF cluster
            which were detected in the background set of genes
 Returns  : Background gene TF binding site cluster nucleotide rate
 Args     : Optional background gene TF binding site cluster nucleotide rate

=cut

sub bg_cluster_rate
{
    my ($self, $rate) = @_;

    if (defined $rate) {
        $self->{-bg_cluster_rate} = $rate;
    }

    return $self->{-bg_cluster_rate};
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
    
1;
