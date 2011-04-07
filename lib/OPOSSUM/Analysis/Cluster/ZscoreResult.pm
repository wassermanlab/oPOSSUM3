=head1 NAME

OPOSSUM::Analysis::Cluster::ZscoreResult.pm - module to hold the result of a
Zscore analysis for a single TFBS cluster

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::ZscoreResult;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $zr = OPOSSUM::Analysis::Cluster::ZscoreResult->new(
                -id             => $id,
                -t_hits         => $t_hits,
                -bg_hits        => $bg_hits,
                -t_rate         => $t_rate,
                -bg_rate        => $bg_rate,
                -t_gene_hits    => $t_gene_hits,
                -bg_gene_hits   => $bg_gene_hits,
                -z_score        => $z_score,
                -p_value        => $p_value);
 Function : Create a new OPOSSUM::Analysis::Cluster::ZscoreResult object.
 Returns  : An OPOSSUM::Analysis::Cluster::ZscoreResult object.
 Args     : id      - cluster ID
            t_hits  - number of times this cluster was detected in the
                      target set of genes
            bg_hits	- number of times this cluster was detected in the
                      background set of genes
            t_rate	- rate that this cluster's nucleotides were detected
                      in the target set of genes
            bg_rate	- rate that this cluster's nucleotides were detected
                      in the background set of genes
            t_gene_hits	- number of genes in the target set for which this
                      cluster was detected
            bg_gene_hits - number of genes in the background set for which
                      this cluster was detected
            z_score	- measure of significance of the target rate
                      versus the background rate
            p_value	- probability that this rate could occur by chance

=cut

sub new
{
    my ($class, %args) = @_;

    my $id = $args{-id};
    if (!$id) {
        carp "must provide ID";
        return;
    }
    my $z_score      = $args{-z_score};
    my $p_value      = $args{-p_value};
    my $t_hits       = $args{-t_hits};
    my $bg_hits      = $args{-bg_hits};
    my $t_rate       = $args{-t_rate};
    my $bg_rate      = $args{-bg_rate};
    my $t_gene_hits  = $args{-t_gene_hits};
    my $bg_gene_hits = $args{-bg_gene_hits};

    my $self = bless {
        -id           => $id,
        -z_score      => $z_score,
        -p_value      => $p_value,
        -t_hits       => $t_hits,
        -bg_hits      => $bg_hits,
        -t_rate       => $t_rate,
        -bg_rate      => $bg_rate,
        -t_gene_hits  => $t_gene_hits,
        -bg_gene_hits => $bg_gene_hits
    }, ref $class || $class;

    return $self;
}

=head2 id

 Title    : id
 Usage    : $id = $zr->id() or $zr->id($id);
 Function : Get/set the cluster ID of this result.
 Returns  : cluster ID string.
 Args     : Optional cluster ID string.

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
        $self->{-id} = $id;
    }

    return $self->{-id};
}

=head2 z_score

 Title    : z_score
 Usage    : $z_score = $zr->z_score() or $zr->z_score($z_score);
 Function : Get/set the z-score of this result.
 Returns  : real z-score.
 Args     : Optional real z-score.

=cut

sub z_score
{
    my ($self, $z_score) = @_;

    if (defined $z_score) {
        $self->{-z_score} = $z_score;
    }

    return $self->{-z_score};
}

=head2 p_value

 Title    : p_value
 Usage    : $p_value = $zr->p_value() or $zr->p_value($p_value);
 Function : Get/set the p-value of this result.
 Returns  : real p-value.
 Args     : Optional real p-value.

=cut

sub p_value
{
    my ($self, $p_value) = @_;

    if (defined $p_value) {
        $self->{-p_value} = $p_value;
    }

    return $self->{-p_value};
}

=head2 t_hits

 Title    : t_hits
 Usage    : $t_hits = $zr->t_hits() or $zr->t_hits($t_hits);
 Function : Get/set the target hits of this result.
 Returns  : integer target hits.
 Args     : Optional integer target hits.

=cut

sub t_hits
{
    my ($self, $t_hits) = @_;

    if (defined $t_hits) {
        $self->{-t_hits} = $t_hits;
    }

    return $self->{-t_hits};
}

=head2 bg_hits

 Title    : bg_hits
 Usage    : $bg_hits = $zr->bg_hits() or $zr->bg_hits($bg_hits);
 Function : Get/set the background hits of this result.
 Returns  : integer background hits.
 Args     : Optional integer background hits.

=cut

sub bg_hits
{
    my ($self, $bg_hits) = @_;

    if (defined $bg_hits) {
        $self->{-bg_hits} = $bg_hits;
    }

    return $self->{-bg_hits};
}

=head2 t_rate

 Title    : t_rate
 Usage    : $t_rate = $zr->t_rate() or $zr->t_rate($t_rate);
 Function : Get/set the target rate of this result.
 Returns  : real target rate.
 Args     : Optional real target rate.

=cut

sub t_rate
{
    my ($self, $t_rate) = @_;

    if (defined $t_rate) {
        $self->{-t_rate} = $t_rate;
    }

    return $self->{-t_rate};
}

=head2 bg_rate

 Title    : bg_rate
 Usage    : $bg_rate = $zr->bg_rate() or $zr->bg_rate($bg_rate);
 Function : Get/set the background rate of this result.
 Returns  : real background rate.
 Args     : Optional real background rate.

=cut

sub bg_rate
{
    my ($self, $bg_rate) = @_;

    if (defined $bg_rate) {
        $self->{-bg_rate} = $bg_rate;
    }

    return $self->{-bg_rate};
}

=head2 t_gene_hits

 Title    : t_gene_hits
 Usage    : $t_gene_hits = $zr->t_gene_hits()
	    or $zr->t_gene_hits($t_gene_hits);
 Function : Get/set the target gene hits of this result.
 Returns  : integer target gene hits.
 Args     : Optional integer target gene hits.

=cut

sub t_gene_hits
{
    my ($self, $t_gene_hits) = @_;

    if (defined $t_gene_hits) {
        $self->{-t_gene_hits} = $t_gene_hits;
    }

    return $self->{-t_gene_hits};
}

=head2 bg_gene_hits

 Title    : bg_gene_hits
 Usage    : $bg_gene_hits = $zr->bg_gene_hits()
	    or $zr->bg_gene_hits($bg_gene_hits);
 Function : Get/set the background gene hits of this result.
 Returns  : integer background gene hits.
 Args     : Optional integer background gene hits.

=cut

sub bg_gene_hits
{
    my ($self, $bg_gene_hits) = @_;

    if (defined $bg_gene_hits) {
        $self->{-bg_gene_hits} = $bg_gene_hits;
    }

    return $self->{-bg_gene_hits};
}

1;
