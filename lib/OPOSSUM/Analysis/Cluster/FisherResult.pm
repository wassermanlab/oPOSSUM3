=head1 NAME

OPOSSUM::Analysis::Cluster::FisherResult.pm - module to hold the result of a Fisher analysis for
a single TFBS cluster

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::FisherResult;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $fr = OPOSSUM::Analysis::Cluster::FisherResult->new(
                -id         => $id,
                -t_hits     => $t_hits,
                -t_no_hits  => $t_no_hits,
                -bg_hits    => $bg_hits,
                -bg_no_hits => $bg_no_hits,
                -p_value    => $p_value
            );
 Function : Create a new OPOSSUM::Analysis::Cluster::FisherResult object.
 Returns  : An OPOSSUM::Analysis::Cluster::FisherResult object.
 Args     : id          - TFBS cluster ID
            t_hits      - number of genes in the target set which had
                          this TFBS cluster
            t_no_hits   - number of genes in the target set which did
                          not have this TFBS cluster
            bg_hits     - number of genes in the background set which
                          had this TFBS cluster
            bg_no_hits  - number of genes in the background set of which
                          did not have this TFBS cluster
            p_value     - probability that the frequency of hits on the
                          target set versus the background set could have
                          occured by chance

=cut

sub new
{
    my ($class, %args) = @_;

    my $id = $args{-id};
    if (!$id) {
        carp "must provide ID";
        return;
    }
    my $t_hits     = $args{-t_hits};
    my $t_no_hits  = $args{-t_no_hits};
    my $bg_hits    = $args{-bg_hits};
    my $bg_no_hits = $args{-bg_no_hits};
    my $p_value    = $args{-p_value};

    my $self = bless {
        -id         => $id,
        -t_hits     => $t_hits,
        -t_no_hits  => $t_no_hits,
        -bg_hits    => $bg_hits,
        -bg_no_hits => $bg_no_hits,
        -p_value    => $p_value,
    }, ref $class || $class;

    return $self;
}

=head2 id

 Title    : id
 Usage    : $id = $fr->id() or $fr->id($id);
 Function : Get/set the TFBS cluster ID of this result.
 Returns  : TFBS cluster ID string.
 Args     : Optional TFBS cluster ID string.

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
        $self->{-id} = $id;
    }

    return $self->{-id};
}

=head2 p_value

 Title    : p_value
 Usage    : $p_value = $fr->p_value() or $fr->p_value($p_value);
 Function : Get/set the p_value of this result.
 Returns  : real p_value.
 Args     : Optional real p_value.

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
 Usage    : $t_hits = $fr->t_hits() or $fr->t_hits($t_hits);
 Function : Get/set the number of genes for which this TFBS cluster was detected
            in the target set of genes.
 Returns  : integer target hits.
 Args     : optional integer target hits.

=cut

sub t_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-t_hits} = $hits;
    }

    return $self->{-t_hits};
}

=head2 t_no_hits

 Title    : t_no_hits
 Usage    : $t_no_hits = $fr->t_no_hits() or $fr->t_no_hits($t_no_hits);
 Function : Get/set the number of genes for which this TFBS cluster was not
            detected in the target set of genes.
 Returns  : integer target non-hits.
 Args     : optional integer target non-hits.

=cut

sub t_no_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-t_no_hits} = $hits;
    }

    return $self->{-t_no_hits};
}

=head2 bg_hits

 Title    : bg_hits
 Usage    : $bg_hits = $fr->bg_hits() or $fr->bg_hits($bg_hits);
 Function : Get/set the number of genes for which this TFBS cluster was detected
            in the background set of genes.
 Returns  : integer background hits.
 Args     : optional integer background hits.

=cut

sub bg_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_hits} = $hits;
    }

    return $self->{-bg_hits};
}

=head2 bg_no_hits

 Title    : bg_no_hits
 Usage    : $bg_no_hits = $fr->bg_no_hits() or $fr->bg_no_hits($bg_no_hits);
 Function : Get/set the number of genes for which this TFBS cluster was not
            detected in the background set of genes.
 Returns  : integer background non-hits.
 Args     : optional integer background non-hits.

=cut

sub bg_no_hits
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_no_hits} = $hits;
    }

    return $self->{-bg_no_hits};
}

1;
