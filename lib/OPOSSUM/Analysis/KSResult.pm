=head1 NAME

OPOSSUM::KSResult.pm - module to hold the result of a Kolmogorov-Smirnov test

=head1 AUTHOR

 Andrew Kwon, based on module by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::KSResult;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $ksr = OPOSSUM::Analysis::KSResult->new(
                -id         => $id,
                -bg_distribution => $bg_distribution,
                -p_value    => $p_value
            );
 Function : Create a new OPOSSUM::Analysis::KSResult object.
 Returns  : An OPOSSUM::Analysis::KSResult object.
 Args     : id          - TFBS ID
            bg_distribution - background distributioin used for KS test
            p_value     - probability that the test set distribution is the same
                            as the background distribution

=cut

sub new
{
    my ($class, %args) = @_;

    my $id = $args{-id};
    if (!$id) {
        carp "must provide ID";
        return;
    }
    my $bg_distribution = $args{-bg_distributon};
    my $p_value    = $args{-p_value};

    my $self = bless {
        -id         => $id,
        -bg_distribution => $bg_distribution,
        -p_value    => $p_value,
    }, ref $class || $class;

    return $self;
}

=head2 id

 Title    : id
 Usage    : $id = $fr->id() or $fr->id($id);
 Function : Get/set the TFBS ID of this result.
 Returns  : TFBS ID string.
 Args     : Optional TFBS ID string.

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

=head2 bg_distribution

 Title    : bg_distribution
 Usage    : $bg_distribution = $ksr->bg_distribution() or $ksr->bg_distribution($bg_distribution);
 Function : Get/set the background distribution used in KS test.
 Returns  : string background distribution (as recognized by R ks.test function).
 Args     : string.

=cut

sub bg_distribution
{
    my ($self, $hits) = @_;

    if (defined $hits) {
        $self->{-bg_distribution} = $hits;
    }

    return $self->{-bg_distribution};
}


1;
