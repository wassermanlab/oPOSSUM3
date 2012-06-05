=head1 NAME

OPOSSUM::Plot::ScoreVsGC.pm - module to call R to plot Z, Fisher or KS scores
vs. GC content of the TFs.

=head1 SYNOPSIS

 $plotter = OPOSSUM::Plot::ScoreVsGC->new();

 $plotter->run($results, $tf_set, $plot_type, $sd_fold, $filename);

=head1 DESCRIPTION

Given an oPOSSUM analysis result set, plot the relevant score vs. TF GC
content using R. Must have R installed to run.

=head1 AUTHORS

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Plot::ScoreVsGC;

use strict;

use Carp;
use Statistics::R;

=head2 new

 Title    : new
 Usage    : $plotter = OPOSSUM::Plot::ScoreVsGC->new();
 Function : Create a new OPOSSUM::Plot::ScoreVsGC object.
 Returns  : An OPOSSUM::Plot::ScoreVsGC object.
 Args     : None

=cut

sub new
{
    my ($class, @args) = @_;

    #
    # Start R
    #
    my $R = Statistics::R->new();

    unless ($R) {
        carp "Error starting R statistics package\n";
        return;
    }


    my $self = bless {
        _R => $R
    }, ref $class || $class;

    return $self;
}

=head2 plot

 Title    : plot
 Usage    : $plotter->plot(
                $results,
                $tf_set,
                $plot_type,
                $sd_fold,
                $filename,
                $error
            );
 Function : Plot the appropriate score from the oPOSSUM result set.
 Returns  : 
 Args     : results     - A listref of OPOSSUM::Analysis::CombinedResult
                          objects.
            tf_set      - An OPOSSUM::TFSet object.
            plot_type   - String informing which score to plot, either
                          'Z', 'Fisher' or 'KS'.
            sd_fold     - Numeric indicating the standard deviation
                          fold used to determine the threshold.
            filename    - The name of the output plot file. Output file is
                          in PNG format.
            error       - A reference to a scalar error messages are
                          passed back to caller

=cut

sub plot
{
    my ($self, $results, $tf_set, $plot_type, $sd_fold, $filename,
        $error) = @_;

    unless ($results) {
        carp "No result set provided";
        return;
    }

    unless ($tf_set) {
        carp "No TF set provided";
        return;
    }

    unless ($plot_type) {
        carp "No plot type provided";
        return;
    }

    unless ($sd_fold) {
        carp "No SD fold provided";
        return;
    }

    unless ($filename) {
        carp "No output plot file name provided";
        return;
    }

    unless (ref $results eq 'ARRAY' && $results->[0]
        && $results->[0]->isa("OPOSSUM::Analysis::CombinedResult"))
    {
        carp "Results is not an array ref of OPOSSUM::Analysis::CombinedResult"
            . " objects";
        return;
    }

    unless (ref $tf_set && $tf_set->isa("OPOSSUM::TFSet")) {
        carp "TF set is not an OPOSSUM::TFSet object";
        return;
    }

    unless ($plot_type eq 'Z' || $plot_type eq 'Fisher' || $plot_type eq 'KS') {
        carp "Provided plot type is not one of 'Z', 'F' or 'KS'";
        return;
    }

    my $title;
    if ($plot_type eq 'Z') {
        $title = 'Z-score vs. TF profile %GC composition';
    } elsif ($plot_type eq 'Fisher') {
        $title = 'Fisher score vs. TF profile %GC composition';
    } elsif ($plot_type eq 'KS') {
        $title = 'KS-score vs. TF profile %GC composition';
    }

    my @gc;
    my @scores;
    my @names;
    foreach my $result (@$results) {
        my $id = $result->id;

        if ($plot_type eq 'Z') {
            next unless defined $result->zscore();
            push @scores, $result->zscore();
        } elsif ($plot_type eq 'Fisher') {
            next unless defined $result->fisher_p_value();
            push @scores, $result->fisher_p_value();
        } elsif ($plot_type eq 'KS') {
            next unless defined $result->ks_p_value();
            push @scores, $result->ks_p_value();
        }

        my $tf = $tf_set->get_tf($id);

        push @gc, $tf->tag('gc_content') * 100;
        push @names, $tf->name;
    }

    my $mean = _compute_mean(\@scores);
    my $sd = _compute_sd(\@scores, $mean);
    my $threshold = _compute_threshold($mean, $sd, $sd_fold);

    my $max_score = -9999;
    my $min_score = 9999;
    my @gc_above;
    my @scores_above;
    my @names_above;
    #my @gc_below;
    #my @scores_below;
    my $num = scalar @scores;
    foreach my $i (0..$num - 1) {
        next unless defined $scores[$i];

        if ($scores[$i] > $max_score) {
            $max_score = $scores[$i];
        }

        if ($scores[$i] < $min_score) {
            $min_score = $scores[$i];
        }

        if ($scores[$i] >= $threshold) {
            push @gc_above, $gc[$i];
            push @scores_above, $scores[$i];
            push @names_above, $names[$i];
        #} else {
        #    push @gc_below, $gc[$i];
        #    push @scores_below, $scores[$i];
        }
    }

    if ($max_score > 50) {
        $max_score = (int($max_score / 100) + 1) * 100;
    } else {
        $max_score = (int($max_score / 10) + 1) * 10;
    }

    if ($min_score < -50) {
        $min_score = (int($min_score / 100) - 1) * 100;
    } else {
        $min_score = (int($min_score / 10) - 1) * 10;
    }

    my $R = $self->R;

    $R->set('scores', \@scores);
    $R->set('gc', \@gc);
    $R->set('gc_above', \@gc_above);
    $R->set('scores_above', \@scores_above);
    $R->set('names_above', \@names_above);
    $R->set('gc', \@gc);
    $R->set('scores', \@scores);
    $R->set('title', $title);
    $R->set('legend', "mean + $sd_fold * sd");
    $R->set('max_score', $max_score);
    $R->set('min_score', $min_score);
    $R->set('threshold', $threshold);
    $R->set('fname', $filename);

    my $out = $R->run(
        q`xlimit = c(0, 100)`,
        q`ylimit = c(min_score, max_score)`,
        q`png(filename=fname, units="px", width=1024, height=1024, res=NA, pointsize=16, bg="white")`,
        q`plot(gc, scores, cex=0.5, cex.axis=0.9, cex.main=0.8, main=title, xlab="TF profile %GC composition", ylab="Z-score", xlim=xlimit, ylim=ylimit, las=1)`,
        q`text(gc_above, scores_above, labels=names_above, offset=0.5, cex=0.8)`,
        q`abline(h=threshold, col="red", lty=2)`,
        q`legend("topright", legend=c(legend), cex=0.8, col="red", lty=2, bty="n")`,
        q`dev.off()`
    );

    if ($out) {
        $$error = "R plotting call returned: $out";
    }

    return 1;
}

sub R
{
    my $self = shift;

    return $self->{_R};
}

sub DESTROY
{
    my $self = shift;

    $self->R->stop();
}

sub _compute_mean
{
    my ($scores) = @_;

    return 0 unless $scores && $scores->[0];

    #my $n = scalar @$scores;
    my $n = 0;
    my $sum = 0;
    foreach my $score (@$scores) {
        next unless defined $score;
        $sum += $score;
        $n++;
    }

    return $sum / $n;
}

sub _compute_sd
{
    my ($scores, $mean) = @_;

    #my $n = scalar @$scores;
    my $sum = 0;
    my $n = 0;
    foreach my $score (@$scores) {
        next unless defined $score;
        $sum += ($score - $mean)**2;
        $n++;
    }

    return sqrt($sum / $n);
}

sub _compute_threshold
{
    my ($mean, $sd, $sd_fold) = @_;

    return $mean + $sd * $sd_fold;
}

1;
