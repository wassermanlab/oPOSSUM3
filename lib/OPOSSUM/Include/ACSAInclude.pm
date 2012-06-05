use strict;

use OPOSSUM::Opt::ACSAOpt;

#
# Find all TFBSs proximal (within $max_dist) to the anchor TFBSs.
#
# TFBSs which overlap anchor TFBSs are excluded.
#
# NOTE: A given TFBS could be proximal to more than one anchor. It is counted
# in combination with each anchor (multiple times).
#
# For memory efficiency, changed to return a list of hashes storing the
# minimal amount of information about sites instead of a TFBS::SiteSet.
#
sub proximal_sites
{
    my ($siteset1, $siteset2, $max_dist) = @_;

    my @prox_sites;
    my @sitepairs;

    my $iter1 = $siteset1->Iterator(-sort_by => 'start');
    while (my $tfbs1 = $iter1->next()) {
        my $iter2 = $siteset2->Iterator(-sort_by => 'start');
        while (my $tfbs2 = $iter2->next()) {
            if ($tfbs2->pattern->ID eq $tfbs1->pattern->ID) {
                # If TF in question is same as anchor TF only count sites where
                # the TF is to the right of the anchor to avoid double counting
                next if $tfbs2->start <= $tfbs1->end;
            }

            my $dist;
            if ($tfbs2->start() > $tfbs1->end()) {
                $dist = $tfbs2->start() - $tfbs1->end() - 1;
            } elsif ($tfbs1->start() > $tfbs2->end()) {
                $dist = $tfbs1->start() - $tfbs2->end() - 1;
            } else {
                # do not include TFBSs which overlap anchor TFBS
                next;
            }

            if ($dist <= $max_dist) {
                my $site1_hash = {
                    start       => $tfbs1->start,
                    end         => $tfbs1->end,
                    strand      => $tfbs1->strand,
                    score       => $tfbs1->score,
                    rel_score   => $tfbs1->rel_score,
                    seq         => $tfbs1->seq->seq
                };

                my $site2_hash = {
                    start       => $tfbs2->start,
                    end         => $tfbs2->end,
                    strand      => $tfbs2->strand,
                    score       => $tfbs2->score,
                    rel_score   => $tfbs2->rel_score,
                    seq         => $tfbs2->seq->seq
                };

                push @prox_sites, $site2_hash;

                #my $key = sprintf "%s|%d",
                #    $tfbs2->pattern->ID(), $tfbs2->start, $tfbs2->end;
                push @sitepairs, {
                    anchor_site => $site1_hash,
                    tf_site     => $site2_hash,
                    distance    => $dist
                };
            }
        }
    }

    return (
        @prox_sites ? \@prox_sites : undef,
        @sitepairs ? \@sitepairs : undef
    );
}


sub logprint_siteset
{
    my ($tag, $siteset, %job_args) = @_;

	my $logger = $job_args{-logger};

    $logger->debug("\n$tag:");

    unless ($siteset && $siteset->size() > 0) {
        $logger->debug("No sites");
        return;
    }

    my $iter = $siteset->Iterator(-sort_by => 'start');

    while (my $site = $iter->next()) {
        $logger->debug(
            sprintf("%15s %7d %7d %.3f",
                $site->pattern->name(),
                $site->start(),
                $site->end(),
                $site->score()
            )
        );
    }
}

sub logprint_sitemap
{
    my ($sitemap, %job_args) = @_;

	my $logger = $job_args{-logger};
    $logger->debug("\nSitemap:");

    unless ($sitemap) {
        $logger->debug("No mapped sites");
        return;
    }

    foreach my $key (sort keys %$sitemap) {
        my $pairs = $sitemap->{$key};
        foreach my $pair (@$pairs) {
            my $site1 = $pair->{anchor_site};
            my $site2 = $pair->{tf_site};
            $logger->debug(
                sprintf("%15s %7d %7d | %15s %7d %7d",
                    $site1->pattern->name(),
                    $site1->start(),
                    $site1->end(),
                    $site2->pattern->name(),
                    $site2->start(),
                    $site2->end()
                )
            );
        }
    }
}

sub logprint_sitepairs
{
    my ($sitepairs, %job_args) = @_;

	my $logger = $job_args{-logger};
    $logger->debug("\nSitemap:");

    unless ($sitepairs) {
        $logger->debug("No anchor/TF site pairs");
        return;
    }

    foreach my $pair (@$sitepairs) {
        my $site1 = $pair->{anchor_site};
        my $site2 = $pair->{tf_site};
        $logger->debug(
            sprintf("%-15s %7d %7d %-15s %7d %7d %7d",
                $site1->pattern->name(),
                $site1->start(),
                $site1->end(),
                $site2->pattern->name(),
                $site2->start(),
                $site2->end(),
                $pair->{distance}
            )
        );
    }
}


1;

