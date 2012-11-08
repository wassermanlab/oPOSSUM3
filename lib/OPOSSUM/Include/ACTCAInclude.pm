use strict;

use OPOSSUM::Include::TCAInclude;


sub find_cluster_sitepairs
{
	my ($anchor_tfbss, $anchor_cluster_id, $tfbss, $cluster_id, $max_dist) = @_;

	my $merged_anchor_tfbss = merge_tfbss($anchor_tfbss, $anchor_cluster_id);
	my $merged_tfbss = merge_tfbss($tfbss, $cluster_id);

	my @sitepairs;
	foreach my $anchor (@$merged_anchor_tfbss) {
		foreach my $tfbs (@$merged_tfbss) {
			if ($tfbs->id eq $anchor->id) {
				# If TF in question is same as anchor TF, only count sites
				# where the TF is to the right of the anchor, to avoid
				# double counting
				next if $tfbs->start <= $anchor->end;
			}

			my $dist;
			if ($tfbs->start() > $anchor->end()) {
				$dist = $tfbs->start() - $anchor->end() - 1;
			} elsif ($anchor->start() > $tfbs->end()) {
				$dist = $anchor->start() - $tfbs->end() - 1;
			} else {
				# do not include TFBSs which overlap anchor TFBS
				next;
			}

			if ($dist <= $max_dist) {
				push @sitepairs, {
					anchor_site		=> $anchor,
					cluster_site	=> $tfbs,
					distance		=> $dist
				};
			}
		}
	}
	
	return @sitepairs ? \@sitepairs : undef;
}

sub merge_tfbss
{
	my ($tfsites, $cluster_id) = @_;

	return if !$tfsites;
	return if scalar(@$tfsites) == 0;

	# I should do some filtering first
	
	my @sorted_sites = sort {$a->start <=> $b->start} @$tfsites;

	my $site1 = $sorted_sites[0];

    $site1->id($cluster_id);

	#
	# Let's not mess around setting and resetting strand info. Strand is
	# meaningless for clusters. Just keep everything on the +ve strand.
	# DJA 2011/06/01
	#
	if ($site1->strand == -1) {
		$site1->strand(1);
		$site1->seq(revcom($site1->seq));
	}

	my @merged_sites;
	push @merged_sites, $site1;

	for (my $i = 1; $i < scalar @sorted_sites; $i++) {
		my $tfsite = $sorted_sites[$i];

        $tfsite->id($cluster_id);

		if ($tfsite->strand == -1) {
			$tfsite->strand(1);
			$tfsite->seq(revcom($tfsite->seq));
		}

		my $prevsite = $merged_sites[$#merged_sites];

		# if overlap, keep the max score
		# merge the two sites
		if (overlap($prevsite, $tfsite)) {
			if ($prevsite->end < $tfsite->end) {
				# merge the sequences
				my $ext_seq = substr(
					$tfsite->seq, $prevsite->end - $tfsite->start + 1
				);

				if ($ext_seq) {
					$prevsite->seq($prevsite->seq . $ext_seq);
				}

				$prevsite->end($tfsite->end);
			}

			if ($tfsite->rel_score > $prevsite->rel_score) {
				$prevsite->rel_score($tfsite->rel_score);
			}

			if ($tfsite->score > $prevsite->score) {
				$prevsite->score($tfsite->score);
			}
		} else {
			$tfsite->id($cluster_id);
			push @merged_sites, $tfsite;
		}
	}

	return @merged_sites ? \@merged_sites : undef;
}

1;
