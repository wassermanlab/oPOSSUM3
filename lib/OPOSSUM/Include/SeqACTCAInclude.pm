use strict;

use OPOSSUM::Opt::ACSAOpt;

use OPOSSUM::Include::SeqInclude;
use OPOSSUM::Include::ACTCAInclude;

use constant BG_COLOR_CLASS => 'bgc_seq_actca';

#
# Search seqs with the anchoring TFBS cluster and all the TFs in the TF set.
#
# Return two hashrefs:
#
# The first hashref refers to a hash of hashes of sites indexed on seq ID
# and TF ID. In the case where there were no TFBSs for a seq/TF combo,
# the value of the hash is undef.
#
# The second hasref refers to a hash of listrefs of sequence IDs, indexed
# on TF ID. Only seq IDs for which there were TFBSs for this TF are stored
# in the list.
#
# i.e.:
# $hashref1->{$seq_id}->{$tf_id} = $siteset
# $hashref2->{$tf_id} = \@seq_ids
#
sub anchored_tf_cluster_set_search_seqs
{
    my (
        $tf_cluster_set, $tf_set, $anchor_cluster, $max_site_dist,
		$seq_id_seqs, $threshold, %job_args
    ) = @_;

	my $logger = $job_args{-logger};

    my $anchor_cluster_tf_ids = $anchor_cluster->tf_ids();
    
    my $cluster_ids = $tf_cluster_set->ids();

    my @seq_ids = keys %$seq_id_seqs;
    
    # hashs for final results
    my %cluster_seq_sites;
    my %cluster_seq_sitepairs;
    
    # Go through each sequence and search for the TFs in the anchor cluster
    # then go through each cluster and search for those TFs
    # merge the site hits based on cluster
    # using the merged sites, find proximal site pairs
    foreach my $seq_id (@seq_ids) {
        $logger->debug("\nSequence: $seq_id\n");

        my $seq = $seq_id_seqs->{$seq_id};
        
        # retrieve the TFs in the anchor cluster
        my $anchor_siteset = TFBS::SiteSet->new();
        foreach my $anchor_tf_id (@$anchor_cluster_tf_ids) {
            my $anchor_matrix = $tf_set->get_matrix($anchor_tf_id);
            next if !$anchor_matrix;
            my $anchor_pwm;
            if ($anchor_matrix->isa("TFBS::Matrix::PFM")) {
                $anchor_pwm = $anchor_matrix->to_PWM();
            } else {
                $anchor_pwm = $anchor_matrix;
            }
            my $anchor_tf_siteset = $anchor_pwm->search_seq(
                -seqobj     => $seq,
                -threshold  => $threshold
            );
            
            # add to anchor_sites
            $anchor_siteset->add_siteset($anchor_tf_siteset);
        }
        
        next if $anchor_siteset->size == 0;
        
        my $merged_anchor_sites = merge_cluster_siteset(
            $anchor_siteset, $seq_id, $anchor_cluster->id);
        
        foreach my $cl_id (@$cluster_ids) {
            #next if $cl_id eq $anchor_cluster->id;

            my $cluster = $tf_cluster_set->get_tf_cluster($cl_id);
            my $cl_tf_ids = $cluster->tf_ids;
            my $cluster_siteset = TFBS::SiteSet->new();
            foreach my $cl_tf_id (@$cl_tf_ids) {
                my $tf_matrix = $tf_set->get_matrix($cl_tf_id);
                if (!$tf_matrix) {
                    next;
                }
                my $tf_pwm;
                if ($tf_matrix->isa("TFBS::Matrix::PFM")) {
                    $tf_pwm = $tf_matrix->to_PWM();
                } else {
                    $tf_pwm = $tf_matrix;
                }
                my $tf_siteset = $tf_pwm->search_seq(
                    -seqobj     => $seq,
                    -threshold  => $threshold
                );
                
                # add to cluster_sites
                $cluster_siteset->add_siteset($tf_siteset);
            }
            
            next if ($cluster_siteset->size == 0);
            
            my $merged_cl_sites = merge_cluster_siteset(
                $cluster_siteset, $seq_id, $cl_id
            );
            
            # now that we have the merged anchor sites and cluster sites,
            # find which ones are proximal pairs
            my ($prox_sites, $sitepairs) = proximal_sites(
                $merged_anchor_sites, $merged_cl_sites, $max_site_dist
            );

            if (DEBUG) {
                logprint_sites("Proximal sites", $prox_sites);
                logprint_sitepairs($sitepairs);
            }

            if ($prox_sites && scalar @$prox_sites > 0) {
                $cluster_seq_sites{$cl_id}->{$seq_id} = $prox_sites;
                $cluster_seq_sitepairs{$cl_id}->{$seq_id} = $sitepairs;
            }
        }
    }

    my $retval1 = %cluster_seq_sites ? \%cluster_seq_sites : undef;
    my $retval2 = %cluster_seq_sitepairs ? \%cluster_seq_sitepairs : undef;

    return ($retval1, $retval2);
}


sub compute_site_cluster_counts
{
    my ($cluster_set, $seq_ids, $cluster_seq_sites) = @_;

    my $cl_ids  = $cluster_set->ids();

    my $counts = OPOSSUM::Analysis::Cluster::Counts->new(
        -gene_ids       => $seq_ids,
        -cluster_ids    => $cl_ids
    );

    foreach my $seq_id (@$seq_ids) {
        foreach my $cl_id (@$cl_ids) {
            my $sites = $cluster_seq_sites->{$cl_id}->{$seq_id};

            if ($sites) {
                # Note set size could be 0.
                $counts->gene_cluster_count($seq_id, $cl_id, scalar(@$sites));
                my $length = 0;
                foreach my $site (@$sites) {
                    $length += length($site->seq); # site => OPOSSUM::ConservedTFBS
                }
                $counts->gene_cluster_length($seq_id, $cl_id, $length);
            } else {
                $counts->gene_cluster_count($seq_id, $cl_id, 0);
                $counts->gene_cluster_length($seq_id, $cl_id, 0);
            }
        }
    }

    return $counts;
}

#
# merges the TFBS::Sites in TFBS::SiteSet based on position
# converted to OPOSSUM::ConservedTFBS objects
#
sub merge_cluster_siteset
{
    my ($tf_siteset, $seq_id, $cluster_id) = @_;
    
    return if !$tf_siteset or $tf_siteset->size == 0;
    
    #my $filtered_siteset = filter_overlapping_sites($tf_siteset);

    #logprint_siteset("Filtered TF siteset", $filtered_siteset);

    #return if !$filtered_siteset || $filtered_siteset->size == 0;

    #my $it = $filtered_siteset->Iterator(-sort_by => 'start');
    my $it = $tf_siteset->Iterator(-sort_by => 'start');
    my $tfsite = $it->next;

    #
    # Let's not mess around setting and resetting strand info. Strand is
    # meaningless for clusters. Just keep everything on the +ve strand.
    # DJA 2011/05/24
    #
    my $seq = $tfsite->seq->seq;
    if ($tfsite->strand == -1) {
        $seq = revcom($seq);
    }

    my $prev_ctfbs = OPOSSUM::ConservedTFBS->new(
        -tf_id      => $cluster_id,
        -gene_id    => $seq_id,
        -start      => $tfsite->start,
        -end        => $tfsite->end,
        -strand     => 1,
        -score      => $tfsite->score,
        -rel_score  => $tfsite->rel_score,
        -seq        => $seq
    );
    
    my @merged_sites;
    push @merged_sites, $prev_ctfbs;
    
    while ($tfsite = $it->next)
    {
        my $seq = $tfsite->seq->seq;
        if ($tfsite->strand == -1) {
            $seq = revcom($seq);
        }

        my $ctfbs = OPOSSUM::ConservedTFBS->new(
            -tf_id      => $cluster_id,
            -gene_id    => $seq_id,
            -start      => $tfsite->start,
            -end        => $tfsite->end,
            -strand     => 1,
            -score      => $tfsite->score,
            -rel_score  => $tfsite->rel_score,
            -seq        => $seq
        );
        
        $prev_ctfbs = $merged_sites[$#merged_sites];
        
        # if overlap, keep the max score
        # merge the two sites
        if (overlap($prev_ctfbs, $ctfbs))
        {
            if ($prev_ctfbs->end < $ctfbs->end) {
                #$prev_ctfbs->end($ctfbs->end);
                #
                # merge the sequences
                # first, check the strands of the sites
                # if negative, reverse complement
                # I should only do this if they are overlapping
                #if ($prev_ctfbs->strand != $ctfbs->strand) {
                #    if ($prev_ctfbs->strand == -1) {
                #        my $seq = Bio::Seq->new(-seq => $prev_ctfbs->seq);
                #        $prev_ctfbs->seq($seq->revcom->seq);
                #    } else {
                #        my $seq = Bio::Seq->new(-seq => $ctfbs->seq);
                #        $ctfbs->seq($seq->revcom->seq);
                #    }
                #}
                
                #my $ext_seq = substr($ctfbs->seq, $prev_ctfbs->end - $ctfbs->start);
                # DJA 2011/05/24
                my $ext_seq = substr(
                    $seq, $prev_ctfbs->end - $ctfbs->start + 1
                );

                #
                # Have to check this. If ends of ctfbs and prev_ctfbs are the
                # same, then ext_seq is undef.
                # DJA 2011/05/24
                #
                if ($ext_seq) {
                    $prev_ctfbs->seq($prev_ctfbs->seq . $ext_seq);
                }

                $prev_ctfbs->end($ctfbs->end);
            }

            if ($ctfbs->rel_score > $prev_ctfbs->rel_score) {
                    $prev_ctfbs->rel_score($ctfbs->rel_score);
            }

            if ($ctfbs->score > $prev_ctfbs->score) {
                    $prev_ctfbs->score($tfsite->score);
            }
        } else {
            push @merged_sites, $ctfbs;
        }
    }
    
    return \@merged_sites;
}


#
# Find all TFBSs proximal (within $max_dist) to the anchor TFBSs.
#
# TFBSs which overlap anchor TFBSs are excluded.
#
# NOTE: A given TFBS could be proximal to more than one anchor. It is counted
# in combination with each anchor (multiple times).
#
sub proximal_sites
{
	my ($sites1, $sites2, $max_dist) = @_;

	return if !$sites1 or !$sites2;
	return if scalar @$sites1 == 0;
	return if scalar @$sites2 == 0;

	my @prox_sites;

	my @sitepairs;

	# $sites1 and $sites2 should have already been sorted
	for (my $i = 0; $i < scalar @$sites1; $i++) {
		my $tfbs1 = $$sites1[$i];

		for (my $j = 0; $j < scalar @$sites2; $j++) {
			my $tfbs2 = $$sites2[$j];

			if ($tfbs2->id eq $tfbs1->id) {
				# If TF in question is the same as the anchor TF, only count
				# sites where the TF is to the right of the anchor to avoid
				# double counting
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
				push @prox_sites, $tfbs2;

				push @sitepairs, {
					anchor_site => $tfbs1,
					cluster_site => $tfbs2,
					distance => $dist
				};
			}
		}
	}

	return (\@prox_sites, \@sitepairs);
}


1;

