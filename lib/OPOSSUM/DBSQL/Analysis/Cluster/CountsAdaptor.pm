
=head1 NAME

OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor - Adaptor for MySQL queries to retrieve
and store gene TFBS cluster counts.

=head1 SYNOPSIS

$aca = $db_adaptor->get_AnalysisCountsAdaptor();

=head1 DESCRIPTION

In order to facilitate fast retrieval of TFBS cluster counts from the oPOSSUM database
several count sets were pre-computed using discrete values for PWM (PSSM)
thresholds, conservation levels, and upstream/downstream search regions.
This adaptor provides an interface to both the pre-computed counts as well as
the raw data used to compute counts and retrieves the information into a
OPOSSUM::Analysis::Cluster::Counts object which can be passed directly to the
OPOSSUM::Analysis::Fisher and OPOSSUM::Analysis::Zscore modules.

=head1 CHANGE HISTORY

 AK 2010/12/02
 - added operon handling routines included in Analysis::CountsAdaptor
 
=head1 AUTHOR

 Andrew Kwon, adapting code by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;
use Bio::Seq;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

sub new
{
    my ($class, $dbobj) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new($dbobj);

    return $self;
}

=head2 fetch_counts

 Title    : fetch_counts
 Usage    : $counts = $aca->fetch_counts(
                -conservation_level     => $clevel,
                -threshold_level        => $thlevel,
                -search_region_level    => $srlevel,
                -gene_ids               => $gids,
                -cluster_ids            => $cids,
				-operon_gene_ids		=> $operon_gids,
				-has_operon				=> $has_operon
            );
 Function : Fetch the pre-computed counts of TFBS hits on the gene
            sequences from the DB and build a
            OPOSSUM::Analysis::Cluster::Counts object.
 Returns  : An OPOSSUM::Analysis::Cluster::Counts object.
 Args	  : conservation_level	- counts are retrieved for TFBSs which
                                  were found in conserved regions with
                                  at least the conservation identified
                                  by this conservation level
            threshold_level	    - counts are retrieved for TFBS hits
                                  which had at least the threshold score
                                  identified by this threshold level
            search_region_level	- counts are retrieved for TFBSs which
                                  fell within the search region identified
                                  by this search region level
            gene_ids	        - optionally restrict the counts to only
                                  those gene identified by these gene IDs
            cluster_ids		        - optionally restrict the counts to only
                                  those TF profiles identified by these
                                  TF IDs
			operon_gene_ids		- hash with key = operon gene id, val = first
								  operon gene id. used to replicate the counts
								  for the first gene with other genes belonging
								  to the same operon
			has_operon			- indicates whether this species has operons
								  optional operon_gene_ids are specified
=cut

sub fetch_counts
{
    my ($self, %args) = @_;

    my $cons_level   = $args{-conservation_level};
    my $thresh_level = $args{-threshold_level};
    my $sr_level     = $args{-search_region_level};
    my $gids         = $args{-gene_ids};
    my $cids         = $args{-cluster_ids};
	my $operon_gids	 = $args{-operon_gene_ids};
	my $has_operon	 = $args{-has_operon};

    if (!$cons_level || !$thresh_level || !$sr_level) {
        carp "must provide conservation_level, threshold_level and"
            . " search_region_level\n";
        return;
    }

    #my $crla = $self->db->get_ConservedRegionLengthAdaptor();
    #if (!$crla) {
    #    carp "error getting ConservedRegionLengthAdaptor\n";
    #    return;
    #}

    my $sql = qq{
		select gene_id, cluster_id, count, sum_length
		from tfbs_cluster_counts
		where conservation_level = $cons_level
		and threshold_level = $thresh_level
		and search_region_level = $sr_level
    };

    #print STDERR $sql ."\n";

    #
    # If no gene IDs are passed in, get all the gene IDs from the database
    # in order to initialize counts object properly, i.e. ensure that
    # counts contain 0's for all gene/TF combinations even if there is
    # no record in the database tfbs_counts table (0 values not stored
    # for efficiency).
    #
	# AK: This one would also need operon treatment.
	#
    unless ($gids) {
        my $ga = $self->db->get_GeneAdaptor();

        if (!$ga) {
            carp "error getting GeneAdaptor\n";
            return;
        }

        $gids = $ga->fetch_gene_ids();
		
		if ($has_operon) {
			my $oa = $self->db->get_OperonAdaptor();
			if (!$oa) {
				carp "error getting OperonAdaptor\n";
				return;
			}
			my $operons = $oa->fetch_where();
			foreach my $operon (@$operons) {
				my $fgene = $operon->fetch_first_gene();
				foreach my $opgene (@{$operon->genes}) {
					$$operon_gids{$opgene->id} = $fgene->id;
					#	if $opgene->id != $fgene->id;
				}
			}
		}
    }
	
    #
	# Set t_gids to all IDs of genes in the input gene IDs list which are
    # either not in an operon or are the first gene in an operon containing
    # an input gene.
    #
	my @t_gids;
    my %gid_inc;
	if ($operon_gids) {
		foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            if ($first_gid) {
                unless ($gid_inc{$first_gid}) {
                    push @t_gids, $first_gid;
                    $gid_inc{$first_gid} = 1;
                }
            } else {
                unless ($gid_inc{$gid}) {
                    push @t_gids, $gid;
                    $gid_inc{$gid} = 1;
                }
            }
		}
	} else {
		@t_gids = @$gids;
	}
        
	$sql .= " and gene_id in (";
    $sql .= join ',', @t_gids;
    $sql .= ")";
	
    if ($cids && $cids->[0]) {
        $sql .= " and cluster_id in ('";
        $sql .= join "','", @$cids;
        $sql .= "')";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TFBS cluster counts with\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TFBS cluster counts with\n$sql\n" . $sth->errstr;
        return;
    }

    my $t_counts = OPOSSUM::Analysis::Cluster::Counts->new(
        -gene_ids => \@t_gids,
        -cluster_ids => $cids
    );

    while (my @row = $sth->fetchrow_array) {
        $t_counts->gene_cluster_count($row[0], $row[1], $row[2]);
		$t_counts->gene_cluster_length($row[0], $row[1], $row[3])
    }
    $sth->finish;
	
	#
	# Now, add the operon gene counts
	#
    my $counts;
	if ($operon_gids) {
        $counts = OPOSSUM::Analysis::Cluster::Counts->new(
            -gene_ids		=> $gids,
            -cluster_ids	=> $cids
        );

        foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            unless ($first_gid) {
                $first_gid = $gid;
            }

            foreach my $cid (@$cids) {
                my $count = $t_counts->gene_cluster_count($first_gid, $cid);
				my $length = $t_counts->gene_cluster_length($first_gid, $cid);
                $counts->gene_cluster_count($gid, $cid, $count);
				$counts->gene_cluster_length($gid, $cid, $length);
            }
        }
    } else {
        $counts = $t_counts;
	}

    $counts->param('conservation_level',  $cons_level);
    $counts->param('threshold_level',     $thresh_level);
    $counts->param('search_region_level', $sr_level);

    return $counts;
}

=head2 fetch_custom_counts

 Title    : fetch_custom_counts
 Usage    : $counts = $aca->fetch_custom_counts(
                -conservation_level     => $clevel,
                -threshold              => $thresh,
                -upstream_bp            => $upbp,
                -downstream_bp          => $downbp,
                -gene_ids               => $gids,
                -cluster_ids                 => $cids);
 Function : Compute the counts of TFBS hits on the gene sequences
            from the DB and build an OPOSSUM::Analysis::Cluster::Counts object.
 Returns  : An OPOSSUM::Analysis::Cluster::Counts object.
 Args	  : conservation_level  - counts are retrieved for TFBSs which
                                  were found in conserved regions with
                                  at least the conservation identified
                                  by this conservation level
            threshold           - counts are retrieved for TFBS hits
                                  which had at least this threshold score
            upstream_bp         - counts are retrieved for TFBSs which
                                  fell within this number of bp upstream
                                  of the TSS
            downstream_bp       - counts are retrieved for TFBSs which
                                  fell within this number of bp downstream
                                  of the TSS
            gene_ids            - optionally restrict the counts to only
                                  those genes identified by
                                  these gene IDs
            cluster_ids              - optionally restrict the counts to only
                                  those TF profiles identified by these
                                  TF IDs
			operon_gene_ids		- hash with key = operon gene id, val = first
								  operon gene id. used to replicate the counts
								  for the first gene with other genes belonging
								  to the same operon
			has_operon			- indicates whether this species has operons
								  optional operon_gene_ids are specified
=cut

sub fetch_custom_counts
{
    my ($self, %args) = @_;

    my $cons_level    = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};
    my $gids          = $args{-gene_ids};
    my $clusters      = $args{-clusters};
	my $operon_gids  = $args{-operon_gene_ids};
	my $has_operon   = $args{-has_operon};

    if (   !$cons_level
        || !defined $threshold
        || !defined $upstream_bp
        || !defined $downstream_bp)
    {
        carp "must provide conservation_level, threshold, upstream_bp"
            . " and downstream_bp\n";
        return;
    }

    #
    # If no gene IDs are passed in, get all the gene IDs from the database
    # in order to initialize counts object properly, i.e. ensure that
    # counts contain 0's for all gene/TF combinations even if there is
    # no record in the database cluster_counts table (0 values not stored
    # for efficiency).
    #
    unless ($gids) {
        my $ga = $self->db->get_GeneAdaptor;
        if (!$ga) {
            carp "error getting GeneAdaptor\n";
            return;
        }
		
        $gids = $ga->fetch_gene_ids();
		if ($has_operon) {
			my $oa = $self->db->get_OperonAdaptor();
			if (!$oa) {
				carp "error getting OperonAdaptor\n";
				return;
			}
			my $operons = $oa->fetch_where();
			foreach my $operon (@$operons) {
				my $fgene = $operon->fetch_first_gene();
				foreach my $opgene (@{$operon->genes}) {
					$$operon_gids{$opgene->id} = $fgene->id;
					#	if $opgene->id != $fgene->id;
				}
			}
		}
    }
	
    #
	# Set t_gids to all IDs of genes in the input gene IDs list which are
    # either not in an operon or are the first gene in an operon containing
    # an input gene.
    #
	my @t_gids;
    my %gid_inc;
	if ($operon_gids) {
		foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            if ($first_gid) {
                unless ($gid_inc{$first_gid}) {
                    push @t_gids, $first_gid;
                    $gid_inc{$first_gid} = 1;
                }
            } else {
                unless ($gid_inc{$gid}) {
                    push @t_gids, $gid;
                    $gid_inc{$gid} = 1;
                }
            }
		}
	} else {
		@t_gids = @$gids;
	}
	
	my @cids;
	foreach my $cluster (@$clusters) {
		push @cids, $cluster->id;
	}
	
    my $t_counts = OPOSSUM::Analysis::Cluster::Counts->new(
        -gene_ids       => \@t_gids,
		-cluster_ids	=> \@cids
    );
	
    my $ctfsa = $self->db->get_ConservedTFBSAdaptor();
    if (!$ctfsa) {
        carp "error getting ConservedTFBSAdaptor\n";
        return;
    }
	
	foreach my $cluster (@$clusters)
	{
		my $tfids = $cluster->tf_ids;
		
		foreach my $gid (@t_gids)
		{
			my $ctfsites = $ctfsa->fetch_by_gene(
				-gene_id			=> $gid,
				-conservation_level	=> $cons_level,
				-threshold			=> $threshold,
				-upstream_bp		=> $upstream_bp,
				-downstream_bp		=> $downstream_bp,
				-tf_ids				=> $tfids
			);
			
			next if !$ctfsites;
		
			my @sorted_sites = sort {$a->start <=> $b->start} @$ctfsites;
			my $merged_sites = merge_cluster_sites(\@sorted_sites, $cluster->id);

			my $count = $merged_sites ? scalar(@$merged_sites) : 0;
			my $length = 0;
			foreach my $site (@$merged_sites) {
				$length += length($site->seq);
			}
			
			$t_counts->gene_cluster_count($gid, $cluster->id, $count);
			$t_counts->gene_cluster_length($gid, $cluster->id, $length);
        }
    }

	#
	# Now, add the operon gene counts
	#
	my $counts;
	if ($operon_gids) {
		
		$counts = OPOSSUM::Analysis::Cluster::Counts->new(
			-gene_ids		=> $gids,
			-cluster_ids	=> \@cids
		);
		
		foreach my $gid (@$gids) {
			my $first_gid = $operon_gids->{$gid};
			unless ($first_gid) {
				$first_gid = $gid;
			}
			
			foreach my $cid (@cids) {
				my $count = $t_counts->gene_cluster_count($first_gid, $cid);
				my $length = $t_counts->gene_cluster_length($first_gid, $cid);
				$counts->gene_cluster_count($gid, $cid, $count);
				$counts->gene_cluster_length($gid, $cid, $length);
			}
		}
	} else {
		$counts = $t_counts;
	}
	
    $counts->param('conservation_level', $cons_level);
    $counts->param('threshold',          $threshold);
    $counts->param('upstream_bp',        $upstream_bp);
    $counts->param('downstream_bp',      $downstream_bp);

    return $counts;
}


=head2 fetch_anchored_counts

 Title    : fetch_anchored_counts

 Usage    : $tf_set = $ctfsa->fetch_anchored_counts(
                -gene_ids           => $gene_ids,
                -anchor_cluster     => $cluster,
                -cluster_ids        => $cluster_ids,
                -distance           => $dist,
                -conservation_level => $level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp,
                -has_operon         => $has_operon,
                -operon_gene_ids    => $operon_gids
            );

 Function : For the given genes, fetch the counts of conserved TFBS clusters
			which are proximal to (within the inter-binding distance of) the
			anchoring TFBS cluster.
 Returns  : A listref of OPOSSUM::TFBSCount objects.
 Args	  : Integer gene ID and integer conservation level.
			operon_gene_ids		- hash with key = operon gene id, val = first
								  operon gene id. used to replicate the counts
								  for the first gene with other genes belonging
								  to the same operon
			has_operon			- indicates whether this species has operons
								  optional operon_gene_ids are specified
 AK's Question: no all gene ids option for this method?
 
=cut

sub fetch_anchored_counts
{
    my ($self, %args) = @_;

    my $gids                = $args{-gene_ids};
    my $anchor_cluster      = $args{-anchor_cluster};
    my $cluster_set         = $args{-cluster_set};
    my $distance            = $args{-distance};
    my $conservation_level  = $args{-conservation_level};
    my $threshold           = $args{-threshold};
    my $upstream_bp         = $args{-upstream_bp};
    my $downstream_bp       = $args{-downstream_bp};
	my $operon_gids  		= $args{-operon_gene_ids};
	my $has_operon   		= $args{-has_operon};
	
    if (!$gids) {
        carp "must provide gene IDs";
        return;
    }

    if (!$anchor_cluster) {
        carp "must provide anchoring TFBS cluster";
        return;
    }

    if (!$cluster_set) {
        carp "must provide the TFClusterSet to be searched";
        return;
    }

    if (!$conservation_level) {
        carp "must provide a conservation level";
        return;
    }

    my $ctfbsa = $self->db->get_ConservedTFBSAdaptor();
    if (!$ctfbsa) {
        carp "error getting ConservedTFBSAdaptor\n";
        return;
    }
	
	my $cids = $cluster_set->cluster_ids;
	
    #
	# Set t_gids to all IDs of genes in the input gene IDs list which are
    # either not in an operon or are the first gene in an operon containing
    # an input gene.
    #
	my @t_gids;
    my %gid_inc;
	if ($operon_gids) {
		foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            if ($first_gid) {
                unless ($gid_inc{$first_gid}) {
                    push @t_gids, $first_gid;
                    $gid_inc{$first_gid} = 1;
                }
            } else {
                unless ($gid_inc{$gid}) {
                    push @t_gids, $gid;
                    $gid_inc{$gid} = 1;
                }
            }
		}
	} else {
		@t_gids = @$gids;
	}
	
	my $t_counts = OPOSSUM::Analysis::Cluster::Counts->new(
        -gene_ids   	=> \@t_gids,
        -cluster_ids    => $cids
    );
	
	my $anchor_tf_ids = $anchor_cluster->tf_ids;
    foreach my $gid (@t_gids) {
        my $anchor_tfbss = $ctfbsa->fetch_by_gene(
            -gene_id                => $gid,
            -tf_ids                 => $anchor_tf_ids,
            -conservation_level     => $conservation_level,
            -threshold              => $threshold,
            -upstream_bp            => $upstream_bp,
            -downstream_bp          => $downstream_bp
        );
		
        unless ($anchor_tfbss) {
            foreach my $cid (@$cids) {
                #next if $cid eq $anchor_cluster->id;
                $t_counts->gene_cluster_count($gid, $cid, 0);
				$t_counts->gene_cluster_length($gid, $cid, 0);
            }
            next;
        }
		
		# merge the anchor_tfbss into cluster sites
		my @sorted_anchor_sites = sort {$a->start <=> $b->start} @$anchor_tfbss;
		my $merged_anchor_sites = merge_cluster_sites(
            \@sorted_anchor_sites, $anchor_cluster->id
        );
        
		foreach my $cid (@$cids) {
			#next if $cid eq $anchor_cluster->id;

			my $cluster = $cluster_set->get_tf_cluster($cid);

            my $ctfsites = $ctfbsa->fetch_by_gene(
                -gene_id                => $gid,
                -tf_ids                 => $cluster->tf_ids,
                -conservation_level     => $conservation_level,
                -threshold              => $threshold,
                -upstream_bp            => $upstream_bp,
                -downstream_bp          => $downstream_bp
            );

			if (!$ctfsites) {
                $t_counts->gene_cluster_count($gid, $cid, 0);
                next;
            }
			
			my @sorted_ctfsites = sort {$a->start <=> $b->start} @$ctfsites;
			my $merged_ctfsites = merge_cluster_sites(\@sorted_ctfsites, $cid);
            
            my $prox_cluster_sites = _proximal_cluster_sites(
                $merged_anchor_sites, $merged_ctfsites, $distance
            );

			my $count = 0;
			$count = scalar @$prox_cluster_sites if $prox_cluster_sites;
			my $length = 0;
			foreach my $site (@$prox_cluster_sites) {
				$length += length($site->seq);
			}
			
            if ($prox_cluster_sites) {
                $t_counts->gene_cluster_count($gid, $cid, $count);
                $t_counts->gene_cluster_length($gid, $cid, $length)
				#$tf_gid_sites{$tf_id}->{$gid} = $prox_tfbss;
            } else {
                $t_counts->gene_cluster_count($gid, $cid, 0);
				$t_counts->gene_cluster_length($gid, $cid, 0);
            }
        }
    }
	
	#
	# Now, set final gene counts
	#
    my $counts;
	if ($operon_gids) {
        $counts = OPOSSUM::Analysis::Cluster::Counts->new(
            -cluster_ids    => $cids,
            -gene_ids   	=> $gids
        );

		foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            unless ($first_gid) {
                $first_gid = $gid;
            }

			foreach my $cid (@$cids) {
                #next if $cid eq $anchor_cluster->id;
				#my $cluster = $cluster_set->get_tf_cluster(id);

                my $count = $t_counts->gene_cluster_count($first_gid, $cid);
				my $length = $t_counts->gene_cluster_length($first_gid, $cid);
                $counts->gene_cluster_count($gid, $cid, $count);
				$counts->gene_cluster_length($gid, $cid, $length);
			}
		}
	} else {
        $counts = $t_counts;
    }
	
    $counts->param('conservation_level', $conservation_level);
    $counts->param('threshold',          $threshold);
    $counts->param('upstream_bp',        $upstream_bp);
    $counts->param('downstream_bp',      $downstream_bp);

    return $counts;
}

#
# Find all TFBS clusters proximal (within $max_dist) to the anchor TFBS cluster.
#
# TFBS clusters which overlap anchor TFBS clusters are excluded.
#
# NOTE: A given TFBS cluster could be proximal to more than one anchor. It is counted
# in combination with each anchor (multiple times).
#
sub _proximal_cluster_sites
{
    my ($anchor_csites, $csites, $max_dist) = @_;

    my @prox_csites;
    foreach my $anchor (@$anchor_csites) {
        foreach my $site (@$csites) {
            if ($site->id eq $anchor->id) {
                # If TF in question is same as anchor TF only count sites where
                # the TF is to the right of the anchor to avoid double counting
                next if $site->start <= $anchor->end;
            }

            my $dist;
            if ($site->start() > $anchor->end()) {
                $dist = $site->start() - $anchor->end() - 1;
            } elsif ($anchor->start() > $site->end()) {
                $dist = $anchor->start() - $site->end() - 1;
            } else {
                # do not include TFBS clusters which overlap anchor TFBS cluster
                next;
            }

            if ($dist <= $max_dist) {
                push @prox_csites, $site;
            }
        }
    }

    return @prox_csites ? \@prox_csites : undef;
}

sub merge_cluster_sites
{
    my ($tfsites, $cluster_id) = @_;
    
    if (!defined $tfsites) {
        return;
    }
    
    my @merged_sites;

    #
    # Let's not mess around setting and resetting strand info. Strand is
    # meaningless for clusters. Just keep everything on the +ve strand.
    # DJA 2011/05/31
    #
    my $tfsite1 = $tfsites->[0];
    if ($tfsite1->strand == -1) {
        $tfsite1->strand(1);
        $tfsite1->seq(revcom($tfsite1->seq));
    }

    push @merged_sites, $tfsite1;

    for (my $i = 1; $i < scalar(@$tfsites); $i++) {        
        my $tfsite = $tfsites->[$i];
        if ($tfsite->strand == -1) {
            $tfsite->strand(1);
            $tfsite->seq(revcom($tfsite->seq));
        }

        my $prevsite = $merged_sites[$#merged_sites];
        $prevsite->id($cluster_id);
        
        # if overlap, keep the max score
        # merge the two sites
        if (overlap($prevsite, $tfsite)) {
            if ($prevsite->end < $tfsite->end) {
                #$prevsite->end($tfsite->end);
                #
                # merge the sequences
                # first, check the strands of the sites
                # if negative, reverse complement
                # I should only do this if they are overlapping
                #if ($prevsite->strand != $tfsite->strand) {
                #    if ($prevsite->strand == -1) {
                #        my $seq = Bio::Seq->new(-seq => $prevsite->seq);
                #        $prevsite->seq($seq->revcom->seq);
                #    } else {
                #        my $seq = Bio::Seq->new(-seq => $tfsite->seq);
                #        $tfsite->seq($seq->revcom->seq);
                #    }
                #}
                
                # DJA 2011/05/31
                my $ext_seq = substr(
                    $tfsite->seq, $prevsite->end - $tfsite->start + 1
                );

                #
                # If ends of current tfbs and prev tfbs are the same, then
                # ext_seq is undef.
                # DJA 2011/05/31
                #
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
    
    return \@merged_sites;
}

sub overlap
{
    my ($tf1, $tf2) = @_;
    
    # DJA 2011/05/31
    #if (($tf1->start <= $tf2->start and $tf1->end > $tf2->start)
    #    or ($tf2->start <= $tf1->start and $tf2->end > $tf1->start))
    if ($tf1->start <= $tf2->end && $tf1->end >= $tf2->start) {
        return 1;
    }
    return 0;
}

sub revcom
{
    my ($seq) = @_;

    my $rc_seq = reverse $seq;

    $rc_seq =~ tr/[acgtACGT]/[tgcaTGCA]/;

    return $rc_seq;
}

1;
