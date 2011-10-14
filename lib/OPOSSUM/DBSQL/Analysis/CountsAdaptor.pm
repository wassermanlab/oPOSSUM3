
=head1 NAME

OPOSSUM::DBSQL::Analysis::CountsAdaptor - Adaptor for MySQL queries to retrieve
and store gene TFBS counts.

=head1 SYNOPSIS

$aca = $db_adaptor->get_AnalysisCountsAdaptor();

=head1 DESCRIPTION

In order to facilitate fast retrieval of TFBS counts from the oPOSSUM database
several count sets were pre-computed using discrete values for PWM (PSSM)
thresholds, conservation levels, and upstream/downstream search regions.
This adaptor provides an interface to both the pre-computed counts as well as
the raw data used to compute counts and retrieves the information into a
OPOSSUM::Analysis::Counts object which can be passed directly to the
OPOSSUM::Analysis::Fisher and WormOPOSSUM::Analysis::Zscore modules.

=head1 CHANGE HISTORY

 AK 010/10/28
 - added arguments -operon_gene_ids and -has_operon to all variations of
   fetch_counts(), fetch_custom_counts(), and fetch_anchored_counts()
 - Now, if input gene ids contain operon genes, the counts from the first
   gene in the operon are duplicated for all other genes belonging to the
   same operon
   
 DJA 2010/10/26
 - added method fetch_anchored_counts

 DJA 2007/01/29
 - modified fetch_custom_counts method to call
   the new ConservedRegionAdaptor->fetch_length_by_upstream_downstream
   method to get the length directly instead of getting all the conserved
   region records with fetch_set_by_upstream_downstream and then calling
   the length method on them (more efficient???)

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

 Modified by Shannan Ho Sui on Dec 21, 2006 to accommodate schema changes

=head1 METHODS

=cut

package OPOSSUM::DBSQL::Analysis::CountsAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Analysis::Counts;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;

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
                -tf_ids                 => $tfids,
				-operon_gene_ids		=> $operon_gids,
				-has_operon				=> $has_operon
            );
 Function : Fetch the pre-computed counts of TFBS hits on the gene
            sequences from the DB and build a
            OPOSSUM::Analysis::Counts object.
 Returns  : An OPOSSUM::Analysis::Counts object.
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
            tf_ids		        - optionally restrict the counts to only
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
    my $tfids        = $args{-tf_ids};
	my $operon_gids  = $args{-operon_gene_ids};
	my $has_operon   = $args{-has_operon};

    if (!$cons_level || !$thresh_level || !$sr_level) {
        carp "must provide conservation_level, threshold_level and"
            . " search_region_level\n";
        return;
    }

    #
    # Keep track of whether gene IDs were passed in initially, i.e. if
    # none were then we are fetching counts for all genes in the oPOSSUM DB
    # (full background set)
    # DJA 2011/04/26
    #
    my $all_genes = 0;
    unless ($gids && $gids->[0]) {
        $all_genes = 1;
    }
	
    #my $crla = $self->db->get_ConservedRegionLengthAdaptor();
    #if (!$crla) {
    #    carp "error getting ConservedRegionLengthAdaptor\n";
    #    return;
    #}

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

    my $sql = qq{
		select gene_id, tf_id, count
		from tfbs_counts
		where conservation_level = $cons_level
		and threshold_level = $thresh_level
		and search_region_level = $sr_level
    };
        
    #
    # For efficiency, if we are fetching counts for all genes in the DB and
    # we are not dealing with operons, then don't add the "in" clause for gene
    # IDs to the SQL query.
    # DJA 2011/04/26
    #
    unless ($all_genes && !$has_operon) {
        $sql .= " and gene_id in (";
        $sql .= join ',', @t_gids;
        $sql .= ")";
    }
	
    if ($tfids && $tfids->[0]) {
        $sql .= " and tf_id in ('";
        $sql .= join "','", @$tfids;
        $sql .= "')";
    }

    #print STDERR "\n\nsql:\n\t$sql\n\n";

    #printf STDERR "%s: preparing fetch TFBS counts\n", scalar localtime;

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch TFBS counts with\n$sql\n" . $self->errstr;
        return;
    }

    #printf STDERR "%s: executing fetch TFBS counts\n", scalar localtime;

    if (!$sth->execute) {
        carp "Error executing fetch TFBS counts with\n$sql\n" . $sth->errstr;
        return;
    }
	
    #printf STDERR "%s: fetching TFBS counts\n", scalar localtime;

    my $raw_counts = $sth->fetchall_arrayref();
    if ($sth->err()) {
        carp "Error fetching gene TFBS counts: " . $sth->errstr . "\n";
        return;
    }

    #printf STDERR "%s: setting analysis counts object\n", scalar localtime;
	
    my $t_counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids   => \@t_gids,
        -tf_ids     => $tfids,
        -counts     => $raw_counts
    );

	#
	# Now, add the operon gene counts
	#
    my $counts;
	if ($operon_gids) {
        $counts = OPOSSUM::Analysis::Counts->new(
            -gene_ids   => $gids,
            -tf_ids     => $tfids
        );

        foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            unless ($first_gid) {
                $first_gid = $gid;
            }

            foreach my $tfid (@$tfids) {
                my $count = $t_counts->gene_tfbs_count($first_gid, $tfid);
                $counts->gene_tfbs_count($gid, $tfid, $count);
            }
        }
    } else {
        $counts = $t_counts;
	}

    #if ($gids) {
    #    $counts->gene_ids($gids);
    #}

    #if ($tfids) {
    #    $counts->tf_ids($tfids);
    #}

    #my $crl_set = $crla->fetch_length_set(
    #    -gene_ids            => $gids,
    #    -conservation_level  => $cons_level,
    #    -search_region_level => $sr_level
    #);
    #
    #$counts->conserved_region_length_set($crl_set);

    $counts->param('conservation_level',  $cons_level);
    $counts->param('threshold_level',     $thresh_level);
    $counts->param('search_region_level', $sr_level);

    #printf STDERR "%s: returning analysis counts object\n", scalar localtime;

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
                -tf_ids                 => $tfids,
				-operon_gene_ids		=> $op_gids,
				-has_operon				=> $has_operon);
 Function : Compute the counts of TFBS hits on the gene sequences
            from the DB and build an OPOSSUM::Analysis::Counts object.
 Returns  : An OPOSSUM::Analysis::Counts object.
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
            tf_ids              - optionally restrict the counts to only
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
    my $tfids         = $args{-tf_ids};
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
    # no record in the database tfbs_counts table (0 values not stored
    # for efficiency).
    #
    #
	# AK: This one would also need operon treatment.
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
	
    my $t_counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids   => \@t_gids,
        -tf_ids     => $tfids
    );

    #my $cra = $self->db->get_ConservedRegionAdaptor;
    #if (!$cra) {
    #    carp "error getting ConservedRegionAdaptor\n";
    #    return;
    #}

    my $ctfsa = $self->db->get_ConservedTFBSAdaptor();
    if (!$ctfsa) {
        carp "error getting ConservedTFBSAdaptor\n";
        return;
    }

    #my $crl_set = OPOSSUM::ConservedRegionLengthSet->new();
    my $count   = 0;
    
	foreach my $gid (@t_gids) {
        #print STDERR "$count\tGene: $gid\n";
        #
        # DJA 2007/01/29
        # Call the new fetch_length_by_upstream_downstream method from
        # ConservedRegionAdaptor to get the length directly instead of
        # getting all the conserved region records with then calling
        # the length method on them (more efficient???)
        #
        #my $cr_set = $cra->fetch_set_by_upstream_downstream(
        #		    $gid, $cons_level, $upstream_bp, $downstream_bp);
        #my $length =
        #    $cra->fetch_length_by_upstream_downstream($gid, $cons_level,
        #    $upstream_bp, $downstream_bp);
        #if (defined($cr_set)) {
        #if ($length) {
            #my $crl = OPOSSUM::ConservedRegionLength->new(
            #    -gene_id            => $gid,
            #    -conservation_level => $cons_level,
            #    -length             => $length
            #);
            #$crl_set->add_conserved_region_length($crl);

            my $gp_counts = $ctfsa->fetch_gene_tfbs_counts(
                -gene_id            => $gid,
                -conservation_level => $cons_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp,
                -tf_ids             => $tfids
            );

            foreach my $gp_count (@$gp_counts) {
                $t_counts->gene_tfbs_count(
                    $gp_count->gene_id,
                    $gp_count->tf_id,
                    $gp_count->count
                );
            }
        #}
        $count++;
    }

	#
	# Now, add the operon gene counts
	#
    my $counts;
	if ($operon_gids) {
        $counts = OPOSSUM::Analysis::Counts->new(
            -gene_ids   => $gids,
            -tf_ids     => $tfids
        );

        foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            unless ($first_gid) {
                $first_gid = $gid;
            }

            foreach my $tfid (@$tfids) {
                my $count = $t_counts->gene_tfbs_count($first_gid, $tfid);
                $counts->gene_tfbs_count($gid, $tfid, $count);
            }
        }
    } else {
        $counts = $t_counts;
	}
	
    #$counts->conserved_region_length_set($crl_set);

    #$counts->gene_ids($gids);

    #if ($tfids) {
    #    $counts->tf_ids($tfids);
    #}

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
                -anchor_tf_id       => $tf_id,
                -tf_ids             => $tf_ids,
                -distance           => $dist,
                -conservation_level => $level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp,
                -has_operon         => $has_operon,
                -operon_gene_ids    => $operon_gids
            );

 Function : For the given genes, fetch the counts of conserved TF binding
            sites which are proximal to (within the inter-binding distance
            of) the anchoring TF.
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
    my $anchor_tf_id        = $args{-anchor_tf_id};
    my $tf_ids              = $args{-tf_ids};
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

    if (!$anchor_tf_id) {
        carp "must provide anchoring TF ID";
        return;
    }

    if (!$tf_ids) {
        carp "must provide a TF IDs";
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
    
	my $t_counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids   => \@t_gids,
        -tf_ids     => $tf_ids
    );
	
    foreach my $gid (@t_gids) {
        my $anchor_tfbss = $ctfbsa->fetch_by_gene(
            -gene_id                => $gid,
            -tf_id                  => $anchor_tf_id,
            -conservation_level     => $conservation_level,
            -threshold              => $threshold,
            -upstream_bp            => $upstream_bp,
            -downstream_bp          => $downstream_bp
        );

        unless ($anchor_tfbss) {
            foreach my $tf_id (@$tf_ids) {
                #next if $tf_id eq $anchor_tf_id;
                $t_counts->gene_tfbs_count($gid, $tf_id, 0);
            }
            next;
        }

        foreach my $tf_id (@$tf_ids) {
            #next if $tf_id eq $anchor_tf_id;

            my $tfbss = $ctfbsa->fetch_by_gene(
                -gene_id                => $gid,
                -tf_id                  => $tf_id,
                -conservation_level     => $conservation_level,
                -threshold              => $threshold,
                -upstream_bp            => $upstream_bp,
                -downstream_bp          => $downstream_bp
            );

            if (!$tfbss) {
                $t_counts->gene_tfbs_count($gid, $tf_id, 0);
                next;
            }

            my $prox_tfbss = _proximal_tfbss(
                $anchor_tfbss, $tfbss, $distance
            );

            if ($prox_tfbss) {
                $t_counts->gene_tfbs_count($gid, $tf_id, scalar @$prox_tfbss);
                #$tf_gid_sites{$tf_id}->{$gid} = $prox_tfbss;
            } else {
                $t_counts->gene_tfbs_count($gid, $tf_id, 0);
            }
        }
    }
	
	#
	# Now, set final gene counts
	#
    my $counts;
	if ($operon_gids) {
        $counts = OPOSSUM::Analysis::Counts->new(
            -tf_ids     => $tf_ids,
            -gene_ids   => $gids
        );

		foreach my $gid (@$gids) {
            my $first_gid = $operon_gids->{$gid};
            unless ($first_gid) {
                $first_gid = $gid;
            }

			foreach my $tfid (@$tf_ids) {
                next if $tfid eq $anchor_tf_id;

                my $count = $t_counts->gene_tfbs_count($first_gid, $tfid);
                $counts->gene_tfbs_count($gid, $tfid, $count);
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
# Find all TFBSs proximal (within $max_dist) to the anchor TFBSs.
#
# TFBSs which overlap anchor TFBSs are excluded.
#
# NOTE: A given TFBS could be proximal to more than one anchor. It is counted
# in combination with each anchor (multiple times).
#
sub _proximal_tfbss
{
    my ($anchor_tfbss, $tfbss, $max_dist) = @_;

    my @prox_tfbss;
    foreach my $anchor (@$anchor_tfbss) {
        foreach my $tfbs (@$tfbss) {
            if ($tfbs->id eq $anchor->id) {
                # If TF in question is same as anchor TF only count sites where
                # the TF is to the right of the anchor to avoid double counting
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
                push @prox_tfbss, $tfbs;
            }
        }
    }

    return @prox_tfbss ? \@prox_tfbss : undef;
}

1;
