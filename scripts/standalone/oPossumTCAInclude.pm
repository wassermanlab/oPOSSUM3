
=head1 NAME

 oPossumTCAInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the TCA
  scripts.

=head1 AUTHOR

  Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca

=cut

use strict;

use oPossumOpt;
use oPossumTCAOpt;

use lib TFBS_CLUSTER_LIB_PATH;

sub tfbs_cluster_db_connect
{
	my ($cl_db) = @_;
	
    my $dbh = TFBSCluster::DBSQL::DBAdaptor->new(
        -host     => TFBS_CLUSTER_DB_HOST,
        -dbname   => TFBS_CLUSTER_DB_NAME,
        -user     => TFBS_CLUSTER_DB_USER,
        -password => TFBS_CLUSTER_DB_PASS
    );

    return $dbh;
}

# for now, the only argument used is -families
# I should add -min_ic and -tax_groups later
sub fetch_tf_cluster_set
{
    my ($cdb, %args) = @_;

    my %matrix_args = %args;

    unless ($matrix_args{-matrixtype}) {
        $matrix_args{-matrixtype} = 'PFM';
    }
    
    my $tfca = $cdb->get_TFClusterAdaptor;

    my $cluster_set = TFBSCluster::TFClusterSet->new();

    my $clusters;
    #print STDERR "Getting cluster set\n";
    if ($matrix_args{-family}) {
        $clusters = $tfca->fetch_by_tf_families($matrix_args{-family});
    } else {
        $clusters = $tfca->fetch_all();    
    }
    #print STDERR "# clusters = " . scalar(@$clusters) . "\n";
    $cluster_set->add_tf_cluster_list($clusters);

    die "Could not fetch TFBS clusters\n"
        if !$cluster_set || $cluster_set->size == 0;

    return $cluster_set;
}

# this needs to be fixed up to take tf_seq_siteset as input
# input: TFBS::SiteSet
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
    my $prev_ctfbs = OPOSSUM::ConservedTFBS->new(
        -tf_id => $cluster_id,
        -gene_id => $seq_id,
        -start => $tfsite->start,
        -end => $tfsite->end,
        -strand => $tfsite->strand,
        -score => $tfsite->score,
        -rel_score => $tfsite->rel_score,
        -seq => $tfsite->seq->seq
    );
    
    my @merged_sites;
    push @merged_sites, $prev_ctfbs;
    
    while ($tfsite = $it->next)
    {
        my $ctfbs = OPOSSUM::ConservedTFBS->new(
            -tf_id => $cluster_id,
            -gene_id => $seq_id,
            -start => $tfsite->start,
            -end => $tfsite->end,
            -strand => $tfsite->strand,
            -score => $tfsite->score,
            -rel_score => $tfsite->rel_score,
            -seq => $tfsite->seq->seq
        );
        
        $prev_ctfbs = $merged_sites[$#merged_sites];
        
        # if overlap, keep the max score
        # merge the two sites
        if (overlap($prev_ctfbs, $ctfbs))
        {
            if ($prev_ctfbs->end < $ctfbs->end) {
                $prev_ctfbs->end($ctfbs->end);
                
                # merge the sequences
                # first, check the strands of the sites
                # if negative, reverse complement
                # I should only do this if they are overlapping
                if ($prev_ctfbs->strand != $ctfbs->strand) {
                    if ($prev_ctfbs->strand == -1) {
                        my $seq = Bio::Seq->new(-seq => $prev_ctfbs->seq);
                        $prev_ctfbs->seq($seq->revcom->seq);
                    } else {
                        my $seq = Bio::Seq->new(-seq => $ctfbs->seq);
                        $ctfbs->seq($seq->revcom->seq);
                    }
                }
                
                my $ext_seq = substr($ctfbs->seq, $prev_ctfbs->end - $ctfbs->start);
                $prev_ctfbs->seq($prev_ctfbs->seq . $ext_seq);
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

# tfsites belong to 1 cluster only
sub merge_cluster_sites
{
    my ($tfsites, $cluster_id) = @_;
    
    return if !$tfsites;
    return if scalar(@$tfsites) == 0;

    my @sorted_sites = sort {$a->start <=> $b->start} @$tfsites;
    
    my @merged_sites;
    push @merged_sites, $sorted_sites[0];
    for (my $i = 1; $i < scalar(@sorted_sites); $i++)
    {        
        my $tfsite = $sorted_sites[$i];
        my $prevsite = $merged_sites[$#merged_sites];
        $prevsite->id($cluster_id);
        
        # if overlap, keep the max score
        # merge the two sites
        if (overlap($prevsite, $tfsite))
        {
            if ($prevsite->end < $tfsite->end) {
                
                # merge the sequences
                # first, check the strands of the sites
                # if negative, reverse complement
                # I should only do this if they are overlapping
				
                if ($prevsite->strand != $tfsite->strand) {
                    if ($prevsite->strand == -1) {
                        my $seq = Bio::Seq->new(-seq => $prevsite->seq);
                        $prevsite->seq($seq->revcom->seq);
                    } else {
                        my $seq = Bio::Seq->new(-seq => $tfsite->seq);
                        $tfsite->seq($seq->revcom->seq);
                    }
				}
				
				my $ext_seq = substr($tfsite->seq, $prevsite->end - $tfsite->start + 1);             
                $prevsite->seq($prevsite->seq . $ext_seq);
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

sub merged_proximal_tfbss
{
    my ($anchor_tfbss, $tfbss, $cluster_id, $max_dist) = @_;
    
    my $merged_anchor_tfbss = merge_cluster_sites($anchor_tfbss, $cluster_id);
    my $merged_tfbss = merge_cluster_sites($tfbss, $cluster_id);
    
    my @sitepairs;
    foreach my $anchor (@$merged_anchor_tfbss) {
        foreach my $tfbs (@$merged_tfbss) {
            my $dist;
            #if (!defined $anchor or !defined $tfbs) {
            #    $logger->error("Cluster ID: $cluster_id");
            #    $logger->error("anchor:\n" . Data::Dumper::Dumper($anchor));
            #    $logger->error("tfbs:\n" . Data::Dumper::Dumper($tfbs));
            #}
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
                    anchor_site     => $anchor,
                    cluster_site    => $tfbs,
                    distance        => $dist
                };
            }
        }
    }

    return @sitepairs ? \@sitepairs : undef;
}

sub overlap
{
    my ($tf1, $tf2) = @_;
    
    if (($tf1->start <= $tf2->start and $tf1->end > $tf2->start)
        or ($tf2->start <= $tf1->start and $tf2->end > $tf1->start))
    {
        return 1;
    }
    return 0;
}

#
# This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SitePairSet where each site pair in the set corresponds to the
# same transcription factor and filter overlapping site pairs such that only
# the highest scoring site pair of any mutually overlapping site pairs is kept.
# In the event that site pairs score equally, the first site pair is kept, i.e.
# bias is towards the site pair with the lowest starting position.
#
sub filter_overlapping_sites
{
    my ($siteset) = @_;

    return if !defined $siteset || $siteset->size == 0;

    my $filtered_set = TFBS::SiteSet->new();

    my $iter = $siteset->Iterator(-sort_by => 'start');
    my $prev_site = $iter->next;
    if ($prev_site) {
        while (my $site = $iter->next) {
            if ($site->overlaps($prev_site)) {
                #
                # Bias is toward the site pair with the lower start
                # site (i.e. if the scores are equal).
                # 
                if ($site->score > $prev_site->score) {
                    $prev_site = $site;
                }
            } else {
                $filtered_set->add_site($prev_site);
                $prev_site = $site;
            }
        }
        $filtered_set->add_site($prev_site);
    }

    return $filtered_set;
}

1;
