=head1 NAME

 OPOSSUM::Include::TCAInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the TFBS Cluster
  type scripts and modules.

=head1 AUTHOR

  Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca

=cut

use strict;

use OPOSSUM::Opt::TCAOpt;
use OPOSSUM::Web::Opt::BaseOpt;

use lib TFBS_CLUSTER_LIB_PATH;

use TFBSCluster::DBSQL::DBAdaptor;


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

    return $cluster_set;
}


# tfsites belong to 1 cluster only
# returns listref of OPOSSUM::ConservedTFBS objects
sub merge_cluster_ctfbs_sites
{
    my ($ctfbss, $cluster_id) = @_;
    
    return if !$ctfbss or @$ctfbss == 0;
    
    #my $it = $tf_siteset->Iterator(-sort_by => 'start');    
    #my $prev_site = $it->next;
	my $prev_site = $$ctfbss[0];

    my @merged_sites;
	push @merged_sites, $prev_site;
    
    #while (my $curr_site = $it->next)
	for (my $i = 0; $i < scalar @$ctfbss; $i++)
    {
		my $curr_site = $$ctfbss[$i];
        my $prev_site = $merged_sites[$#merged_sites];
        $prev_site->id($cluster_id);
        
        # if overlap, keep the max score
        # merge the two sites
		if (overlap($prev_site, $curr_site))
        {
			if ($prev_site->end < $curr_site->end) {
				
				$prev_site->end($curr_site->end);
                
                # merge the sequences
                # first, check the strands of the sites
                # if negative, reverse complement
                # I should only do this if they are overlapping
				if ($prev_site->strand != $curr_site->strand) {
                    if ($prev_site->strand == -1) {
                        my $seq = Bio::Seq->new(-seq => $prev_site->seq);
                        $prev_site->seq($seq->revcom->seq);
                    } else {
                        my $seq = Bio::Seq->new(-seq => $curr_site->seq);
                        $curr_site->seq($seq->revcom->seq);
                    }
				}
                
				my $ext_seq = substr($curr_site->seq, $prev_site->end - $curr_site->start);
				$prev_site->seq($prev_site->seq . $ext_seq);
            }

			if ($curr_site->score > $prev_site->score) {
				$prev_site->score($curr_site->score);
			}
			if ($curr_site->rel_score > $prev_site->rel_score) {
				$prev_site->rel_score($curr_site->rel_score);
			}

        } else {
            $curr_site->id($cluster_id);
			push @merged_sites, $curr_site;
        }
    }
    
	return @merged_sites ? \@merged_sites : undef;
}


# tfsites belong to 1 cluster only
sub merge_cluster_tfbs_siteset
{
    my ($tfsiteset, $cluster_id, $seq_id) = @_;
    
    return if !$tfsiteset or $tfsiteset->size == 0;
	
	my @sites = $tfsiteset->all_sites;
	my $ctfbss = tfbss_to_conserved_tfbss(\@sites, $cluster_id, $seq_id);

	return merge_cluster_ctfbs_sites($ctfbss, $cluster_id);
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

1;
