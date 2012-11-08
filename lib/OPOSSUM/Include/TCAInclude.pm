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

    #
    # The code below will only work if the incoming TFBSs are already sorted
    # DJA 2012/11/01
    #
    @$ctfbss = sort {$a->start <=> $b->start} @$ctfbss;

	my $prev_site = $ctfbss->[0];
    $prev_site->id($cluster_id);
    
    #
    # Let's not mess around setting and resetting strand info. Strand is
    # meaningless for clusters. Just keep everything on the +ve strand.
    # DJA 2012/11/02
    #
    if ($prev_site->strand == -1) {
        $prev_site->strand(1);
        $prev_site->seq(revcom($prev_site->seq));
    }

    my @merged_sites;
	push @merged_sites, $prev_site;
    
	for (my $i = 0; $i < scalar @$ctfbss; $i++) {
		my $curr_site = $ctfbss->[$i];

        $curr_site->id($cluster_id);

        if ($curr_site->strand == -1) {
            $curr_site->strand(1);
            $curr_site->seq(revcom($curr_site->seq));
        }

        my $prev_site = $merged_sites[$#merged_sites];
        
        # if overlap (or adjacent - DJA), keep the max score
        # merge the two sites
		if (overlap($prev_site, $curr_site, 1)) {
			if ($prev_site->end < $curr_site->end) {
                # merge the sequences
				my $ext_seq = substr(
                    $curr_site->seq, $prev_site->end - $curr_site->start + 1
                );
				
                if ($ext_seq) {
                    $prev_site->seq($prev_site->seq . $ext_seq);
                }

				$prev_site->end($curr_site->end);
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

#
# Re-wrote this to simplify and also check for adjacent sites (we want to
# merge adjacent sites in clusters).
# DJA 2012/11/02
#
sub overlap
{
    my ($tf1, $tf2, $adjacent_ok) = @_;
    
    if ($adjacent_ok) {
        #
        # Sites can also be directly adjacent to one another
        #
        if ($tf1->start <= ($tf2->end + 1) && $tf1->end >= ($tf2->start - 1)) {
            return 1;
        }
    } else {
        #
        # The sites must overlap by at least 1 bp
        #
        if ($tf1->start <= $tf2->end && $tf1->end >= $tf2->start) {
            return 1;
        }
    }

    return 0;
}

1;
