#package oPossumTCAWebInclude;


# This module should be included in all the oPossum*TCAWeb.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to all the oPossum TCA variants.
#

use oPossumWebOpt;
use oPossumTCAWebOpt;

use lib TFBS_CLUSTER_LIB_PATH;

use Data::Dumper;    # for debugging only

#use OPOSSUM::Analysis::Cluster::Counts;


use TFBSCluster::DBSQL::DBAdaptor;
use TFBSCluster::TFCluster;
use TFBSCluster::TFClusterSet;
use TFBS::DB::JASPAR5;

#use Template;
#use CGI::Carp qw(carpout);    # fatalsToBrowser;

use strict;


sub cldba
{
    my $self = shift;

    if (@_) {
        $self->{-cldba} = shift;
    }

    return $self->{-cldba};
}

sub write_cluster_results
{
    my ($self, $filename, $results, $cluster_set) = @_;

    my $text = "TFBS Cluster Name\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget Cluster hits\tTarget Cluster nucleotide rate\tBackground Cluster hits\tBackground Cluster nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $cl = $cluster_set->get_tf_cluster($result->id());

        $text .= sprintf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $cl->name,
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_cluster_hits() || 0,
            defined $result->t_cluster_rate()
                ? sprintf("%.3f", $result->t_cluster_rate()) : 'N/A',
            $result->bg_cluster_hits() || 0,
            defined $result->bg_cluster_rate()
                ? sprintf("%.3f", $result->bg_cluster_rate()) : 'N/A',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'N/A',
            defined $result->fisher_p_value()
                ? sprintf("%.3g", $result->fisher_p_value()) : 'N/A';

    }

    unless (open(FH, ">$filename")) {
        $self->_error("Unable to create results text file $filename - $!");
        return;
    }

    print FH $text;
    close(FH);
}


sub tfbs_cluster_info
{
    my $self = shift;

    #print STDERR "input\n";
    
    my $q = $self->query;
    #print STDERR "input query:\n"
    #		    . Data::Dumper::Dumper($q);

    my $state = $self->state();
	#my $species = $state->species;
	my $cluster_id = $q->param('cluster_id');

    #
    # Connect to TFBSCluster DB and retrieve TFCluster info
    # Connect to JASPAR DB
	#
	$self->jaspar_db_connect;
    my $jdbh = $self->jdbh();
	$self->opossum_cluster_db_connect;
	my $cldba = $self->cldba();
	my $tfca = $cldba->get_TFClusterAdaptor;
    my $tfc = $tfca->fetch_by_cluster_id($cluster_id);
	
	my @tfc_tfs;
	my %collections;
	my %tax_groups;
	my %tf_ic;
	foreach my $tfid (@{$tfc->tf_ids})
	{
		my $tf = $jdbh->get_Matrix_by_ID($tfid);
		push @tfc_tfs, $tf;
		$collections{$tf->ID} = $tf->tag('collection');
		$tax_groups{$tf->ID} = $tf->tag('tax_group');
		$tf_ic{$tf->ID} = sprintf "%.2f", $tf->to_ICM->total_ic;
		#print STDERR $tf->ID . "\t" . $tf->tag('collection') . "\t";
		#print STDERR $tf->tag('tax_group') . "\n";
		
	}
	
    my $vars = {
        abs_htdocs_path    => ABS_HTDOCS_PATH,
        rel_htdocs_path    => REL_HTDOCS_PATH,
        abs_cgi_bin_path   => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path   => REL_CGI_BIN_PATH,
        bg_color_class     => $state->bg_color_class(),
        title              => $state->title(),
        heading            => $state->heading(),
        section            => 'TFBS Cluster Information',
        version            => VERSION,
        devel_version      => DEVEL_VERSION,
		jaspar_url         => JASPAR_URL,
        #sid                => $state->sid(),
        #species            => $species,
        #db_info            => $state->db_info(),
		tf_cluster         => $tfc,
		cluster_tfs        => \@tfc_tfs,
		collections        => \%collections,
		tax_groups         => \%tax_groups,
		tf_ic              => \%tf_ic,
        var_template       => "tfbs_cluster_info.html"
    };

    my $output = $self->process_template('master.html', $vars);
    #print STDERR "input results:\n"
    #		    . Data::Dumper::Dumper($output);

    return $output;
}

sub opossum_cluster_db_connect
{
    my $self = shift;

    my $dbh = TFBSCluster::DBSQL::DBAdaptor->new(
        -host     => TFBS_CLUSTER_DB_HOST,
        -dbname   => TFBS_CLUSTER_DB_NAME,
        -user     => TFBS_CLUSTER_DB_USER,
        -password => TFBS_CLUSTER_DB_PASS
    );

    if (!$dbh) {
        $self->_error("Could not connect to oPOSSUM_cluster database "
                      . TFBS_CLUSTER_DB_NAME);
    }

    $self->cldba($dbh);
}

sub fetch_cluster_analysis_counts
{
    my ($self, $analysis_type, %args) = @_;

    #printf STDERR "fetch_cluster_analysis_counts ($analysis_type) args:\n"
    #    . Data::Dumper::Dumper(%args) . "\n\n";

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $aca = $opdba->get_AnalysisClusterCountsAdaptor();
    if (!$aca) {
        $self->_error("Could not get AnalysisClusterCountsAdaptor");
    }

    my $counts;
    if ($analysis_type eq 'default') {
        printf STDERR "\nfetch_cluster_analysis_counts: fetching pre-computed counts\n";
        $counts = $aca->fetch_counts(%args);
    } elsif ($analysis_type eq 'custom') {
        printf STDERR "\nfetch_cluster_analysis_counts: fetching custom counts\n";
        $counts = $aca->fetch_custom_counts(%args);
    } else {
        $self->_error("fetch_cluster_analysis_counts: unknown analysis type");
    }

    if (!$counts) {
        $self->_error("Could not fetch analysis counts");
    }

    #$self->state->search_region_level_hash($srlh);

    return $counts;
}

sub fetch_tf_cluster_set
{
    my ($self, %args) = @_;

    my %matrix_args = %args;

    unless ($matrix_args{-matrixtype}) {
        $matrix_args{-matrixtype} = 'PFM';
    }

    #printf STDERR "fetch_tf_cluster_set: matrix_args = \n"
    #    . Data::Dumper::Dumper(%matrix_args) . "\n";

    my $cldba = $self->cldba();
    unless ($cldba) {
        $cldba = $self->opossum_cluster_db_connect();
    }
    
    my $tfca = $cldba->get_TFClusterAdaptor;

    my $cluster_set = TFBSCluster::TFClusterSet->new();

    my $clusters;
    #print STDERR "Getting cluster set\n";
    if ($matrix_args{-families}) {
        $clusters = $tfca->fetch_by_tf_families($matrix_args{-families});
    } elsif ($matrix_args{-family}) {
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

# tfsites belong to 1 cluster only
sub merge_cluster_sites
{
    my ($tfsites, $cluster_id) = @_;
    
    if (!defined $tfsites or scalar(@$tfsites) == 0) {
        print STDERR "merge_cluster_sites: No tf sites provided\n";
        return;
    }
    
    my @merged_sites;
    push @merged_sites, $$tfsites[0];
    for (my $i = 1; $i < scalar(@$tfsites); $i++)
    {        
        my $tfsite = $$tfsites[$i];
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
                #if ($i < 20) {
                #    print STDERR "prev: " . $prevsite->seq . "\tcurrent: " . $tfsite->seq;
                #    print STDERR "\text: $ext_seq\n";
                #    print STDERR $prevsite->seq . $ext_seq . "\n";
                #}                
				
				$prevsite->end($tfsite->end);                
                $prevsite->seq($prevsite->seq . $ext_seq);
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
    
    if (($tf1->start <= $tf2->start and $tf1->end > $tf2->start)
        or ($tf2->start <= $tf1->start and $tf2->end > $tf1->start))
    {
        return 1;
    }
    return 0;
}

1;