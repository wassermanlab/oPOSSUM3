=head1 NAME

OPOSSUM::Analysis::Cluster::Counts - Object to store the count of cluster cluster
hits on the gene sequences

=head1 SYNOPSIS

 my $aca = $db_adaptor->get_AnalysisClusterCountsAdaptor();

 my $counts = $aca->fetch_counts(
     -conservation_level     => 2,
     -threshold_level        => 3,
     -search_region_level    => 1
 );

=head1 DESCRIPTION

This object stores a count of the number of times each Cluster cluster was
found on each gene. These counts can be retrieved from the database
by the OPOSSUM::DBSQL::Analysis::Cluster::CountsAdaptor. This object can be
passed to the OPOSSUM::Analysis::Cluster::Fisher and
OPOSSUM::Analysis::Cluster::Zscore modules.

Note: for cluster analysis, simply counting hits is not enough.
One must make sure that overlapping cluster belonging to the same cluster should not be
read more than once - so have a separate tfbs_cluster_counts table

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

 Modified by Shannan Ho Sui on Dec 21, 2006 to accommodate schema changes

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::Counts;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $counts = OPOSSUM::Analysis::Cluster::Counts->new();
 Function : Create a new OPOSSUM::Analysis::Cluster::Counts object.
 Returns  : An OPOSSUM::Analysis::Cluster::Counts object.

=cut

sub new
{
    my ($class, %args) = @_;

    if ($args{-conserved_region_length_set}) {
        carp "Support for ConservedRegionLengthSet is deprecated.\n";
    }

    my $gene_ids = $args{-gene_ids};
    my $cluster_ids   = $args{-cluster_ids};

    my $self = bless {
        -gene_ids           => undef,
        -cluster_ids             => undef,
        -gene_exists        => undef,
        -cluster_exists          => undef,
        _gene_cluster_counts   => {},
		_gene_cluster_lengths	=> {},
        _cluster_gene_exists   => {},
        _params             => {}
    }, ref $class || $class;

    #
    # If gene ID and Cluster ID list provided, initialize counts with 0's
    #
    if ($gene_ids && $cluster_ids) {
        foreach my $gene_id (@$gene_ids) {
            foreach my $cluster_id (@$cluster_ids) {
                $self->gene_cluster_count($gene_id, $cluster_id, 0);
				$self->gene_cluster_length($gene_id, $cluster_id, 0);
            }
        }
    }

    return $self;
}

=head2 param

 Title    : param
 Usage    : $val = $counts->param($param)
	        or $counts->param($param, $value);
 Function : Get/set a value of a counts parameter.
 Returns  : The value of the names parameter.
 Args     : The name of the parameter to get/set.
            optionally the value of the parameter.

=cut

sub param
{
    my ($self, $param, $value) = @_;

    if ($param) {
        if (defined $value) {
            $self->{_params}->{$param} = $value;
        }
        return $self->{_params}->{$param};
    }
    return keys %{$self->{_params}};
}

=head2 gene_ids

 Title    : gene_ids
 Usage    : $gids = $counts->gene_ids()
 Function : Get the list of gene IDs stored in the counts object.
 Returns  : A reference to a list of gene IDs.
 Args     : None.

=cut

sub gene_ids
{
    my $self = shift;

    return $self->{-gene_ids};
}

=head2 get_all_gene_ids

 Title    : get_all_gene_ids
 Usage    : $gids = $counts->get_all_gene_ids()
 Function : Get the list of gene IDs. This is a synonym for the
            get variant of the gene_ids method.
 Returns  : A reference to a list of gene IDs.
 Args     : None.

=cut

sub get_all_gene_ids
{
    my $self = shift;

    return $self->gene_ids();
}

=head2 cluster_ids

 Title    : cluster_ids
 Usage    : $cids = $counts->cluster_ids()
 Function : Get the list of Cluster IDs stored in the counts object.
 Returns  : A reference to a list of Cluster IDs.
 Args     : None.

=cut

sub cluster_ids
{
    my $self = shift;

    return $self->{-cluster_ids}
}

=head2 get_all_cluster_ids

 Title    : get_all_cluster_ids
 Usage    : $tfids = $counts->get_all_cluster_ids()
 Function : Get the list of Cluster IDs. This is a synonym for the get
            variant of the cluster_ids method.
 Returns  : A reference to a list of Cluster IDs.
 Args     : Optionally a reference to a list of Cluster IDs.

=cut

sub get_all_cluster_ids
{
    my $self = shift;

    return $self->cluster_ids();
}

=head2 conserved_region_length_set

 Title    : conserved_region_length_set
 Usage    : This method is deprecated.

=cut

sub conserved_region_length_set
{
    carp "conserved_region_length_set() is deprecated\n";
}

=head2 get_conserved_region_length

 Title    : get_conserved_region_length
 Usage    : This method is deprecated.

=cut

sub get_conserved_region_length
{
    carp "get_conserved_region_length() is deprecated\n";
}

=head2 num_genes

 Title    : num_genes
 Usage    : $num = $counts->num_genes()
 Function : Get the number of genes/promoters in the counts object
 Returns  : An integer.
 Args     : None.

=cut

sub num_genes
{
    my $self = shift;

    my $num_genes = 0;
    if ($self->gene_ids()) {
        $num_genes = scalar @{$self->gene_ids()};
    }

    return $num_genes;
}

=head2 num_clusters

 Title    : num_clusters
 Usage    : $num = $counts->num_clusters()
 Function : Get the number of Clusters in the counts object
 Returns  : An integer.
 Args     : None.

=cut

sub num_clusters
{
    my $self = shift;

    my $num_tfs = 0;
    if ($self->cluster_ids()) {
        $num_tfs = scalar @{$self->cluster_ids()};
    }

    return $num_tfs;
}

=head2 gene_exists

 Title    : gene_exists
 Usage    : $bool = $counts->gene_exists($id)
 Function : Return whether the gene with the given ID exists in the
            counts object.
 Returns  : Boolean.
 Args     : Gene/promoter ID.

=cut

sub gene_exists
{
    my ($self, $id) = @_;

    return $self->{-gene_exists}->{$id};
}

=head2 cluster_exists

 Title    : cluster_exists
 Usage    : $bool = $counts->cluster_exists($id)
 Function : Return whether sites for the Cluster with the given ID exists in
            the counts object.
 Returns  : Boolean.
 Args     : Cluster ID.

=cut

sub cluster_exists
{
    my ($self, $id) = @_;

    return $self->{-cluster_exists}->{$id};
}

=head2 exists

 Title    : exists
 Usage    : $bool = $counts->exists($gene_id, $cluster_id)
 Function : Return whether the gene/Cluster pair with the given
            IDs exist in the counts object.
 Returns  : Boolean.
 Args     : A gene ID and a Cluster ID.

=cut

sub exists
{
    my ($self, $gene_id, $cluster_id) = @_;

    return $self->gene_exists($gene_id) && $self->cluster_exists($cluster_id);
}

=head2 total_cr_length

 Title    : total_cr_length
 Usage    : This method is deprecated

=cut

sub total_cr_length
{
    carp "total_cr_length() is deprecated\n";
}

=head2 gene_cluster_count

 Title    : gene_cluster_count
 Usage    : $count = $counts->gene_cluster_count($gene_id, $cluster_id);
            $counts->gene_cluster_count($gene_id, $cluster_id, $count);
 Function : Get/set the count of the number of times sites for the given
            Cluster were detected for the given gene/sequence.
 Returns  : An integer.
 Args     : A gene/sequence ID,
            A Cluster ID, 
            Optionally a new count for this gene/Cluster pair and the sum of
			the lengths covered by the cluster sites

=cut

sub gene_cluster_count
{
    my ($self, $gene_id, $cluster_id, $count) = @_;

    return if !defined $gene_id || !defined $cluster_id;

    if (defined $count) {
        $self->{_gene_cluster_counts}->{$gene_id}->{$cluster_id} = $count;

        if ($count > 0) {
            $self->{_cluster_gene_exists}->{$cluster_id}->{$gene_id} = 1;
        }

        unless ($self->gene_exists($gene_id)) {
            $self->_add_gene($gene_id);
        }

        unless ($self->cluster_exists($cluster_id)) {
            $self->_add_cluster($cluster_id);
        }
    }

    if ($self->{_gene_cluster_counts}->{$gene_id}) {
        return $self->{_gene_cluster_counts}->{$gene_id}->{$cluster_id} || 0;
    }

    return 0;
}

=head2 gene_cluster_length

 Title    : gene_cluster_length
 Usage    : $length = $counts->gene_cluster_length($gene_id, $cluster_id);
            $counts->gene_cluster_length($gene_id, $cluster_id, $length);
 Function : Get/set the length covered by the given cluster sites for the given
			gene/sequence. 
 Returns  : An integer.
 Args     : A gene/sequence ID,
            A Cluster ID, 
            Optionally a new length for this gene/Cluster and the sum of
			the lengths covered by the cluster sites

=cut

sub gene_cluster_length
{
    my ($self, $gene_id, $cluster_id, $length) = @_;

    return if !defined $gene_id || !defined $cluster_id;

	if (!$self->gene_exists($gene_id)) {
		carp "Gene $gene_id must first be defined with an appropriate count\n";
		return;
	}
	
	if (!$self->cluster_exists($cluster_id)) {
		carp "Cluster $cluster_id must first be defined\n";
		return;
	}
	
    if (defined $length) {
		$self->{_gene_cluster_lengths}->{$gene_id}->{$cluster_id} = $length;
    }

    if ($self->{_gene_cluster_lengths}->{$gene_id}) {
        return $self->{_gene_cluster_lengths}->{$gene_id}->{$cluster_id} || 0;
    }

    return 0;
}

=head2 cluster_gene_count

 Title    : cluster_gene_count
 Usage    : $count = $counts->cluster_gene_count($cluster_id)
 Function : Get the count of the number of sequences/genes for which
            sites for the given Cluster were detected.
 Returns  : An integer.
 Args     : A Cluster ID. 

=cut

sub cluster_gene_count
{
    my ($self, $cluster_id) = @_;

    return if !$cluster_id;

    my $gene_count = 0;

    if ($self->{_cluster_gene_exists}->{$cluster_id}) {
        $gene_count = scalar(keys %{$self->{_cluster_gene_exists}->{$cluster_id}}) || 0;
    }


    return $gene_count;
}

=head2 cluster_gene_ids

 Title    : cluster_gene_ids
 Usage    : $ids = $counts->cluster_gene_ids($cluster_id)
 Function : Get the list of gene/sequence IDs for which sites for 
            the given Cluster were detected.
 Returns  : A ref to a list of gene/sequence IDs.
 Args     : A Cluster ID. 

=cut

sub cluster_gene_ids
{
    my ($self, $cluster_id) = @_;

    return if !$cluster_id;

    my @gene_ids;
    if ($self->{_cluster_gene_exists}->{$cluster_id}) {
        @gene_ids = keys %{$self->{_cluster_gene_exists}->{$cluster_id}};
    }

    return @gene_ids ? \@gene_ids : undef;
}

=head2 cluster_count

 Title    : cluster_count
 Usage    : $count = $counts->cluster_count($cluster_id)
 Function : For the given Cluster, return the total number of cluster which appear
            for all the genes/sequences in the counts object.
 Returns  : An integer.
 Args     : A Cluster ID. 

=cut

sub cluster_count
{
    my ($self, $cluster_id) = @_;

    return if !$cluster_id;

    my $count = 0;
    if ($self->{_cluster_gene_exists}->{$cluster_id}) {
        my @gene_ids = keys %{$self->{_cluster_gene_exists}->{$cluster_id}};
        foreach my $gene_id (@gene_ids) {
            $count += $self->gene_cluster_count($gene_id, $cluster_id);
        }
    }

    return $count;
}

=head2 cluster_length

 Title    : cluster_length
 Usage    : $length = $counts->cluster_length($cluster_id)
 Function : For the given Cluster, return the total length of the cluster sites
			which appear on all the genes/sequences in the counts object.
 Returns  : An integer.
 Args     : A Cluster ID. 

=cut

sub cluster_length
{
    my ($self, $cluster_id) = @_;

    return if !$cluster_id;

    my $length = 0;
    if ($self->{_cluster_gene_exists}->{$cluster_id}) {
        my @gene_ids = keys %{$self->{_cluster_gene_exists}->{$cluster_id}};
        foreach my $gene_id (@gene_ids) {
            $length += $self->gene_cluster_length($gene_id, $cluster_id);
        }
    }

    return $length;
}

=head2 missing_gene_ids

 Title    : missing_gene_ids
 Usage    : $ids = $counts->missing_gene_ids()
 Function : Get a list of missing gene/sequence IDs. For convenience,
            the counts object allows storage of genes which may have
            been entered for analysis but could not be found in the
            database.
 Returns  : A ref to a list of gene/sequence IDs.
 Args     : None.

=cut

sub missing_gene_ids
{
    $_[0]->{-missing_gene_ids};
}

=head2 missing_cluster_ids

 Title    : missing_cluster_ids
 Usage    : $ids = $counts->missing_cluster_ids()
 Function : Get a list of missing Cluster IDs. For convenience, the counts
            object allows storage of Clusters which may have been entered
            for analysis but could not be found in the database.
 Returns  : A ref to a list of Cluster IDs.
 Args     : None.

=cut

sub missing_cluster_ids
{
    $_[0]->{-missing_cluster_ids};
}

=head2 subset

 Title    : subset
 Usage    : $subset = $counts->subset(
				    -gene_ids   => $gene_ids,
				    -cluster_ids     => $cluster_ids
            );

            OR

            $subset = $counts->subset(
                    -gene_start => $gene_start,
                    -gene_end   => $gene_end,
                    -cluster_start   => $tf_start,
                    -cluster_end     => $tf_end
            );

            OR some combination of above.

 Function : Get a new counts object from the current one with only a
            subset of the Clusters and/or genes. The subset of Clusters and 
            genes may be either specified with explicit ID lists or by
            starting and/or ending IDs.
 Returns  : A new OPOSSUM::Analysis::Cluster::Counts object.
 Args     : gene_ids    - optionally a list ref of gene IDs,
            cluster_ids      - optionally a list ref of Cluster IDs,
            gene_start  - optionally the starting gene ID,
            gene_end    - optionally the ending gene ID,
            cluster_start    - optionally the starting cluster ID,
            cluster_end      - optionally the ending cluster ID.

=cut

sub subset
{
    my ($self, %args) = @_;

    my $gene_start = $args{-gene_start};
    my $gene_end   = $args{-gene_end};
    my $gene_ids   = $args{-gene_ids};
    my $cl_start   = $args{-cluster_start};
    my $cl_end     = $args{-cluster_end};
    my $cluster_ids     = $args{-cluster_ids};

    my $all_gene_ids = $self->gene_ids;
    my $all_cluster_ids   = $self->cluster_ids;

    my $subset_gene_ids;
    my $subset_cluster_ids;

    my @missing_gene_ids;
    my @missing_cluster_ids;

    if (!defined $gene_ids) {
        if (!$gene_start && !$gene_end) {
            $subset_gene_ids = $all_gene_ids;
        } else {
            if (!$gene_start) {
                $gene_start = $all_gene_ids->[0];
            }
            if (!$gene_end) {
                $gene_end = $all_gene_ids->[scalar @$all_gene_ids - 1];
            }
            foreach my $gene_id (@$all_gene_ids) {
                if ($gene_id ge $gene_start && $gene_id le $gene_end) {
                    push @$subset_gene_ids, $gene_id;
                }
            }
        }
    } else {
        foreach my $gene_id (@$gene_ids) {
            if (grep(/^$gene_id$/, @$all_gene_ids)) {
                push @$subset_gene_ids, $gene_id;
            } else {
                carp "warning: gene ID $gene_id not in super set,"
                    . " omitting from subset";
                push @missing_gene_ids, $gene_id;
            }
        }
    }

    if (!defined $cluster_ids) {
        if (!$cl_start && !$cl_end) {
            $subset_cluster_ids = $all_cluster_ids;
        } else {
            if (!$cl_start) {
                $cl_start = $all_cluster_ids->[0];
            }
            if (!$cl_end) {
                $cl_end = $all_cluster_ids->[scalar @$all_cluster_ids - 1];
            }
            foreach my $cluster_id (@$all_cluster_ids) {
                if ($cluster_id ge $cl_start && $cluster_id le $cl_end) {
                    push @$subset_cluster_ids, $cluster_id;
                }
            }
        }
    } else {
        foreach my $cluster_id (@$cluster_ids) {
            if (grep(/^$cluster_id$/, @$all_cluster_ids)) {
                push @$subset_cluster_ids, $cluster_id;
            } else {
                carp "warning: Cluster ID $cluster_id not in super set,"
                    . " omitting from subset";
                push @missing_cluster_ids, $cluster_id;
            }
        }
    }

    my $subset = OPOSSUM::Analysis::Cluster::Counts->new();

    return if !$subset;

    foreach my $gene_id (@$subset_gene_ids) {
        foreach my $cluster_id (@$subset_cluster_ids) {
            $subset->gene_cluster_count(
                $gene_id,
                $cluster_id,
                $self->gene_cluster_count($gene_id, $cluster_id)
            );
			$subset->gene_cluster_length(
				$gene_id,
				$cluster_id,
				$self->gene_cluster_length($gene_id, $cluster_id)
			);
        }
    }

    $subset->{-missing_gene_ids} =
        @missing_gene_ids ? \@missing_gene_ids : undef;
    $subset->{-missing_cluster_ids} = @missing_cluster_ids ? \@missing_cluster_ids : undef;

    return $subset;
}

sub _add_gene
{
    my ($self, $gene_id) = @_;
    
    return if !defined $gene_id;

    unless ($self->{-gene_exists}->{$gene_id}) {
        push @{$self->{-gene_ids}}, $gene_id;
        $self->{-gene_exists}->{$gene_id} = 1;   
    }
}

sub _add_cluster
{
    my ($self, $cluster_id) = @_;
    
    return if !defined $cluster_id;

    unless ($self->{-cluster_exists}->{$cluster_id}) {
        push @{$self->{-cluster_ids}}, $cluster_id;
        $self->{-cluster_exists}->{$cluster_id} = 1;   
    }
}

1;
