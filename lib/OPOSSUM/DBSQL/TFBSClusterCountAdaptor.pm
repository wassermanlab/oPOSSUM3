
=head1 NAME

OPOSSUM::DBSQL::TFBSClusterCountsAdaptor - Adaptor for MySQL queries to retrieve
and store TFBS cluster counts.

=head1 SYNOPSIS

$tfbsca = $db_adaptor->get_TFBSClusterCountsAdaptor();

=head1 DESCRIPTION

In order to facilitate fast retrieval of TFBS cluster counts from the oPOSSUM database
several count sets were pre-computed using discrete values for PWM (PSSM)
thresholds, conservation levels, and upstream/downstream search regions.
This adaptor provides an interface to these pre-computed counts stored in
the tfbs_cluster_counts table.

=head1 AUTHOR

 Andrew Kwon, adapting code by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::TFBSClusterCountAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::TFBSClusterCount;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

sub new
{
    my ($class, $dbobj, $ext) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new($dbobj);

    $self->{-extended_table_name} = $ext;

    return $self;
}

=head2 fetch_where

 Title    : fetch_where
 Usage    : $counts = $tfbsca->fetch_where($where);
 Function : Fetch a single value or list of counts using a where clause
 Returns  : A single or ref to a list of OPOSSUM::TFBSClusterCount objects
 Args	  : $where = where clause

=cut

sub fetch_where
{
    my ($self, $where) = @_;

    if ($where and !$where =~ /^\s*where/) {
        $where = "where $where";
    }

    my $sql = qq{
        select gene_id, cluster_id, conservation_level, threshold_level,
        search_region_level, count, sum_length from tfbs_cluster_counts};
	$sql .= $where if $where;
	$sql .= "order by gene_id, cluster_id";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch TFBS cluster counts with\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch TFBS cluster counts with\n$sql\n" . $self->errstr;
        return;
    }

    #
    # Maybe we should return some sort of TFBSClusterCountSet?
    #
    my @counts;
    while (my @row = $sth->fetchrow_array) {
        push @counts, OPOSSUM::TFBSClusterCount->new(
            -gene_id             => $row[0],
            -cluster_id          => $row[1],
            -conservation_level  => $row[2],
            -threshold_level     => $row[3],
            -search_region_level => $row[4],
            -count               => $row[5],
			-sum_length			 => $row[6]
        );
    }
    $sth->finish;

    return undef if !@counts || !$counts[0];

    if (scalar @counts > 1) {
        return \@counts;
    } else {
        return $counts[0];
    }
}

=head2 fetch_values_where

 Title    : fetch_values_where
 Usage    : $count(s) = $tfbsca->fetch_values_where($where);
 Function : Fetch count(s).
 Returns  : If a single value is retrieved (count for one gene and one TF
            cluster), return a single integer count value. If counts for a
            single gene and one or more clusters are retrieved, return a ref to
            a hash with cluster ID as key and count as value. If counts for
            multiple genes and one or more clusters are retrieved, return a
            ref to a hash of hashes. i.e. $ref->{$gene_id}->{$cluster_id}->$count

 Args	  : $where = where clause

=cut

sub fetch_values_where
{
    my ($self, $where) = @_;

    if (!$where =~ /^\s*where/) {
        $where = "where $where";
    }

    my $sql = qq{select gene_id, cluster_id, count from tfbs_cluster_counts $where};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch TFBS cluster counts with\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch TFBS cluster counts with\n$sql\n" . $self->errstr;
        return;
    }

    my %counts;
    while (my @row = $sth->fetchrow_array) {
        $counts{$row[0]}->{$row[1]} = $row[2];
    }
    $sth->finish;

    return undef if !%counts;

    my @gids = keys %counts;
    if (scalar @gids > 1) {
        # Counts for multiple genes fetched.
        # Return ref to (gene) hash of (cluster) hashes of counts
        return \%counts;
    } else {
        # Counts for single genes fetched.
        my $gid = $gids[0];
        my @clusterids = keys %{$counts{$gid}};
        if (scalar @clusterids > 1) {
            # Counts for multiple TFs fetched.
            # Return ref to TF hash.
            return $counts{$gid};
        } else {
            # Counts for single cluster fetched. Return count value.
            my $cid = $clusterids[0];
            return $counts{$gid}->{$cid};
        }
    }
}

=head2 fetch_cluster_ids

 Title    : fetch_cluster_ids
 Usage    : $ids = $tfbsca->fetch_cluster_ids();
 Function : Fetch all the TF cluster IDs in the tfbs_cluster_counts table
 Returns  : A ref to a list of TF cluster IDs
 Args	  : None

=cut

sub fetch_cluster_ids
{
    my ($self) = @_;

    my $sql = qq{select distinct cluster_id from tfbs_cluster_counts};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TF cluster IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TF cluster IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }

    return @ids ? \@ids : undef;
}

=head2 fetch_gene_ids

 Title    : fetch_gene_ids
 Usage    : $ids = $tfbsca->fetch_gene_ids();
 Function : Fetch all the gene IDs in the tfbs_cluster_counts table
 Returns  : A ref to a list of gene IDs
 Args	  : None

=cut

sub fetch_gene_ids
{
    my ($self) = @_;

    my $sql = qq{select distinct gene_id from tfbs_cluster_counts};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching gene IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching gene IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }

    return @ids ? \@ids : undef;
}

=head2 fetch_tfbs_cluster_count

 Title    : fetch_tfbs_cluster_count
 Usage    : $count = $tfbsca->fetch_tfbs_cluster_count(
				-gene_id		        => $gid,
				-cluster_id			        => $cid,
 				-conservation_level	    => $clevel,
				-threshold_level	    => $thlevel,
				-search_region_level	=> $srlevel);
 Function : Fetch a TFBS cluster count for a given gene and TF cluster at the specified
            conservation, threshold and search region levels
 Returns  : An OPOSSUM::TFBSClusterCount object
 Args	  : gene_id             - gene ID
            cluster_id		    - TF cluster ID
            conservation_level	- conservation level
            threshold_level	    - threshold level
            search_region_level	- search region level

=cut

sub fetch_tfbs_cluster_count
{
    my ($self, %args) = @_;

    my $gid          = $args{-gene_id};
    my $cid         = $args{-cluster_id};
    my $cons_level   = $args{-conservation_level};
    my $thresh_level = $args{-threshold_level};
    my $sr_level     = $args{-search_region_level};

    if (!$gid || !$cid || !$cons_level || !$thresh_level || !$sr_level) {
        carp "must provide gene_id, cluster_id, conservation_level,"
            . " threshold_level and search_region_level\n";
        return;
    }

    my $where = qq{
        gene_id = $gid
        and cluster_id = $cid
        and conservation_level = $cons_level
        and threshold_level = $thresh_level
        and search_region_level = $sr_level
    };

    return $self->fetch_where($where);
}

=head2 fetch_gene_tfbs_cluster_counts

 Title    : fetch_gene_tfbs_cluster_counts
 Usage    : $counts = $tfbsca->fetch_gene_tfbs_cluster_counts(
				-gene_id                => $gids,
 				-conservation_level	    => $clevel,
				-threshold_level        => $thlevel,
				-search_region_level    => $srlevel,
				-cluster_ids                 => $clusterids);
 Function : Fetch TFBS cluster counts for a given gene at a specified conservation,
            threshold and search region level for a set of TFs
 Returns  : A ref to a list of OPOSSUM::TFBSClusterCount objects
 Args	  : gene_id	            - gene ID
            conservation_level	- conservation level
            threshold_level	    - threshold level
            search_region_level	- search region level
            cluster_ids         - optionally restrict the TFBS cluster counts to
                                  only those TF clusters identified by
                                  these TF cluster IDs

=cut

sub fetch_gene_tfbs_cluster_counts
{
    my ($self, %args) = @_;

    my $gid          = $args{-gene_id};
    my $cons_level   = $args{-conservation_level};
    my $thresh_level = $args{-threshold_level};
    my $sr_level     = $args{-search_region_level};
    my $clusterids   = $args{-cluster_ids};

    if (!$gid || !$cons_level || !$thresh_level || !$sr_level) {
        carp "must provide gene ID, conservation_level, threshold_level"
            . " and search_region_level\n";
        return;
    }

    my $where = qq{
        where gene_id = $gid
        and conservation_level = $cons_level
        and threshold_level = $thresh_level
        and search_region_level = $sr_level
    };
 
    if ($clusterids && $clusterids->[0]) {
        $where .= " and cluster_id in (";
        $where .= join(',', @$clusterids);
        $where .= ")";
    }

    return $self->fetch_where($where);
}

=head2 fetch_cluster_tfbs_cluster_counts

 Title    : fetch_tf_tfbs_cluster_counts
 Usage    : $counts = $tfbsca->fetch_cluster_tfbs_cluster_counts(
				-gene_id                => $gids,
 				-conservation_level	    => $clevel,
				-threshold_level        => $thlevel,
				-search_region_level    => $srlevel,
				-cluster_ids                 => $clusterids);
 Function : Fetch TFBS cluster counts for a given TF cluster at a specified conservation,
            threshold and search region level for a set of genes
 Returns  : A ref to a list of OPOSSUM::TFBSClusterCount objects
 Args	  : cluster_id	            - TF cluster ID
            conservation_level	- conservation level
            threshold_level	    - threshold level
            search_region_level	- search region level
            gene_ids            - optionally restrict the TFBS cluster counts to
                                  only those genes identified by these
                                  gene IDs

=cut

sub fetch_cluster_tfbs_cluster_counts
{
    my ($self, %args) = @_;

    my $cid          = $args{-cluster_id};
    my $cons_level   = $args{-conservation_level};
    my $thresh_level = $args{-threshold_level};
    my $sr_level     = $args{-search_region_level};
    my $gids         = $args{-gene_ids};

    if (!$cid || !$cons_level || !$thresh_level || !$sr_level) {
        carp "must provide TF cluster ID, conservation_level, threshold_level"
            . " and search_region_level\n";
        return;
    }

    my $where = qq{
        cluster_id = $cid
        and conservation_level = $cons_level
        and threshold_level = $thresh_level
        and search_region_level = $sr_level
    };
 
    if ($gids && $gids->[0]) {
        $where .= " and gene_id in (";
        $where .= join(',', @$gids);
        $where .= ")";
    }

    return $self->fetch_where($where);
}

=head2 fetch_tfbs_cluster_counts

 Title    : fetch_tfbs_cluster_counts
 Usage    : $counts = $tfbsca->fetch_tfbs_cluster_counts(
 				-conservation_level	    => $clevel,
				-threshold_level	    => $thlevel,
				-search_region_level	=> $srlevel,
				-gene_ids		        => $gids,
				-cluster_ids			        => $clusterids);
 Function : Fetch TFBS cluster counts at a specified conservation, threshold
            and search region level for a set of genes and TF clusters
 Returns  : A ref to a list of OPOSSUM::TFBSClusterCount objects
 Args	  : conservation_level	- conservation level
            threshold_level	    - threshold level
            search_region_level	- search region level
            gene_ids	        - optionally restrict the TFBS cluster counts to
                                  genes identified by these gene IDs
            cluster_ids         - optionally restrict the TFBS cluster counts to
                                  only TF clusters identified by these TF cluster IDs

=cut

sub fetch_tfbs_cluster_counts
{
    my ($self, %args) = @_;

    my $cons_level   = $args{-conservation_level};
    my $thresh_level = $args{-threshold_level};
    my $sr_level     = $args{-search_region_level};
    my $gids         = $args{-gene_ids};
    my $clusterids   = $args{-cluster_ids};

    if (!$cons_level || !$thresh_level || !$sr_level) {
        carp "must provide conservation_level, threshold_level"
            . " and search_region_level\n";
        return;
    }

    my $where = qq{
        conservation_level = $cons_level
        and threshold_level = $thresh_level
        and search_region_level = $sr_level
    };

    if ($gids && $gids->[0]) {
        $where .= " and gene_id in (";
        $where .= join(',', @$gids);
        $where .= ")";
    }

    if ($clusterids && $clusterids->[0]) {
        $where .= " and cluster_id in (";
        $where .= join(',', @$clusterids);
        $where .= ")";
    }

    return $self->fetch_where($where);
}

=head2 store

 Title   : store
 Usage   : $tfbsca->store($count);
 Function: Store TFBS cluster count in the database.
 Args    : An OPOSSUM::TFBSClusterCount object
 Returns : True on success, false otherwise.

=cut

sub store
{
    my ($self, $tfbsc) = @_;

    my $sql = qq{
        insert into tfbs_cluster_counts
        (gene_id, cluster_id, conservation_level, threshold_level,
         search_region_level, count)
        values (?, ?, ?, ?, ?, ?)
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert tfbs_cluster_counts statement\n" . $self->errstr;
        return;
    }

    if (!$sth->execute(
            $tfbsc->gene_id,             $tfbsc->cluster_id,
            $tfbsc->conservation_level,  $tfbsc->threshold_level,
            $tfbsc->search_region_level, $tfbsc->count
        )
        )
    {
        carp sprintf(
            "Error inserting TFBS cluster count for Gene ID %d TF cluster ID %d: "
            . $self->errstr,
            $tfbsc->gene_id, $tfbsc->cluster_id
        );
        return 0;
    }

    return 1;
}

=head2 store_list

 Title   : store_list
 Usage   : $tfbsca->store_list($counts);
 Function: Store TFBS cluster counts in the database.
 Args    : A listref of OPOSSUM::TFBSClusterCount objects
 Returns : True on success, false otherwise.

=cut

sub store_list
{
    my ($self, $tfbs_cluster_counts) = @_;

    my $sql = qq{
        insert into tfbs_cluster_counts
        (gene_id, cluster_id, conservation_level, threshold_level,
         search_region_level, count, sum_length)
        values (?, ?, ?, ?, ?, ?, ?)
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert tfbs_cluster_counts statement\n" . $self->errstr;
        return;
    }

    my $ok = 1;
    foreach my $tfbsc (@$tfbs_cluster_counts) {
        if (!$sth->execute(
                $tfbsc->gene_id,             $tfbsc->cluster_id,
                $tfbsc->conservation_level,  $tfbsc->threshold_level,
                $tfbsc->search_region_level, $tfbsc->count,
				$tfbsc->sum_length
            )
            )
        {
            carp sprintf(
                "Error inserting TFBS cluster count for Gene ID %d TF cluster ID %d: "
                . $self->errstr,
                $tfbsc->gene_id, $tfbsc->cluster_id);
            $ok = 0;
        }
    }

    return $ok;
}

1;
