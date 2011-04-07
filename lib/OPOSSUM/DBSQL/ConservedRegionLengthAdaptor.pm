=head1 NAME

OPOSSUM::DBSQL::ConservedRegionLengthAdaptor - Adaptor for MySQL queries to
retrieve and store conserved region lengths.

=head1 SYNOPSIS

$cra = $db_adaptor->get_ConservedRegionLengthAdaptor();

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 CHANGE HISTORY

 DJA 2010/11/24
 - cleaned up operon handling in _fetch_length_default and _fetch_length_custom

 AK 2010/10/28
 - added arguments -operon_gene_ids and -has_operon to all variations of
   fetch_counts(), fetch_custom_counts(), and fetch_anchored_counts()
 - Now, if input gene ids contain operon genes, the counts from the first
   gene in the operon are duplicated for all other genes belonging to the
   same operon
 
 Shannan Ho Sui on Dec 29, 2006

=head1 METHODS

=cut
package OPOSSUM::DBSQL::ConservedRegionLengthAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::ConservedRegionLengthSet;
use OPOSSUM::ConservedRegionLength;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title    : new
 Usage    : $crla = OPOSSUM::DBSQL::ConservedRegionLengthAdaptor->new(@args);
 Function : Create a new ConservedRegionLengthAdaptor.
 Returns  : A new OPOSSUM::DBSQL::ConservedRegionLengthAdaptor object.
 Args	  : An OPOSSUM::DBSQL::DBConnection object.

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_gene_ids

 Title    : fetch_gene_ids
 Usage    : $ids = $crla->fetch_gene_ids();
 Function : Fetch a list of all the gene IDs for the conserved
 	    region lengths in the current DB.
 Returns  : A reference to a list of gene IDs.
 Args	  : None.

=cut

sub fetch_gene_ids
{
    my ($self) = @_;

    my $sql = qq{select distinct gene_id from conserved_region_lengths};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error fetching gene IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error fetching gene IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }

    return @ids ? \@ids : undef;
}

=head2 fetch_where

 Title    : fetch_where
 Usage    : $lengths = $crla->fetch_where($where);
 Function : Fetch a list of conserved region lengths based on the where
            clause.
 Returns  : A reference to a list of oPOSSUM::ConservedRegionLength objects.
 Args	  : Optional where clause.

=cut

sub fetch_where
{
    my ($self, $where) = @_;

    my $sql = qq{
		select gene_id, conservation_level, search_region_level, length
		from conserved_region_lengths
    };

    if ($where) {
        unless ($where =~ /^\s*where /) {
            $where = 'where ' . $where;
        }

        $sql .= " $where";
    }

    $sql .= " order by gene_id";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch conserved region lengths with\n$sql\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute()) {
        carp "Error executing fetch conserved region lengths with\n$sql\n"
            . $self->errstr;
        return;
    }

    my @lengths;
    while (my @row = $sth->fetchrow_array()) {
        push @lengths, OPOSSUM::ConservedRegionLength->new(
            -gene_id                => $row[0],
            -conservation_level     => $row[1],
            -search_region_level    => $row[2],
            -length                 => $row[3]
        );
    }
    $sth->finish;

    return @lengths ? \@lengths : undef;
}

=head2 fetch_length_list

 Title    : fetch_length_list
 Usage    : $lengths = $crla->fetch_length_list(
 				-gene_ids               => $gpids,
				-conservation_level     => $clevel,
				-search_region_level    => $srlevel
            );
 Function : Fetch a list of the the conserved region lengths for the
            given list of gene IDs at the given level of conservation
            and search region level.
 Returns  : A reference to a list of OPOSSUM::ConservedRegionLength objects.
 Args	  : Optional list ref of gene IDS;
            Optional integer conservation level;
            Optional integer search region level.

=cut

sub fetch_length_list
{
    my ($self, %args) = @_;

    my $gpids = $args{-gene_ids};
    my $cons_level = $args{-conservation_level};
    my $sr_level = $args{-search_region_level};

    my $where = "";
    my $clause = "where";
    if (defined $cons_level) {
    	$where .= " $clause conservation_level = $cons_level";
        $clause = "and";
    }

    if (defined $sr_level) {
        $where .= " $clause search_region_level = $sr_level";
        $clause = "and";
    }

    if ($gpids && $gpids->[0]) {
        $where .= " $clause gene_id in (";
        $where .= join(',', @$gpids);
        $where .= ")";
    }

    return $self->fetch_where($where);
}

=head2 fetch_length_set

 Title    : fetch_length_set
 Usage    : $len_set = $crla->fetch_length_set(
                -gene_ids               => $gpids,
                -conservation_level     => $clevel,
                -search_region_level    => $srlevel
            );
 Function : Fetch the set of conserved region lengths for the given list
            of gene IDs at the given level of conservation and
            search region level.
 Returns  : An OPOSSUM::ConservedRegionLengthSet.
 Args	  : Optional list ref of gene IDS;
            Optional integer conservation level;
            Optional integer search region level.

=cut

sub fetch_length_set
{
    my ($self, %args) = @_;

    my $list = $self->fetch_length_list(%args);

    return if !$list;

    my $set = OPOSSUM::ConservedRegionLengthSet->new();
    foreach my $crl (@$list) {
        $set->add_conserved_region_length($crl);
    }

    if (defined $args{-conservation_level}) {
        $set->param('conservation_level', $args{-conservation_level});
    }

    if (defined $args{-search_region_level}) {
        $set->param('search_region_level', $args{-search_region_level});
    }

    return $set;
}

=head2 fetch_total_length

 Title    : fetch_total_length
 Usage    : $length = $crla->fetch_total_length(
 				-gene_ids               => $gpids,
				-conservation_level     => $clevel,
				-search_region_level    => $srlevel
            );

            $length = $crla->fetch_total_length(
 				-gene_ids               => $gpids,
				-conservation_level     => $clevel,
				-upstream_bp            => $upstream_bp,
				-downstream_bp          => $downstream_bp
            );
 Function : Fetch the total conserved region length for the given list of
            gene IDs at the given level of conservation and search region
            defined by either a search region level or upstream/downstream
            bp.
 Returns  : The total conserved region length.
 Args	  : Optional list ref of gene IDS;
            Optional integer conservation level;
            Optional integer search region level.
			Optional hash ref of operon gene id mapping to first operon gene
			Optional has_operon
=cut

sub fetch_total_length
{
    my ($self, %args) = @_;

    my $gene_ids = $args{-gene_ids};
	
    my $sr_level = $args{-search_region_level};
    my $cons_level = $args{-conservation_level};
    my $upstream_bp = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};
	
	my $operon_gids  = $args{-operon_gene_ids};
	my $has_operon   = $args{-has_operon};
	
    unless ($cons_level) {
        carp "Must specify conservation level\n";
        return;
    }

    if ($sr_level) {
        return $self->_fetch_length_default(
			$cons_level, $sr_level, $gene_ids, $operon_gids, $has_operon
		);
    } elsif (defined $upstream_bp || defined $downstream_bp) {
        return $self->_fetch_length_custom(
            $cons_level, $upstream_bp, $downstream_bp, $gene_ids, $operon_gids, $has_operon
        );
    } else {
        carp  "Must specify either search region level or"
            . " upstream/downstream bp\n";
        return;
    }
}

#
# Fetch total length of conserved regions for input genes at given conservation
# and search region level from precomputed lengths in conserged_region_lengths
# table.
#
sub _fetch_length_default
{
    my ($self, $cons_level, $sr_level, $gids, $operon_gids, $has_operon) = @_;

    if ($has_operon && !$gids) {
        $operon_gids = {}; # reset this.

        my $ga = $self->db->get_GeneAdaptor();
        if (!$ga) {
            carp "error getting GeneAdaptor\n";
            return;
        }

        $gids = $ga->fetch_gene_ids();
        
        my $oa = $self->db->get_OperonAdaptor();
        if (!$oa) {
            carp "error getting OperonAdaptor\n";
            return;
        }

        my $operons = $oa->fetch_where();
        
        foreach my $operon (@$operons) {
            my $fgene = $operon->fetch_first_gene();
            foreach my $opgene (@{$operon->genes()}) {
                $operon_gids->{$opgene->id} = $fgene->id;
                #	if $opgene->id != $fgene->id;
            }
        }
    }

    my $t_gids;
    if ($operon_gids) {
        my %gid_inc;
        foreach my $gid (@$gids) {
            my $fgid = $operon_gids->{$gid};
            if ($fgid) {
                unless ($gid_inc{$fgid}) {
                    push @$t_gids, $fgid;
                    $gid_inc{$fgid} = 1;
                }
            } else {
                push @$t_gids, $gid;
            }
        }
    } else {
        if ($gids) {
            $t_gids = $gids;
        }
    }

    my $sql = qq{select gene_id, length from conserved_region_lengths};

    my $where =   " where conservation_level = $cons_level"
                . " and search_region_level = $sr_level";

    if ($t_gids && $t_gids->[0]) {
        $where .= " and gene_id in (";
        $where .= join(',', @$t_gids);
        $where .= ")";
	}
	
	$sql .= $where;
	
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch total conserved region length with\n$sql\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute()) {
        carp "Error executing fetch total conserved region length with\n$sql\n"
            . $self->errstr;
        return;
    }

    my %gid_lengths;
    while (my ($gid, $length) = $sth->fetchrow_array()) {
        $gid_lengths{$gid} = $length;
    }

    $sth->finish();

    my $length = 0;
    if ($gids && $gids->[0]) {
        foreach my $gid (@$gids) {
            if ($operon_gids) {
                my $fgid = $operon_gids->{$gid};
                if ($fgid) {
                    $length += $gid_lengths{$fgid} if $gid_lengths{$fgid};
                } else {
                    $length += $gid_lengths{$gid} if $gid_lengths{$gid};
                }
            } else {
                $length += $gid_lengths{$gid} if $gid_lengths{$gid};
            }
        }
    } else {
        foreach my $gid (keys %gid_lengths) {
            $length += $gid_lengths{$gid} if $gid_lengths{$gid};
        }
    }
	
    return $length;
}

#
# Fetch total length of conserved regions for input genes by reading from
# the conserved_regions table and calculating the amount of sequence which
# falls within the specified search region at the given level of conservation.
#
sub _fetch_length_custom
{
    my ($self, $conservation_level, $upstream_bp, $downstream_bp, $gids,
        $operon_gids, $has_operon) = @_;

    unless ($gids) {
        my $ga = $self->db()->get_GeneAdaptor();
        
        unless ($ga) {
            carp "Error getting GeneAdaptor";
            return;
        }

        $gids = $ga->fetch_gene_ids();

        unless ($gids) {
            carp "Error fetching gene IDs";
            return;
        }
		
		if ($has_operon) {
            $operon_gids = {}; # reset this.
            
			my $oa = $self->db->get_OperonAdaptor();
			if (!$oa) {
				carp "error getting OperonAdaptor\n";
				return;
			}

			my $operons = $oa->fetch_where();
			foreach my $operon (@$operons) {
				my $fgene = $operon->fetch_first_gene();
				foreach my $opgene (@{$operon->genes}) {
					$operon_gids->{$opgene->id} = $fgene->id;
					#	if $opgene->id != $fgene->id;
				}
			}
		}
    }

	#
	# If operon_gids is defined, go through the list and make sure non-operon
    # gene and first gene of any operon gene IDs are added to the input gene
    # IDs
	#
	my $t_gids;
	if ($operon_gids) {
        my %gid_inc;
		foreach my $gid (@$gids) {
            my $fgid = $operon_gids->{$gid};

            if ($fgid) {
                unless ($gid_inc{$fgid}) {
                    push @$t_gids, $fgid;
                    $gid_inc{$fgid} = 1;
                }
            } else {
                push @$t_gids, $gid;
            }
		}
	} else {
		$t_gids = $gids;
	}
	
    my $gid_lengths = $self->_fetch_custom_lengths_by_gene_ids(
        $conservation_level, $t_gids, $upstream_bp, $downstream_bp
    );

    my $length = 0;
    foreach my $gid (@$gids) {
        if ($operon_gids) {
            my $fgid = $operon_gids->{$gid};
			my $fgid_length = $gid_lengths->{$fgid};
            if ($fgid) {
                $length += $fgid_length if $fgid_length;
            } else {
                $length += $fgid_length if $fgid_length;
            }
        } else {
            $length += $gid_lengths->{$gid} if $gid_lengths->{$gid};
        }
    }

    return $length;
}

#
# Given a list of gene ids, upstream_bp, downstream_bp, return a hash ref
# of key = gene id, val = length
#
sub _fetch_custom_lengths_by_gene_ids
{
	my ($self, $conservation_level, $gids, $upstream_bp, $downstream_bp) = @_;
	
	my %lengths;

    my $ga = $self->db()->get_GeneAdaptor() || carp "Error getting GeneAdaptor";
	
	foreach my $gid (@$gids) {
        my $gene = $ga->fetch_by_id($gid);
        unless ($gene) {
            carp "Error fetching gene $gid";
            next;
        }

        my $psrs = $gene->promoter_search_regions($upstream_bp, $downstream_bp);
        unless ($psrs) {
            carp "Error fetching promoter search regions for gene $gid";
            next;
        }

        my $crs = $gene->conserved_regions($conservation_level);
        unless ($crs) {
            carp "No conserved regions for gene $gid at conservation level"
                . " $conservation_level";
            next;
        }
		
		my $length = 0;
        foreach my $psr (@$psrs) {
            my $psr_start = $psr->start();
            my $psr_end   = $psr->end();
            foreach my $cr (@$crs) {
                my $cr_start = $cr->start();
                my $cr_end   = $cr->end();

                next if $cr_end < $psr_start;
                next if $cr_start > $psr_end;

                #
                # Only include part of conserved region which falls within
                # promoter search region.
                #
                $cr_start = $psr_start if $cr_start < $psr_start;
                $cr_end   = $psr_end if $cr_end > $psr_end;

                $length += $cr_end - $cr_start + 1;
            }
        }
		$lengths{$gid} = $length;
    }
	
	return %lengths ? \%lengths : undef;
}

=head2 store

 Title    : store
 Usage    : $crla->store($crl);
 Function : Store a ConservedRegionLength object in the DB.
 Returns  : True on success, false otherwise.
 Args	  : An OPOSSUM::ConservedRegionLength

=cut

sub store
{
    my ($self, $crl) = @_;

    my $sql = qq{insert into conserved_region_lengths
		    (gene_id, conservation_level, search_region_level,
		     length)
		values (?, ?, ?, ?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
    	carp "Error preparing insert conserved_region_lengths statement\n"
		. $self->errstr;
	return;
    }

    if (!$sth->execute($crl->gene_id, $crl->conservation_level,
    			$crl->search_region_level, $crl->length))
    {
	carp sprintf("Error inserting conserved region length for gene ID"
		    . " %d " . $self->errstr,
		    $crl->gene_id);
	return 0;
    }

    return 1;
}

1;
