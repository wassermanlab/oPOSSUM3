=head1 NAME

OPOSSUM::DBSQL::ConservedTFBSAdaptor - Adaptor for MySQL queries to retrieve
and store conserved TFBSs.

=head1 SYNOPSIS

 $ctfbsa = $db_adaptor->get_ConservedTFBSAdaptor();

=head1 MODIFICATIONS
 2010/03/04
 - This adaptor now retrieves TFBSs as TFBS::Site objects instead of
   OPOSSUM::ConservedTFBS objects.
 2010/08/27 by Andrew Kwon
 - added -tf_id argument to fetch_set_by_gene_id so that you can retrieve the sites for
   that TF only
 - added fetch_set_by_gene_tf_id, which is a synonym of fetch_set_by_gene_id
   (more intuitive name if you want to retrieve sites for a specific TF and gene)
   
=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ConservedTFBSAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::TFBSCount;
use OPOSSUM::ConservedTFBS;
use OPOSSUM::ConservedTFBSSet;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;
use TFBS::Site;
use TFBS::SiteSet;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title    : new
 Usage    : $ctfsa = OPOSSUM::DBSQL::ConservedTFBSAdaptor->new(@args);
 Function : Create a new ConservedTFBSAdaptor.
 Returns  : A new OPOSSUM::DBSQL::ConservedTFBSAdaptor object.
 Args	  : An OPOSSUM::DBSQL::DBConnection object.

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_tf_ids

 Title    : fetch_tf_ids
 Usage    : $ids = $ctfsa->fetch_tf_ids();
 Function : Fetch a list of all the distinct TF IDs for all the conserved
            TFBSs in the current DB.
 Returns  : A list ref of integer TF IDs.
 Args	  : None.

=cut

sub fetch_tf_ids
{
    my ($self) = @_;

    my $sql = qq{select distinct tf_id from conserved_tfbss};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching TF IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching TF IDs\n" . $self->errstr;
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
 Usage    : $ids = $ctfsa->fetch_gene_ids();
 Function : Fetch a list of all the distinct gene IDs for all
            the conserved TFBSs in the current DB.
 Returns  : A list ref of integer gene IDs.
 Args	  : None.

=cut

sub fetch_gene_ids
{
    my ($self) = @_;

    my $sql = qq{select distinct gene_id from conserved_tfbss};

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

=head2 fetch_gene_tfbs_count

 Title    : fetch_gene_tfbs_count
 Usage    : $count = $ctfsa->fetch_gene_tfbs_count(
                 -gene_id		        => $gid,
                 -tf_id			        => $tfid,
                 -conservation_level	=> $clevel,
                 -threshold		        => $threshold,
                 -upstream_bp		    => $upstream_bp,
                 -downstream_bp         => $downstream_bp
            );
 Function : Fetch a count of the number of conserved TFBSs for the
            given TF found for the given gene at the given level of
            conservation, optional TFBS score threshold and optional
            search region.
 Returns  : An integer count of the total number of conserved binding sites
            of the given TF which fall within the promoter regions of the
            given gene.
 Args	  : Gene ID,
            TF ID,
            conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp

=cut

sub fetch_gene_tfbs_count
{
    my ($self, %args) = @_;

    my $gid           = $args{-gene_id};
    my $tfid          = $args{-tf_id};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    if (!$tfid) {
        carp "must provide TF ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_gene_id($gid);
    if (!$gene) {
        carp "error fetching Gene $gid";
        return;
    }

    my $strand = $gene->strand;

    my $sql = qq{
        select count(tf_id) from conserved_tfbss
        where gene_id = $gid
        and conservation_level >= $clevel
        and tf_id = '$tfid'
    };

    if ($threshold) {
        if ($threshold =~ /(.+)%$/) {
            $threshold = $1 / 100;
        }
            
        $sql .= " and rel_score >= $threshold";
    }

    my $search_regions = $gene->promoter_search_regions(
        $upstream_bp, $downstream_bp
    );

    if ($search_regions) {
        my $first_sr = 1;
        foreach my $sr (@$search_regions) {
            my $sr_start = $sr->start;
            my $sr_end   = $sr->end;
            if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
                $sql .= " and ((start >= $sr_start and end <= $sr_end)";
                $first_sr = 0;
            } else {
                $sql .= " or (start >= $sr_start and end <= $sr_end)";
            }
        }

        if (!$first_sr) {
            $sql .= ")";
        }
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch conserved TF site count with:\n$sql\n"
            . $self->errstr;
        return;
    }
    
    if (!$sth->execute) {
        carp "error executing fetch conserved TF site count with:\n$sql\n"
            . $self->errstr;
        return;
    }
    
    my $count = 0;
    if (my @row = $sth->fetchrow_array) {
        $count = $row[0];
    }

    return $count;
}

=head2 fetch_gene_tfbs_counts

 Title    : fetch_gene_tfbs_counts
 Usage    : $counts = $ctfsa->fetch_gene_tfbs_counts(
                -gene_id            => $gid,
                -tf_ids             => $tfids,
                -conservation_level => $clevel,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );
 Function : Fetch counts of the number of conserved binding sites for a
            given gene at a given level of conservation. Optionally limit
            to a given set of TFs, TFBS score threshold and search region.
 Returns  : A listref of OPOSSUM::TFBSCount objects.
 Args	  : Gene ID,
            Conservation level,
            OPTIONAL TF IDs,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp

=cut

sub fetch_gene_tfbs_counts
{
    my ($self, %args) = @_;

    my $gid           = $args{-gene_id};
    my $tfids         = $args{-tf_ids};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_gene_id($gid);
    if (!$gene) {
        carp "error fetching Gene $gid";
        return;
    }

    my $strand = $gene->strand;

    my $sql = qq{
        select tf_id, count(tf_id) from conserved_tfbss
        where gene_id = $gid
        and conservation_level >= $clevel
    };

    if ($threshold) {
        if ($threshold =~ /(.+)%$/) {
            $threshold = $1 / 100;
        }
            
        $sql .= " and rel_score >= $threshold";
    }

    my $search_regions = $gene->promoter_search_regions(
        $upstream_bp, $downstream_bp
    );

    if ($search_regions) {
        my $first_sr = 1;
        foreach my $sr (@$search_regions) {
            my $sr_start = $sr->start;
            my $sr_end   = $sr->end;
            if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
                $sql .= " and ((start >= $sr_start and end <= $sr_end)";
                $first_sr = 0;
            } else {
                $sql .= " or (start >= $sr_start and end <= $sr_end)";
            }
        }

        if (!$first_sr) {
            $sql .= ")";
        }
    }

    if ($tfids && $tfids->[0]) {
        $sql .= " and tf_id in ('";
        $sql .= join("','", @$tfids);
        $sql .= "')";
    }

    $sql .= " group by tf_id";

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch conserved TF site count with:\n$sql\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch conserved TF site counts with:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @counts;
    while (my @row = $sth->fetchrow_array) {
        push @counts, OPOSSUM::TFBSCount->new(
            -gene_id            => $gid,
            -conservation_level => $clevel,
            -tf_id              => $row[0],
            -count              => $row[1]
        );
    }

    return @counts ? \@counts : undef;
}

###  fetch_by_gene added by SHS on Dec 27,2006, modified on Jan 10, 2007
# modified by AK on Sept. 19, 2010
# can now take an array ref for -tf_ids parameter
=head2 fetch_by_gene

 Title    : fetch_by_gene
 Usage    : $count = $ctfbsa->fetch_by_gene(
                -gene_id            => $gid,
                -tf_id              => $tfid,
                -tf_ids             => $tfids,
                -conservation_level => $clevel,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );
 Function : Fetch a list of conserved TFBSs for the given TF found for
            the given gene at the given conservation level. Optionally
            limit to only those sites which score greater than or equal
            to the given threshold and/or fall within the given bounds
            specified by the number of upstream and downstream bp.
 Returns  : A list ref of OPOSSUM::ConservedTFBS objects.
 Args     : Gene ID,
            TF ID,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
=cut

sub fetch_by_gene
{
    my ($self, %args) = @_;

    my $gid           = $args{-gene_id};
    my $tfid          = $args{-tf_id};
    my $tfids         = $args{-tf_ids};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    if (!$tfid and !$tfids) {
        carp "must provide TF ID(s)";
        return;
    }
    
    if ($tfid and $tfids) {
        carp "only 1 of -tf_id or -tf_ids must be used";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_gene_id($gid);
    if (!$gene) {
        carp "error fetching Gene $gid";
        return;
    }

    my $sql = qq{
        select tf_id, start, end, strand, score, rel_score, seq,
		    conservation_level, conservation
		from conserved_tfbss 
		where gene_id = $gid
		    and conservation_level >= $clevel
    };
    
    if ($tfid) {
        $sql .= " and tf_id = '$tfid'";
    } elsif ($tfids) {
        $sql .= " and tf_id in ('";
        $sql .= join "','", @$tfids;
        $sql .= "')";
    }

    if ($threshold) {
        if ($threshold =~ /(.+)%$/) {
            $threshold = $1 / 100;
        }
            
        $sql .= " and rel_score >= $threshold";
    }

    my $search_regions = $gene->promoter_search_regions(
        $upstream_bp, $downstream_bp
    );

    if ($search_regions) {
        my $first_sr = 1;
        foreach my $sr (@$search_regions) {
            my $sr_start = $sr->start;
            my $sr_end   = $sr->end;
            if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
                $sql .= " and ((start >= $sr_start and end <= $sr_end)";
                $first_sr = 0;
            } else {
                $sql .= " or (start >= $sr_start and end <= $sr_end)";
            }
        }

        if (!$first_sr) {
            $sql .= ")";
        }
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch conserved TF sites with:\n$sql\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch conserved TF sites with:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @tfbss;
    while (my @row = $sth->fetchrow_array) {
        push @tfbss, OPOSSUM::ConservedTFBS->new(
            -gene_id            => $gid,
            -tf_id              => $row[0],
            -start              => $row[1],
            -end                => $row[2],
            -strand             => $row[3],
            -score              => $row[4],
            -rel_score          => $row[5],
            -seq                => $row[6],
            -conservation_level => $row[7],
            -conservation       => $row[8]
        );
    }
    $sth->finish;

    return @tfbss ? \@tfbss : undef;
}

=head2 fetch_all_by_gene_id

 Title    : fetch_all_by_gene_id
 Usage    : $tfbss = $ctfbsa->fetch_all_by_gene_id(
                 -gene_id               => $gid,
                 -conservation_level    => $clevel,
                 -threshold             => $threshold,
                 -upstream_bp           => $upstream_bp,
                 -downstream_bp         => $downstream_bp
            );
 Function : Fetch a list of all conserved TFBSs for the given gene
            ID at the given conservation level.
            Optionally limit to only those sites which score
            greater than or equal to the given threshold and/or fall
            within the given bounds specified by the number of upstream
            and downstream bp.
 Returns  : A list ref of OPOSSUM::ConservedTFBS objects.
 Args     : Gene ID,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
=cut

sub fetch_all_by_gene_id
{
    my ($self, %args) = @_;

    my $gid           = $args{-gene_id};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $pa = $self->db->get_PromoterAdaptor;
    if (!$pa) {
        carp "error getting PromoterAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_gene_id($gid);
    if (!$gene) {
        carp "error fetching Gene $gid";
        return;
    }

    my $sql = qq{
        select tf_id, start, end, strand, score, rel_score, seq,
		    conservation_level, conservation
		from conserved_tfbss 
		where gene_id = $gid and conservation_level >= $clevel
    };

    if ($threshold) {
        if ($threshold =~ /(.+)%$/) {
            $threshold = $1 / 100;
        }
            
        $sql .= " and rel_score >= $threshold";
    }

    my $search_regions = $gene->promoter_search_regions(
        $upstream_bp, $downstream_bp
    );

    if ($search_regions) {
        my $first_sr = 1;
        foreach my $sr (@$search_regions) {
            my $sr_start = $sr->start;
            my $sr_end   = $sr->end;
            if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
                $sql .= " and ((start >= $sr_start and end <= $sr_end)";
                $first_sr = 0;
            } else {
                $sql .= " or (start >= $sr_start and end <= $sr_end)";
            }
        }

        if (!$first_sr) {
            $sql .= ")";
        }
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch conserved TF sites with:\n$sql\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch conserved TF sites with:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @tfbss;
    while (my @row = $sth->fetchrow_array) {
        push @tfbss, OPOSSUM::ConservedTFBS->new(
            -tf_id              => $row[0],
            -gene_id            => $gid,
            -start              => $row[1],
            -end                => $row[2],
            -strand             => $row[3],
            -score              => $row[4],
            -rel_score          => $row[5],
            -seq                => $row[6],
            -conservation_level => $row[7],
            -conservation       => $row[8]
        );
    }
    $sth->finish;

    return @tfbss ? \@tfbss : undef;
}

=head2 fetch_set_by_gene_id

 Title    : fetch_set_by_gene_id
 Usage    : $tf_set = $ctfsa->fetch_set_by_gene_id($id, $level);
 Function : Fetch the set of conserved TF sites for the given gene
            ID at the given level of conservation.
 Returns  : An OPOSSUM::ConservedTFBSSet object.
 Args	  : Integer gene ID and integer conservation level.

=cut

sub fetch_set_by_gene_id
{
    my ($self, %args) = @_;

    my $gid           = $args{-gene_id};
    my $tfid          = $args{-tf_id};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $pa = $self->db->get_PromoterAdaptor;
    if (!$pa) {
        carp "error getting PromoterAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_gene_id($gid);
    if (!$gene) {
        carp "error fetching Gene $gid";
        return;
    }

    my $sql = qq{
        select tf_id, start, end, strand, score, rel_score, seq,
		    conservation_level, conservation
		from conserved_tfbss 
		where gene_id = $gid and conservation_level >= $clevel
    };

    if ($threshold) {
        if ($threshold =~ /(.+)%$/) {
            $threshold = $1 / 100;
        }
            
        $sql .= " and rel_score >= $threshold";
    }
    
    if ($tfid) {
        $sql .= " and tf_id = '$tfid'";
    }

    my $search_regions = $gene->promoter_search_regions(
        $upstream_bp, $downstream_bp
    );

    if ($search_regions) {
        my $first_sr = 1;
        foreach my $sr (@$search_regions) {
            my $sr_start = $sr->start;
            my $sr_end   = $sr->end;
            if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
                $sql .= " and ((start >= $sr_start and end <= $sr_end)";
                $first_sr = 0;
            } else {
                $sql .= " or (start >= $sr_start and end <= $sr_end)";
            }
        }

        if (!$first_sr) {
            $sql .= ")";
        }
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error preparing fetch conserved TF sites with:\n$sql\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error executing fetch conserved TF sites with:\n$sql\n"
            . $self->errstr;
        return;
    }

    my $tf_site_set = OPOSSUM::ConservedTFBSSet->new();
    while (my @row = $sth->fetchrow_array) {
        $tf_site_set->add_tf_site(
            OPOSSUM::ConservedTFBS->new(
                -tf_id              => $row[0],
                -gene_id            => $gid,
                -start              => $row[1],
                -end                => $row[2],
                -strand             => $row[3],
                -score              => $row[4],
                -rel_score          => $row[5],
                -seq                => $row[6],
                -conservation_level => $row[7],
                -conservation       => $row[8]
            )
        );
    }
    $sth->finish;
    $tf_site_set->param('gene_id', $gid);
    $tf_site_set->param('conservation_level', $clevel);

    return $tf_site_set;
}

=head2 fetch_by_gene_id

 Title    : fetch_by_gene_id
 Usage    : $tf_set = $ctfsa->fetch_by_gene_id($id, $level);
 Function : Fetch the set of conserved TF sites for the given gene
            ID at the given level of conservation.
            Synonym of fetch_set_by_gene_id
 Returns  : An OPOSSUM::ConservedTFBSSet object.
 Args	  : Integer gene ID and integer conservation level.

=cut

sub fetch_by_gene_id
{
    my ($self, $gene_id, $level) = @_;
    
    return $self->fetch_set_by_gene_id($gene_id, $level);
}

=head2 fetch_set_by_gene_tf_id

 Title    : fetch_set_by_gene_tf_id
 Usage    : $tf_set = $ctfsa->fetch_set_by_gene_tf_id();
 Function : Fetch the set of conserved TF sites for the given gene
            ID and TF ID at the given level of conservation.
 Returns  : An OPOSSUM::ConservedTFBSSet object.
 Args	  : Integer gene ID string TF ID and integer conservation level.

=cut

sub fetch_set_by_gene_tf_id
{
    my ($self, %args) = @_;
    
    return $self->fetch_set_by_gene_id(%args);
}

=head2 fetch_all_by_tf_id

 Title    : fetch_all_by_tf_id
 Usage    : $tfbss = $ctfbsa->fetch_all_by_tf_id(
                                 -tf_id                 => $tfid,
                                 -conservation_level    => $clevel,
                                 -threshold             => $threshold,
                                 -upstream_bp           => $upstream_bp,
                                 -downstream_bp         => $downstream_bp,
                                 -gene_ids         => gene_ids);
 Function : Fetch a list of conserved TFBSs for the given TF ID
            at the given conservation level. Optionally limit to only
            those sites which score greater than or equal to the given
            threshold and/or fall within the given bounds specified by
            the number of upstream and downstream bp and that provided
            list of gene IDs.
 Returns  : A list ref of OPOSSUM::ConservedTFBS objects.
 Args     : TF ID,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
            OPTIONAL gene ID list
=cut

sub fetch_all_by_tf_id
{
    my ($self, %args) = @_;

    my $tfid          = $args{-tf_id};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};
    my $gene_ids      = $args{-gene_ids};

    if (!$tfid) {
        carp "must provide TF ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $base_sql = qq{
        select tf_id, start, end, strand, score, rel_score, seq,
            conservation_level, conservation
        from conserved_tfbss 
        where tf_id = '$tfid'
            and conservation_level >= $clevel
            and gene_id = ?
    };

    if ($threshold) {
        if ($threshold =~ /(.+)%$/) {
            $threshold = $1 / 100;
        }
            
        $base_sql .= " and rel_score >= $threshold";
    }

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    $gene_ids = $ga->fetch_gene_ids() if !$gene_ids;

    my @tfbss;
    foreach my $gene_id (@$gene_ids) {
        my $sql = $base_sql;

        my $gene = $ga->fetch_by_gene_id($gene_id);
        if (!$gene) {
            carp "error fetching Gene by Gene ID $gene_id";
            next;
        }

        my $search_regions = $gene->promoter_search_regions(
            $upstream_bp, $downstream_bp
        );

        if ($search_regions) {
            my $first_sr = 1;
            foreach my $sr (@$search_regions) {
                my $sr_start = $sr->start;
                my $sr_end   = $sr->end;
                if ($first_sr) {
                    # TFBS must be completely contained in search region
                    # (not just overlapping)
                    $sql .= " and ((start >= $sr_start and end <= $sr_end)";
                    $first_sr = 0;
                } else {
                    $sql .= " or (start >= $sr_start and end <= $sr_end)";
                }
            }

            if (!$first_sr) {
                $sql .= ")";
            }
        }

        my $sth = $self->prepare($sql);
        if (!$sth) {
            carp "error preparing fetch conserved TF sites with:\n$sql\n"
                . $self->errstr;
            return;
        }

        if (!$sth->execute($gene_id)) {
            carp "error executing fetch conserved TF sites with:\n$sql\n"
                . $self->errstr;
            return;
        }

        while (my @row = $sth->fetchrow_array) {
            push @tfbss, OPOSSUM::ConservedTFBS->new(
                -gene_id            => $gene_id,
                -tf_id              => $row[0],
                -start              => $row[1],
                -end                => $row[2],
                -strand             => $row[3],
                -score              => $row[4],
                -rel_score          => $row[5],
                -seq                => $row[6],
                -conservation_level => $row[7],
                -conservation       => $row[8]
            );
        }
        $sth->finish;
    }

    return @tfbss ? \@tfbss : undef;
}

=head2 fetch_anchored_by_tf_id

 Title    : fetch_anchored_by_tf_id
 Usage    : $tfbss = $ctfbsa->fetch_anchored_by_tf_id(
                 -tf_id                 => $tfid,
                 -gene_id               => $gene_id,
                 -anchor_tf_id          => $anchor_tfid,
                 -distance              => $distance,
                 -conservation_level    => $clevel,
                 -threshold             => $threshold,
                 -upstream_bp           => $upstream_bp,
                 -downstream_bp         => $downstream_bp
            );
 Function : Fetch a list of conserved TFBSs for the given TF anchored
            by (proximal to) conserved TFBSs of the anchoring TF for the
            given gene at the given conservation level. Optionally limit
            to only those sites which score greater than or equal to the
            given threshold and/or fall within the given bounds specified
            by the number of upstream and downstream bp.
 Returns  : A list ref to a list ref of OPOSSUM::ConservedTFBS objects,
            and a list ref of hashes containing site pair info.
 Args     : TF ID,
            Gene ID,
            Anchor TF ID,
            Inter-binding distance,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
=cut

sub fetch_anchored_by_tf_id
{
    my ($self, %args) = @_;

    my $tfid          = $args{-tf_id};
    my $gene_id       = $args{-gene_id};
    my $anchor_tfid   = $args{-anchor_tf_id};
    my $distance      = $args{-distance};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$tfid) {
        carp "must provide TF ID";
        return;
    }

    if (!$gene_id) {
        carp "must provide gene ID";
        return;
    }

    if (!$anchor_tfid) {
        carp "must provide anchoring TF ID";
        return;
    }

    if (!$distance) {
        carp "must provide inter-binding distance";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    #if ($tfid eq $anchor_tfid) {
    #    carp "anchoring TF and search TF are the same";
    #    return;
    #}

    my $anchor_tfbss = $self->fetch_by_gene(
        -gene_id                => $gene_id,
        -tf_id                  => $anchor_tfid,
        -conservation_level     => $clevel,
        -threshold              => $threshold,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp
    );

    return if !$anchor_tfbss;

    my $tfbss = $self->fetch_by_gene(
        -gene_id                => $gene_id,
        -tf_id                  => $tfid,
        -conservation_level     => $clevel,
        -threshold              => $threshold,
        -upstream_bp            => $upstream_bp,
        -downstream_bp          => $downstream_bp
    );

    return if !$tfbss;

    my $sitepairs = _proximal_tfbss(
        $anchor_tfbss, $tfbss, $distance
    );

    return $sitepairs;
}

=head2 store

 Title   : store
 Usage   : $ctfsa->store($tfbs);
 Function: Store conserved TFBS in the database.
 Args    : An OPOSSUM::ConservedTFBS object
 Returns : True on success, false otherwise.

=cut

sub store
{
    my ($self, $tfbs) = @_;

    return if !$tfbs;

    if (!ref $tfbs || !$tfbs->isa("OPOSSUM::ConservedTFBS")) {
        carp "Not an OPOSSUM::ConservedTFBS object";
        return;
    }

    my $sql = qq{
        insert into conserved_tfbss (
		    gene_id,
		    tf_id,
		    start,
		    end,
		    strand,
		    score,
		    rel_score,
		    seq,
		    conservation_level,
		    conservation)
		values (?,?,?,?,?,?,?,?,?,?)
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert conserved TFBS statement\n"
            . $self->errstr;
        return;
    }

    if (
        !$sth->execute(
            $tfbs->gene_id,             $tfbs->tf_id,
            $tfbs->start,               $tfbs->end,
            $tfbs->strand,              $tfbs->score,
            $tfbs->rel_score,           $tfbs->seq,
            $tfbs->conservation_level,  $tfbs->conservation
        )
        )
    {
        carp "Error inserting conserved TFBS";
        return 0;
    }

    return 1;
}

=head2 store_list

 Title   : store_list
 Usage   : $ctfsa->store_list($tfbss);
 Function: Store conserved TFBSs in the database.
 Args    : A listref of OPOSSUM::ConservedTFBS objects
 Returns : True on success, false otherwise.

=cut

sub store_list
{
    my ($self, $tfbss) = @_;

    return if !$tfbss || !$tfbss->[0];

    my $sql = qq{insert into conserved_tfbss (
		    gene_id,
		    tf_id,
		    start,
		    end,
		    strand,
		    score,
		    rel_score,
		    seq,
		    conservation_level,
		    conservation)
		values (?,?,?,?,?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert conserved TFBS statement\n"
            . $self->errstr;
        return;
    }

    my $ok = 1;
    foreach my $tfbs (@$tfbss) {
        if (
            !$sth->execute(
                $tfbs->gene_id,            $tfbs->tf_id,
                $tfbs->start,              $tfbs->end,
                $tfbs->strand,             $tfbs->score,
                $tfbs->rel_score,          $tfbs->seq,                
                $tfbs->conservation_level, $tfbs->conservation
            )
            )
        {
            carp "Error inserting conserved TFBS";
            # keep trying to store TFBSs but return error status...
            $ok = 0;
        }
    }

    return $ok;
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

    my @sitepairs;
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
                push @sitepairs, {
                    anchor_site => $anchor,
                    tf_site     => $tfbs,
                    distance    => $dist
                };
            }
        }
    }

    return @sitepairs ? \@sitepairs : undef;
}

1;
