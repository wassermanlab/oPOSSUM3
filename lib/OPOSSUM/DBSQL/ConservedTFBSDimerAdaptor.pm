=head1 NAME

OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor - Adaptor for MySQL queries to
retrieve and store conserved TFBS dimers.

=head1 SYNOPSIS

$cra = $db_adaptor->get_ConservedTFBSDimerAdaptor();

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::ConservedTFBSDimer;
use OPOSSUM::ConservedTFBSSet;
use OPOSSUM::TFBSCount;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title    : new
 Usage    : $ctfsa = OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor->new(@args);
 Function : Create a new ConservedTFBSDimerAdaptor.
 Returns  : A new OPOSSUM::DBSQL::ConservedTFBSDimerAdaptor object.
 Args	  : An OPOSSUM::DBSQL::DBConnection object.

=cut

sub new
{
    my ($class, $dbobj, $ext) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new($dbobj);

    $self->{-extended_table_name} = $ext;

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

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{select distinct tf_id from $table};

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

=head2 fetch_gene_pair_ids

 Title    : fetch_gene_pair_ids
 Usage    : $ids = $ctfsa->fetch_gene_pair_ids();
 Function : Fetch a list of all the distinct gene pair IDs for all
	    the conserved TFBSs in the current DB.
 Returns  : A list ref of integer gene pair IDs.
 Args	  : None.

=cut

sub fetch_gene_pair_ids
{
    my ($self) = @_;

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{select distinct gene_pair_id from $table};

    my $sth = $self->prepare($sql);
    if (!$sth) {
	carp "error fetching gene pair IDs\n" . $self->errstr;
	return;
    }

    if (!$sth->execute) {
	carp "error fetching gene pair IDs\n" . $self->errstr;
	return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
	push @ids, $id;
    }

    return @ids ? \@ids : undef;
}

=head2 fetch_tfbs_count

 Title    : fetch_tfbs_count
 Usage    : $count = $ctfsa->fetch_tfbs_count(
				 -gene_pair_id		=> $gpid,
				 -tf_id			=> $tfid,
				 -conservation_level	=> $clevel,
				 -threshold		=> $threshold,
				 -upstream_bp		=> $upstream_bp,
				 -downstream_bp		=> $downstream_bp);
 Function : Fetch a count of the number of conserved TFBSs for the
	    given TF ID found for the given gene pair ID at the given
	    level of conservation, optional TFBS score threshold and
	    optional search region.
 Returns  : An integer count of the total number of conserved binding sites
	    of the given TF which fall within the promoter regions of the given
	    gene.
 Args	  : Gene pair ID,
 	    TF ID,
	    Conservation level,
	    OPTIONAL TFBS score threshold (range 0.0 - 1.0),
	    OPTIONAL upstream bp,
	    OPTIONAL downstream bp

=cut

sub fetch_tfbs_count
{
    my ($self, %args) = @_;

    my $gpid = $args{-gene_pair_id};
    my $tfid = $args{-tf_id};
    my $clevel = $args{-conservation_level};
    my $threshold = $args{-threshold};
    my $upstream_bp = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gpid) {
        carp "must provide TF ID";
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

    my $gpa = $self->db->get_GenePairAdaptor;
    if (!$gpa) {
    	carp "error getting GenePairAdaptor";
	return;
    }

    my $ppa = $self->db->get_PromoterPairAdaptor;
    if (!$ppa) {
    	carp "error getting PromoterPairAdaptor";
	return;
    }

    my $gene_pair = $gpa->fetch_by_gene_pair_id($gpid);
    if (!$gene_pair) {
    	carp "error fetching GenePair $gpid";
	return;
    }

    my $promoter_pairs = $ppa->fetch_by_gene_pair_id($gpid);
    if (!$promoter_pairs) {
    	carp "error fetching PromoterPairs for GenePair ID $gpid";
	return;
    }

    my $strand = $gene_pair->strand1;

    #
    # We want to find all the TFBSs within all the promoters for the given gene
    # (optionally limited to specified number of upstream and downstream bp of
    # the TSS). These promoter regions may overlap and we do not want to double
    # count so combine overlapping promoter regions into single search regions.
    #
    my @search_regions = _define_search_regions($promoter_pairs, $upstream_bp,
    				$downstream_bp, $strand);
    @search_regions = _combine_search_regions(@search_regions);

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{
		    select count(tf_id) from $table
		    where gene_pair_id = $gpid and tf_id = $tfid
		    and conservation_level >= $clevel
		    and start1 >= ? and end1 <= ?
		};

    if ($threshold) {
    	$threshold /= 100 if $threshold > 1;
    	$sql .= " and rel_score1 >= $threshold and rel_score2 >= $threshold";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
	carp "error preparing fetch conserved TF site count with:\n$sql\n"
		. $self->errstr;
	return;
    }

    my $count = 0;
    foreach my $sr (@search_regions) {
	my $sr_start = $sr->start;
	my $sr_end = $sr->end;

	#
	# TFBS must be completely contained in search region
	# (not just overlapping)
	#
	if (!$sth->execute($sr_start, $sr_end)) {
	    carp "error executing fetch conserved TF site count with:\n$sql\n"
		    . $self->errstr;
	    return;
	}

	if (my @row = $sth->fetchrow_array) {
	    $count += $row[0];
	}
    }

    return $count;
}

=head2 fetch_tfbs_counts

 Title    : fetch_tfbs_counts
 Usage    : $counts = $ctfsa->fetch_tfbs_counts(
				 -gene_pair_id		=> $gpid,
				 -tf_ids		=> $tfids,
				 -conservation_level	=> $clevel,
				 -threshold		=> $threshold,
				 -upstream_bp		=> $upstream_bp,
				 -downstream_bp		=> $downstream_bp);
 Function : Fetch counts of the number of conserved binding sites for the
	    given TFs (all TFs if not specified) found for the given gene pair
	    ID at the given level of conservation, optional TFBS score
	    threshold and optional search region.
 Returns  : A listref of OPOSSUM::TFBSCount objects.
 Args	  : Gene pair ID,
	    Conservation level,
 	    OPTIONAL TF IDs,
	    OPTIONAL TFBS score threshold (range 0.0 - 1.0),
	    OPTIONAL upstream bp,
	    OPTIONAL downstream bp

=cut

sub fetch_tfbs_counts
{
    my ($self, %args) = @_;

    my $gpid = $args{-gene_pair_id};
    my $tfids = $args{-tf_ids};
    my $clevel = $args{-conservation_level};
    my $threshold = $args{-threshold};
    my $upstream_bp = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gpid) {
        carp "must provide gene pair ID";
	return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
	return;
    }

    my $gpa = $self->db->get_GenePairAdaptor;
    if (!$gpa) {
    	carp "error getting GenePairAdaptor";
	return;
    }

    my $ppa = $self->db->get_PromoterPairAdaptor;
    if (!$ppa) {
    	carp "error getting PromoterPairAdaptor";
	return;
    }

    my $gene_pair = $gpa->fetch_by_gene_pair_id($gpid);
    if (!$gene_pair) {
    	carp "error fetching GenePair $gpid";
	return;
    }

    my $promoter_pairs = $ppa->fetch_by_gene_pair_id($gpid);
    if (!$promoter_pairs) {
    	carp "error fetching PromoterPairs for GenePair ID $gpid";
	return;
    }

    my $strand = $gene_pair->strand1;

    #
    # We want to find all the TFBSs within all the promoters for the given gene
    # (optionally limited to specified number of upstream and downstream bp of
    # the TSS). These promoter regions may overlap and we do not want to double
    # count so combine overlapping promoter regions into single search regions.
    #
    my @search_regions = _define_search_regions($promoter_pairs, $upstream_bp,
    				$downstream_bp, $strand);
    @search_regions = _combine_search_regions(@search_regions);

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{
		    select tf_id, count(tf_id) from $table
		    where gene_pair_id = $gpid
		    and conservation_level >= $clevel
		};

    if ($threshold) {
    	$threshold /= 100 if $threshold > 1;
    	$sql .= " and rel_score1 >= $threshold and rel_score2 >= $threshold";
    }

    my $first_sr = 1;
    foreach my $sr (@search_regions) {
        my $sr_start = $sr->start;
        my $sr_end = $sr->end;
        if ($first_sr) {
            # TFBS must be completely contained in search region
            # (not just overlapping)
            $sql .= " and ((start1 >= $sr_start and end1 <= $sr_end)";
            $first_sr = 0;
        } else {
            $sql .= " or (start1 >= $sr_start and end1 <= $sr_end)";
        }
    }

    if (!$first_sr) {
        $sql .= ")";
    }

    if ($tfids && $tfids->[0]) {
        $sql .= " and tf_id in (";
	$sql .= join(',', @$tfids);
        $sql .= ")";
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
			    -gene_pair_id	=> $gpid,
			    -conservation_level	=> $clevel,
			    -tf_id		=> $row[0],
			    -count		=> $row[1]);
    }

    return @counts ? \@counts : undef;
}

###  fetch_by_gene_pair added by SHS on Dec 27,2006, modified on Jan 10, 2007

=head2 fetch_by_gene_pair

 Title    : fetch_by_gene_pair
 Usage    : $count = $ctfbsa->fetch_by_gene_pair(
                                 -gene_pair_id          => $gpid,
                                 -tf_id                 => $tfid,
                                 -conservation_level    => $clevel,
                                 -threshold             => $threshold,
                                 -upstream_bp           => $upstream_bp,
				 -downstream_bp         => $downstream_bp);
 Function : Fetch a list of conserved TFBSs for the given TF ID
            found for the given gene pair ID at the given conservation
            level. Optionally limit to only those sites which score
            greater than or equal to the given threshold and/or fall
            within the given bounds specified by the number of upstream
            and downstream bp.
 Returns  : A list ref of OPOSSUM::ConservedTFBS objects.
 Args     : Gene pair ID,
            TF ID,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
=cut

sub fetch_by_gene_pair
{
    my ($self, %args) = @_;

    my $gpid = $args{-gene_pair_id};
    my $tfid = $args{-tf_id};
    my $clevel = $args{-conservation_level};
    my $threshold = $args{-threshold};
    my $upstream_bp = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gpid) {
        carp "must provide gene pair ID";
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

    my $gpa = $self->db->get_GenePairAdaptor;
    if (!$gpa) {
        carp "error getting GenePairAdaptor";
        return;
    }

    my $ppa = $self->db->get_PromoterPairAdaptor;
    if (!$ppa) {
        carp "error getting PromoterPairAdaptor";
        return;
    }

    my $gene_pair = $gpa->fetch_by_gene_pair_id($gpid);
    if (!$gene_pair) {
        carp "error fetching GenePair $gpid";
        return;
    }

    my $promoter_pairs = $ppa->fetch_by_gene_pair_id($gpid);
    if (!$promoter_pairs) {
        carp "error fetching PromoterPairs for GenePair ID $gpid";
        return;
    }

    my $strand1 = $gene_pair->strand1;
    my $strand2 = $gene_pair->strand2;

    # 
    # We want to find all the TFBSs within all the promoters for the given gene
    # (optionally limited to specified number of upstream and downstream bp of   
    # the TSS). These promoter regions may overlap and we do not want to double
    # count so combine overlapping promoter regions into single search regions.
    #                                                   
    my @tf_sites;

    my @search_regions;
    foreach my $pp (@$promoter_pairs) {
        my $pp_tss1 = $pp->tss1;
        my $pp_start1 = $pp->start1;
        my $pp_end1 = $pp->end1;

        my ($tss_start, $tss_end) = _define_search_region( $pp_tss1, $pp_start1,
							   $pp_end1, $strand1, $upstream_bp,
							   $downstream_bp);

	push @search_regions, Bio::SeqFeature::Generic->new(
							    -start  => $tss_start,
							    -end    => $tss_end);
    }
    
    @search_regions = _combine_search_regions(@search_regions);

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{select tf_id, start1, end1, rel_start1, rel_end1, strand1,
                   score1, rel_score1, seq1, start2, end2, rel_start2,
                   rel_end2, strand2, score2, rel_score2, seq2,
                   conservation_level, conservation, repeat_type1, repeat_type2
                 from $table 
                 where gene_pair_id = $gpid
                   and conservation_level >= $clevel
                   and tf_id = $tfid};

    if ($threshold) {
        $threshold /= 100 if $threshold > 1;
        $sql .= " and rel_score1 >= $threshold and rel_score2 >= $threshold";
    }

    my $first_sr = 1;
    foreach my $sr (@search_regions) {
	my $sr_start = $sr->start;
	my $sr_end = $sr->end;
	if ($first_sr) {
            # TFBS must be completely contained in search region
            # (not just overlapping)
	    $sql .= " and ((start1 >= $sr_start and end1 <= $sr_end)";
	    $first_sr = 0;
	} else {
	    $sql .= " or (start1 >= $sr_start and end1 <= $sr_end)";
	}
    }

    if (!$first_sr) {
	$sql .= ")";
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
	
    while (my @row = $sth->fetchrow_array) {
	push @tf_sites, OPOSSUM::ConservedTFBSDimer->new(
				    -id             => $row[0],
				    -start1         => $row[1],
				    -end1           => $row[2],
				    -rel_start1		=> $row[3],
				    -rel_end1		=> $row[4],
				    -strand1        => $row[5],
                    -score1         => $row[6],
				    -rel_score1		=> $row[7],
                    -seq1           => $row[8],
                    -start2         => $row[9],
                    -end2           => $row[10],
				    -rel_start2		=> $row[11],
				    -rel_end2		=> $row[12],
                    -strand2        => $row[13],
                    -score2         => $row[14],
				    -rel_score2		=> $row[15],
                    -seq2           => $row[16],
                    -conservation_level => $row[17],
				    -conservation   => $row[18],
                    -repeat_type1   => $row[19],
                    -repeat_type2   => $row[20]);
    }
    $sth->finish;

    return @tf_sites ? \@tf_sites : undef;
}

=head2 fetch_all_by_gene_pair_id

 Title    : fetch_all_by_gene_pair_id
 Usage    : $tfbss = $ctfbsa->fetch_all_by_gene_pair_id(
                                 -gene_pair_id          => $gpid,
                                 -conservation_level    => $clevel,
                                 -threshold             => $threshold,
                                 -upstream_bp           => $upstream_bp,
                                 -downstream_bp         => $downstream_bp);
 Function : Fetch a list of all conserved TFBSs for the given gene pair
            ID at the given conservation level.
            Optionally limit to only those sites which score
            greater than or equal to the given threshold and/or fall
            within the given bounds specified by the number of upstream
            and downstream bp.
 Returns  : A list ref of OPOSSUM::ConservedTFBSDimer objects.
 Args     : Gene pair ID,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
=cut

sub fetch_all_by_gene_pair_id
{
    my ($self, %args) = @_;

    my $gpid = $args{-gene_pair_id};
    my $clevel = $args{-conservation_level};
    my $threshold = $args{-threshold};
    my $upstream_bp = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};

    if (!$gpid) {
        carp "must provide gene pair ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $gpa = $self->db->get_GenePairAdaptor;
    if (!$gpa) {
        carp "error getting GenePairAdaptor";
        return;
    }

    my $ppa = $self->db->get_PromoterPairAdaptor;
    if (!$ppa) {
        carp "error getting PromoterPairAdaptor";
        return;
    }

    my $gene_pair = $gpa->fetch_by_gene_pair_id($gpid);
    if (!$gene_pair) {
        carp "error fetching GenePair $gpid";
        return;
    }

    my $promoter_pairs = $ppa->fetch_by_gene_pair_id($gpid);
    if (!$promoter_pairs) {
        carp "error fetching PromoterPairs for GenePair ID $gpid";
        return;
    }

    my $strand1 = $gene_pair->strand1;
    my $strand2 = $gene_pair->strand2;

    # 
    # We want to find all the TFBSs within all the promoters for the given gene
    # (optionally limited to specified number of upstream and downstream bp of   
    # the TSS). These promoter regions may overlap and we do not want to double
    # count so combine overlapping promoter regions into single search regions.
    #                                                   
    my @tf_sites;

    my @search_regions;
    foreach my $pp (@$promoter_pairs) {
        my $pp_tss1 = $pp->tss1;
        my $pp_start1 = $pp->start1;
        my $pp_end1 = $pp->end1;

        my ($tss_start, $tss_end) = _define_search_region( $pp_tss1, $pp_start1,
							   $pp_end1, $strand1, $upstream_bp,
							   $downstream_bp);

	push @search_regions, Bio::SeqFeature::Generic->new(
							    -start  => $tss_start,
							    -end    => $tss_end);
    }
    
    @search_regions = _combine_search_regions(@search_regions);

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{select tf_id, start1, end1, rel_start1, rel_end1, strand1,
		    score1, rel_score1, seq1, start2, end2, rel_start2,
		    rel_end2, strand2, score2, rel_score2, seq2,
		    conservation_level, conservation, repeat_type1, repeat_type2
		from $table 
		where gene_pair_id = $gpid
		    and conservation_level >= $clevel
		};

    if ($threshold) {
        $threshold /= 100 if $threshold > 1;
        $sql .= " and rel_score1 >= $threshold and rel_score2 >= $threshold";
    }

    my $first_sr = 1;
    foreach my $sr (@search_regions) {
        my $sr_start = $sr->start;
        my $sr_end = $sr->end;
        if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
            $sql .= " and ((start1 >= $sr_start and end1 <= $sr_end)";
            $first_sr = 0;
        } else {
            $sql .= " or (start1 >= $sr_start and end1 <= $sr_end)";
        }
    }

    if (!$first_sr) {
        $sql .= ")";
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
	
    while (my @row = $sth->fetchrow_array) {
        push @tf_sites, OPOSSUM::ConservedTFBSDimer->new(
				    -id                 => $row[0],
				    -start1             => $row[1],
				    -end1               => $row[2],
				    -rel_start1		    => $row[3],
				    -rel_end1		    => $row[4],
				    -strand1            => $row[5],
                    -score1             => $row[6],
				    -rel_score1		    => $row[7],
                    -seq1               => $row[8],
                    -start2             => $row[9],
                    -end2               => $row[10],
				    -rel_start2		    => $row[11],
				    -rel_end2		    => $row[12],
                    -strand2            => $row[13],
                    -score2             => $row[14],
				    -rel_score2		    => $row[15],
                    -seq2               => $row[16],
                    -conservation_level => $row[17],
				    -conservation       => $row[18],
                    -repeat_type1   => $row[19],
                    -repeat_type2   => $row[20]);
    }
    $sth->finish;

    return @tf_sites ? \@tf_sites : undef;
}

=head2 fetch_by_gene_pair_id

 Title    : fetch_by_gene_pair_id
 Usage    : $tf_set = $ctfsa->fetch_by_gene_pair_id($id, $level);
 Function : Fetch the set of conserved TF sites for the given gene
 	    pair ID at the given level of conservation.
 Returns  : An OPOSSUM::ConservedTFBSSet object.
 Args	  : Integer gene pair ID and integer conservation level.

=cut

sub fetch_by_gene_pair_id
{
    my ($self, $gene_pair_id, $level) = @_;

    return if !$gene_pair_id;

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{select tf_id, start1, end1, rel_start1, rel_end1, strand1,
		    score1, rel_score1, seq1, start2, end2, rel_start2,
		    rel_end2, strand2, score2, rel_score2, seq2,
		    conservation_level, conservation, repeat_type1, repeat_type2
		from $table
		where gene_pair_id = $gene_pair_id
	    };

    if ($level) {
	$sql .= " and level >= $level";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
	carp "error fetching conserved TF site with gene_pair_id ="
		. " $gene_pair_id\n" . $self->errstr;
	return;
    }

    if (!$sth->execute) {
	carp "error fetching conserved TF site with gene_pair_id ="
		. " $gene_pair_id\n" . $self->errstr;
	return;
    }

    my $tf_site_set = OPOSSUM::ConservedTFBSSet->new();
    while (my @row = $sth->fetchrow_array) {
	$tf_site_set->add_tf_site(
			OPOSSUM::ConservedTFBSDimer->new(
				    -id                 => $row[0],
				    -start1             => $row[1],
				    -end1               => $row[2],
				    -rel_start1		=> $row[3],
				    -rel_end1		=> $row[4],
				    -strand1            => $row[5],
                                    -score1             => $row[6],
				    -rel_score1		=> $row[7],
                                    -seq1               => $row[8],
                                    -start2             => $row[9],
                                    -end2               => $row[10],
				    -rel_start2		=> $row[11],
				    -rel_end2		=> $row[12],
                                    -strand2            => $row[13],
                                    -score2             => $row[14],
				    -rel_score2		=> $row[15],
                                    -seq2               => $row[16],
                                    -conservation_level => $row[17],
				    -conservation       => $row[18],
                    -repeat_type1   => $row[19],
                    -repeat_type2   => $row[20]));
    }
    $sth->finish;
    $tf_site_set->param('gene_pair_id', $gene_pair_id);
    $tf_site_set->param('conservation_level', $level);

    return $tf_site_set;
}

=head2 fetch_all_by_tf_id

 Title    : fetch_all_by_tf_id
 Usage    : $count = $ctfbsa->fetch_all_by_tf_id(
                                 -tf_id                 => $tfid,
                                 -conservation_level    => $clevel,
                                 -threshold             => $threshold,
                                 -upstream_bp           => $upstream_bp,
                                 -downstream_bp         => $downstream_bp,
                                 -gene_pair_ids         => gene_pair_ids);
 Function : Fetch a list of conserved TFBSs for the given TF ID
            at the given conservation level. Optionally limit to only
            those sites which score greater than or equal to the given
            threshold and/or fall within the given bounds specified by
            the number of upstream and downstream bp and that provided
            list of gene pair IDs.
 Returns  : A list ref of OPOSSUM::ConservedTFBSDimer objects.
 Args     : TF ID,
            Conservation level,
            OPTIONAL TFBS score threshold (range 0.0 - 1.0),
            OPTIONAL upstream bp,
            OPTIONAL downstream bp
            OPTIONAL gene pair ID list
=cut

sub fetch_all_by_tf_id
{
    my ($self, %args) = @_;

    my $tfid          = $args{-tf_id};
    my $clevel        = $args{-conservation_level};
    my $threshold     = $args{-threshold};
    my $upstream_bp   = $args{-upstream_bp};
    my $downstream_bp = $args{-downstream_bp};
    my $gene_pair_ids = $args{-gene_pair_ids};

    if (!$tfid) {
        carp "must provide TF ID";
        return;
    }

    if (!$clevel) {
        carp "must provide a conservation level";
        return;
    }

    my $gpa = $self->db->get_GenePairAdaptor;
    if (!$gpa) {
        carp "error getting GenePairAdaptor";
        return;
    }

    my $ppa = $self->db->get_PromoterPairAdaptor;
    if (!$ppa) {
        carp "error getting PromoterPairAdaptor";
        return;
    }

    $gene_pair_ids = $gpa->fetch_gene_pair_ids() if !$gene_pair_ids;

    my @tf_sites;
    foreach my $gpid (@$gene_pair_ids) {
        my $gene_pair = $gpa->fetch_by_gene_pair_id($gpid);
        if (!$gene_pair) {
            carp "error fetching GenePair by GenePair ID $gpid";
            next;
        }

        my $promoter_pairs = $ppa->fetch_by_gene_pair_id($gpid);
        if (!$promoter_pairs) {
            carp "error fetching PromoterPairs for GenePair ID $gpid";
            next;
        }

        my $strand1 = $gene_pair->strand1;
        my $strand2 = $gene_pair->strand2;

        # 
        # We want to find all the TFBSs within all the promoters for the given
        # gene (optionally limited to specified number of upstream and
        # downstream bp of the TSS). These promoter regions may overlap and
        # we do not want to double count so combine overlapping promoter
        # regions into single search regions.
        #                                                   

        my @search_regions;
        foreach my $pp (@$promoter_pairs) {
            my $pp_tss1 = $pp->tss1;
            my $pp_start1 = $pp->start1;
            my $pp_end1 = $pp->end1;

            my ($tss_start, $tss_end) = _define_search_region(
                                            $pp_tss1, $pp_start1,
                                            $pp_end1, $strand1, $upstream_bp,
                                            $downstream_bp);

            push @search_regions, Bio::SeqFeature::Generic->new(
                                    -start  => $tss_start,
                                    -end    => $tss_end);
        }
    
        @search_regions = _combine_search_regions(@search_regions);

        my $table = "conserved_tfbss";
        $table .= "_" . $self->{-extended_table_name}
                    if $self->{-extended_table_name};

        my $sql = qq{select gene_pair_id, tf_id, start1, end1, rel_start1,
                        rel_end1, strand1, score1, rel_score1, seq1, start2,
                        end2, rel_start2, rel_end2, strand2, score2,
                        rel_score2, seq2, conservation_level, conservation,
                        repeat_type1, repeat_type2
                     from $table 
                     where gene_pair_id = $gpid
                        and conservation_level >= $clevel
                        and tf_id = $tfid};

        if ($threshold) {
            $threshold /= 100 if $threshold > 1;
            $sql .= " and rel_score1 >= $threshold
                and rel_score2 >= $threshold";
        }

        my $first_sr = 1;
        foreach my $sr (@search_regions) {
            my $sr_start = $sr->start;
            my $sr_end = $sr->end;
            if ($first_sr) {
                # TFBS must be completely contained in search region
                # (not just overlapping)
                $sql .= " and ((start1 >= $sr_start and end1 <= $sr_end)";
                $first_sr = 0;
            } else {
                $sql .= " or (start1 >= $sr_start and end1 <= $sr_end)";
            }
        }

        if (!$first_sr) {
            $sql .= ")";
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
	
        while (my @row = $sth->fetchrow_array) {
            push @tf_sites, OPOSSUM::ConservedTFBSDimer->new(
				    -gene_pair_id       => $row[0],
				    -id                 => $row[1],
				    -start1             => $row[2],
				    -end1               => $row[3],
				    -rel_start1		    => $row[4],
				    -rel_end1		    => $row[5],
				    -strand1            => $row[6],
                    -score1             => $row[7],
				    -rel_score1		    => $row[8],
                    -seq1               => $row[9],
                    -start2             => $row[10],
                    -end2               => $row[11],
				    -rel_start2		    => $row[12],
				    -rel_end2		    => $row[13],
                    -strand2            => $row[14],
                    -score2             => $row[15],
				    -rel_score2		    => $row[16],
                    -seq2               => $row[17],
                    -conservation_level => $row[18],
				    -conservation       => $row[19],
                    -repeat_type1       => $row[20],
                    -repeat_type2       => $row[21]);
        }
        $sth->finish;
    }

    return @tf_sites ? \@tf_sites : undef;
}

=head2 store

 Title   : store
 Usage   : $ctfsa->store($tfbs);
 Function: Store conserved TFBS in the database.
 Args    : An OPOSSUM::ConservedTFBSDimer object
 Returns : True on success, false otherwise.

=cut

sub store
{
    my ($self, $tfbs) = @_;

    return if !$tfbs;

    if (!ref $tfbs || !$tfbs->isa("OPOSSUM::ConservedTFBSDimer")) {
        carp "not an OPOSSUM::ConservedTFBSDimer object";
        return;
    }

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{insert into $table (
		    gene_pair_id,
		    tf_id,
		    start1,
		    end1,
		    rel_start1,
		    rel_end1,
		    strand1,
		    score1,
		    rel_score1,
		    seq1,
		    start2,
		    end2,
		    rel_start2,
		    rel_end2,
		    strand2,
		    score2,
		    rel_score2,
		    seq2,
		    conservation_level,
		    conservation,
            repeat_type1,
            repeat_type2)
		values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
    	carp "Error preparing insert conserved TFBS statement\n"
		. $self->errstr;
        return;
    }

    if (!$sth->execute(
			$tfbs->gene_pair_id,
			$tfbs->tf_id,
			$tfbs->start1,
			$tfbs->end1,
			$tfbs->rel_start1,
			$tfbs->rel_end1,
			$tfbs->strand1,
			$tfbs->score1,
			$tfbs->rel_score1,
			$tfbs->seq1,
			$tfbs->start2,
			$tfbs->end2,
			$tfbs->rel_start2,
			$tfbs->rel_end2,
			$tfbs->strand2,
			$tfbs->score2,
			$tfbs->rel_score2,
			$tfbs->seq2,
			$tfbs->conservation_level,
			$tfbs->conservation,
            $tfbs->repeat_type1,
            $tfbs->repeat_type2))
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
 Args    : A listref of OPOSSUM::ConservedTFBSDimer objects
 Returns : True on success, false otherwise.

=cut

sub store_list
{
    my ($self, $tfbss) = @_;

    return if !$tfbss || !$tfbss->[0];

    my $table = "conserved_tfbss";
    $table .= "_" . $self->{-extended_table_name}
    			if $self->{-extended_table_name};
    my $sql = qq{insert into $table (
		    gene_pair_id,
		    tf_id,
		    start1,
		    end1,
		    rel_start1,
		    rel_end1,
		    strand1,
		    score1,
		    rel_score1,
		    seq1,
		    start2,
		    end2,
		    rel_start2,
		    rel_end2,
		    strand2,
		    score2,
		    rel_score2,
		    seq2,
		    conservation_level,
		    conservation,
            repeat_type1,
            repeat_type2)
		values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
    	carp "Error preparing insert conserved TFBS statement\n"
		. $self->errstr;
        return;
    }

    my $ok = 1;
    foreach my $tfbs (@$tfbss) {
        if (!$sth->execute(
                    $tfbs->gene_pair_id,
                    $tfbs->tf_id,
                    $tfbs->start1,
                    $tfbs->end1,
                    $tfbs->rel_start1,
                    $tfbs->rel_end1,
                    $tfbs->strand1,
                    $tfbs->score1,
                    $tfbs->rel_score1,
                    $tfbs->seq1,
                    $tfbs->start2,
                    $tfbs->end2,
                    $tfbs->rel_start2,
                    $tfbs->rel_end2,
                    $tfbs->strand2,
                    $tfbs->score2,
                    $tfbs->rel_score2,
                    $tfbs->seq2,
                    $tfbs->conservation_level,
                    $tfbs->conservation,
                    $tfbs->repeat_type1,
                    $tfbs->repeat_type2))
        {
            carp "Error inserting conserved TFBS";
            # keep trying to store TFBSs but return error status...
            $ok = 0;
        }
    }

    return $ok;
}

sub _define_search_region 
{
    my ($pp_tss, $pp_start, $pp_end, $strand, $upstream_bp,
    	$downstream_bp) = @_;

    my $tss_start;
    my $tss_end;
    if ($strand == 1) {
	if (defined $upstream_bp) {
	    $tss_start = $pp_tss - $upstream_bp;
	    $tss_start = $pp_start if $pp_start > $tss_start;
	} else {
	    $tss_start = $pp_start;
	}
	
	if (defined $downstream_bp) {
	    $tss_end = $pp_tss + $downstream_bp - 1;
	    $tss_end = $pp_end if $pp_end < $tss_end;
	} else {
	    $tss_end = $pp_end;
	}
    } elsif ($strand == -1) {
	if (defined $upstream_bp) {
	    $tss_end = $pp_tss + $upstream_bp;
	    $tss_end = $pp_end if $pp_end < $tss_end;
	} else {
	    $tss_end = $pp_end;
	}
	
	if (defined $downstream_bp) {
	    $tss_start = $pp_tss - $downstream_bp + 1;
	    $tss_start = $pp_start if $pp_start > $tss_start;
	} else {
	    $tss_start = $pp_start;
	}
    } else {
	carp "error determining gene pair strand";
	return;
    }
    return ($tss_start, $tss_end);
}

sub _combine_search_region_pairs 
{
    my (@regs) = @_;
    
    @regs = sort {$a->feature1->start <=> $b->feature1->start} @regs;

    my $num_regs = scalar @regs;
    for (my $i = 0; $i < $num_regs; $i++) {
        my $reg1 = $regs[$i] if exists($regs[$i]);
        if ($reg1) {
            for (my $j = $i+1; $j < $num_regs; $j++) {
                my $reg2 = $regs[$j] if exists ($regs[$j]);
                if ($reg2) {
                    if (_feature_combine($reg1->feature1, $reg2->feature1)) {
                        if ($reg2->feature1->start < $reg1->feature1->start) {
                            $reg1->feature1->start($reg2->feature1->start);
			    if($reg2->feature2->start < $reg1->feature2->start){
				$reg1->feature2->start($reg2->feature2->start);
			    }
                        }

                        if ($reg2->feature1->end > $reg1->feature1->end) {
                            $reg1->feature1->end($reg2->feature1->end);
			    if ($reg2->feature2->end > $reg1->feature2->end) {
				$reg1->feature2->end($reg2->feature2->end);
			    }			    
                        }
                        delete $regs[$j];
                    } else {
                        last;
                    }
                }
            }
        }
    }
    my @unique_regs;
    foreach my $reg (@regs) {
        if (defined $reg) {
            push @unique_regs, $reg;
        }
    }

    return @unique_regs;
}

sub _define_search_regions
{
    my ($promoter_pairs, $upstream_bp, $downstream_bp, $strand) = @_;

    my @search_regions;
    foreach my $pp (@$promoter_pairs) {
	my $pp_tss1 = $pp->tss1;
	my $pp_start1 = $pp->start1;
	my $pp_end1 = $pp->end1;

	my $tss_start;
	my $tss_end;
	if ($strand == 1) {
	    if (defined $upstream_bp) {
		$tss_start = $pp_tss1 - $upstream_bp;
		$tss_start = $pp_start1 if $pp_start1 > $tss_start;
	    } else {
		$tss_start = $pp_start1;
	    }

	    if (defined $downstream_bp) {
		$tss_end = $pp_tss1 + $downstream_bp - 1;
		$tss_end = $pp_end1 if $pp_end1 < $tss_end;
	    } else {
		$tss_end = $pp_end1;
	    }
	} elsif ($strand == -1) {
	    if (defined $upstream_bp) {
		$tss_end = $pp_tss1 + $upstream_bp;
		$tss_end = $pp_end1 if $pp_end1 < $tss_end;
	    } else {
		$tss_end = $pp_end1;
	    }

	    if (defined $downstream_bp) {
		$tss_start = $pp_tss1 - $downstream_bp + 1;
		$tss_start = $pp_start1 if $pp_start1 > $tss_start;
	    } else {
		$tss_start = $pp_start1;
	    }
	} else {
	    carp "error determining promoter pair strand";
	    return;
	}

	push @search_regions, Bio::SeqFeature::Generic->new(
					-start	=> $tss_start,
					-end	=> $tss_end);
    }

    return @search_regions;
}

sub _combine_search_regions
{
    my (@regs) = @_;

    @regs = sort {$a->start <=> $b->start} @regs;

    my $num_regs = scalar @regs;
    for (my $i = 0; $i < $num_regs; $i++) {
	my $reg1 = $regs[$i] if exists($regs[$i]);
	if ($reg1) {
	    for (my $j = $i+1; $j < $num_regs; $j++) {
		my $reg2 = $regs[$j] if exists ($regs[$j]);
		if ($reg2) {
		    if (_feature_combine($reg1, $reg2)) {
			if ($reg2->start < $reg1->start) {
			    $reg1->start($reg2->start);
			}

			if ($reg2->end > $reg1->end) {
			    $reg1->end($reg2->end);
			}
			delete $regs[$j];
		    } else {
			last;
		    }
		}
	    }
	}
    }

    my @unique_regs;
    foreach my $reg (@regs) {
	if (defined $reg) {
	    push @unique_regs, $reg;
	}
    }

    return @unique_regs;
}

sub _feature_combine
{
    my ($feat1, $feat2) = @_;

    my $combine = 1;
    $combine = 0 if $feat1->start > $feat2->end + 1
			|| $feat1->end < $feat2->start - 1;

    return $combine;
}

1;
