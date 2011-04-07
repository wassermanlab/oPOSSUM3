
=head1 NAME

OPOSSUM::DBSQL::ConservedRegionAdaptor - Adaptor for MySQL queries to retrieve
and store conserved regions.

=head1 SYNOPSIS

$cra = $db_adaptor->get_ConservedRegionAdaptor();

=head1 CHANGE HISTORY

 AK 2010/09/03
 - Modified to reflect the addition of the gc_content column to the schema
 - added fetch_gc_content_by_gene_id and fetch_gc_content_by_upstream_downstream
 
 DJA 2010/02/16
 - Modified to reflect new oPOSSUM_2010 schema
 - The functionality provided by _define_search_regions has now been
   implemented in the OPOSSUM::Gene::promoter_search_regions() method.
   e.g. OPOSSUM::Gene->promoter_search_regions($upstream, $downstream);

 DJA 2007/01/29
 - modified fetch_set_by_upstream_downstream and
   fetch_list_by_upstream_downstream methods to fetch and truncate
   boundaries of conserved regions on a per search region basis.
 - added new method fetch_length_by_upstream_downstream
 - created new internal routine _define_search_regions

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ConservedRegionAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::ConservedRegion;
use OPOSSUM::ConservedRegionSet;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title    : new
 Usage    : $cra = OPOSSUM::DBSQL::ConservedRegionAdaptor->new(@args);
 Function : Create a new ConservedRegionAdaptor.
 Returns  : A new OPOSSUM::DBSQL::ConservedRegionAdaptor object.
 Args	  : An OPOSSUM::DBSQL::DBConnection object.

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_by_gene_id

 Title    : fetch_by_gene_id
 Usage    : $regions = $cra->fetch_by_gene_id($gid, $cons_level);
 Function : Alias for fetch_list_by_gene_id.
 Returns  : A listref of OPOSSUM::ConservedRegion objects.
 Args	  : Gene ID,
            Conservation level

=cut

sub fetch_by_gene_id
{
    my ($self, $gid, $cons_level) = @_;

    return $self->fetch_list_by_gene_id($gid, $cons_level);
}

=head2 fetch_list_by_gene_id

 Title    : fetch_list_by_gene_id
 Usage    : $regions = $cra->fetch_list_by_gene_id($gid, $cons_level);
 Function : Fetch the list of conserved regions for a given gene
            at the given level of conservation from the DB.
 Returns  : A listref of OPOSSUM::ConservedRegion objects.
 Args	  : Gene pair ID,
            Conservation level

=cut

sub fetch_list_by_gene_id
{
    my ($self, $gid, $cons_level) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $sql = qq{select start, end, conservation_level, conservation, gc_content
		from conserved_regions where gene_id = $gid
		and conservation_level = $cons_level order by start};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        my $error = "ERROR preparing fetch conserved region for gene $gid"
            . " at conservation level $cons_level";
        carp $error . "\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        my $error = "ERROR executing fetch conserved region for Gene $gid"
            . " at conservation level $cons_level";
        carp $error . "\n" . $self->errstr;
        return;
    }

    my @regions;
    while (my @row = $sth->fetchrow_array) {
        push @regions, OPOSSUM::ConservedRegion->new(
            -adaptor            => $self,
            -gene_id            => $gid,
            -start              => $row[0],
            -end                => $row[1],
            -conservation_level => $row[2],
            -conservation       => $row[3],
			-gc_content			=> $row[4]
        );
    }
    $sth->finish;

    return @regions ? \@regions : undef;
}

=head2 fetch_set_by_gene_id

 Title    : fetch_set_by_gene_id
 Usage    : $cr_set = $cra->fetch_set_by_gene_id($gid, $cons_level);
 Function : Fetch the set of conserved regions for a given gene
            at the given level of conservation from the DB.
 Returns  : An OPOSSUM::ConservedRegionSet object.
 Args	  : Integer gene ID and integer conservation level.

=cut

sub fetch_set_by_gene_id
{
    my ($self, $gid, $cons_level) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_list = $self->fetch_list_by_gene_id($gid, $cons_level);
    return if !$cr_list;

    my $cr_set = OPOSSUM::ConservedRegionSet->new();
    $cr_set->param('gene_id',            $gid);
    $cr_set->param('conservation_level', $cons_level);

    foreach my $cr (@$cr_list) {
        $cr_set->add_conserved_region($cr);
    }

    return $cr_set;
}

=head2 fetch_by_upstream_downstream

 Title    : fetch_by_upstream_downstream
 Usage    : $regions = $cra->fetch_by_upstream_downstream(
				$gid, $cons_level, $upstream_bp,
				$downstream_bp);
 Function : Alias for fetch_list_by_upstream_downstream.
 Returns  : A listref of OPOSSUM::ConservedRegion objects.
 Args	  : Gene ID,
            Conservation level,
            OPTIONAL upstream BP,
            OPTIONAL downstream BP

=cut

sub fetch_by_upstream_downstream
{
    my ($self, $gid, $cons_level, $upstream_bp, $downstream_bp) = @_;

    return $self->fetch_list_by_upstream_downstream($gid, $cons_level,
        $upstream_bp, $downstream_bp);
}

=head2 fetch_list_by_upstream_downstream

 Title    : fetch_list_by_upstream_downstream
 Usage    : $regions = $cra->fetch_list_by_upstream_downstream(
				$gid, $cons_level, $upstream_bp,
				$downstream_bp);
 Function : Fetch the list of conserved regions for a given gene
            at the given level of conservation from the DB. Limit
            the conserved regions to those which fall within
            the specified amount of upstream/downstream bp of any
            promoter (TSS) of the gene. If upstream/downstream bp are
            not specified, conserved regions are limited to the boundaries
            of the promoter regions.
 Returns  : A listref of OPOSSUM::ConservedRegion objects.
 Args	  : Gene ID,
            Conservation level,
            OPTIONAL upstream BP,
            OPTIONAL downstream BP

=cut

sub fetch_list_by_upstream_downstream
{
    my ($self, $gid, $cons_level, $upstream_bp, $downstream_bp) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_list = $self->fetch_list_by_gene_id($gid, $cons_level);

    return if !$cr_list;

    if (defined $upstream_bp || defined $downstream_bp) {
        $cr_list = $self->_filter_by_promoter_search_regions(
            $gid, $cr_list, $upstream_bp, $downstream_bp
        );
    }

    return $cr_list;
}

=head2 fetch_set_by_upstream_downstream

 Title    : fetch_set_by_upstream_downstream
 Usage    : $regions = $cra->fetch_set_by_upstream_downstream(
				$gid,
                $cons_level,
                $upstream_bp,
                $downstream_bp
            );
 Function : Fetch the set of conserved regions for a given gene
            at the given level of conservation from the DB. Limit
            the conserved regions to those which fall within
            the specified amount of upstream/downstream bp of any
            promoter (TSS) of the gene. If upstream/downstream bp are
            not specified, conserved regions are limited to the boundaries
            of the promoter regions.
 Returns  : An OPOSSUM::ConservedRegionSet object.
 Args	  : Gene pair ID,
            Conservation level,
            OPTIONAL upstream BP,
            OPTIONAL downstream BP

=cut

sub fetch_set_by_upstream_downstream
{
    my ($self, $gid, $cons_level, $upstream_bp, $downstream_bp) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_list = $self->fetch_list_by_upstream_downstream(
        $gid, $cons_level, $upstream_bp, $downstream_bp
    );

    return if !$cr_list;

    my $cr_set = OPOSSUM::ConservedRegionSet->new;
    foreach my $cr (@$cr_list) {
        $cr_set->add_conserved_region($cr);
    }

    $cr_set->param('gene_id',            $gid);
    $cr_set->param('conservation_level', $cons_level);

    return $cr_set;
}

=head2 fetch_length_by_gene_id

 Title    : fetch_length_by_gene_id
 Usage    : $len = $cra->fetch_length_by_gene_id($gid, $cons_level);
 Function : Fetch the total length of the conserved regions for a
            given gene at the given level of conservation.
 Returns  : An integer length.
 Args	  : Integer gene ID and integer conservation level.

=cut

sub fetch_length_by_gene_id
{
    my ($self, $gid, $cons_level) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_set = $self->fetch_set_by_gene_id($gid, $cons_level);

    my $length = 0;

    $length = $cr_set->total_length if $cr_set;

    return $length;
}

=head2 fetch_gc_content_by_gene_id

 Title    : fetch_gc_content_by_gene_id
 Usage    : $len = $cra->fetch_gc_content_by_gene_id($gid, $cons_level);
 Function : Fetch the total gc_content of the conserved regions for a
            given gene at the given level of conservation.
 Returns  : A real gc_content.
 Args	  : Integer gene ID and integer conservation level.

=cut

sub fetch_gc_content_by_gene_id
{
    my ($self, $gid, $cons_level) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_set = $self->fetch_set_by_gene_id($gid, $cons_level);

    my $gc_content = 0;

    $gc_content = $cr_set->total_gc_content if $cr_set;

    return $gc_content;
}

# New method added by DJA on 2007/01/29

=head2 fetch_length_by_upstream_downstream

 Title    : fetch_length_by_upstream_downstream
 Usage    : $length = $cra->fetch_length_by_upstream_downstream(
				$gid, $cons_level, $upstream_bp, $downstream_bp);
 Function : Fetch the length of conserved regions for a given gene
            at the given level of conservation from the DB. Limit
            the conserved regions to those which fall within
            the specified amount of upstream/downstream bp of any
            promoter (TSS) of the gene. If upstream/downstream bp are
            not specified, conserved regions are limited to the boundaries
            of the promoter regions.
 Returns  : The total length of the conserved regions.
 Args	  : Gene ID,
            Conservation level,
            OPTIONAL upstream BP,
            OPTIONAL downstream BP

=cut

sub fetch_length_by_upstream_downstream
{
    my ($self, $gid, $cons_level, $upstream_bp, $downstream_bp) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_set = $self->fetch_set_by_upstream_downstream(
        $gid, $cons_level, $upstream_bp, $downstream_bp
    );

    return 0 if !$cr_set || $cr_set->size() == 0;

    return $cr_set->total_length();
}

=head2 fetch_gc_content_by_upstream_downstream

 Title    : fetch_gc_content_by_upstream_downstream
 Usage    : $gc_content = $cra->fetch_gc_content_by_upstream_downstream(
				$gid, $cons_level, $upstream_bp, $downstream_bp);
 Function : Fetch the gc_content of conserved regions for a given gene
            at the given level of conservation from the DB. Limit
            the conserved regions to those which fall within
            the specified amount of upstream/downstream bp of any
            promoter (TSS) of the gene. If upstream/downstream bp are
            not specified, conserved regions are limited to the boundaries
            of the promoter regions.
 Returns  : The total gc_content of the conserved regions.
 Args	  : Gene ID,
            Conservation level,
            OPTIONAL upstream BP,
            OPTIONAL downstream BP

=cut

sub fetch_gc_content_by_upstream_downstream
{
    my ($self, $gid, $cons_level, $upstream_bp, $downstream_bp) = @_;

    if (!$gid) {
        carp "must provide gene ID";
        return;
    }

    $cons_level = 1 if !$cons_level;

    my $cr_set = $self->fetch_set_by_upstream_downstream(
        $gid, $cons_level, $upstream_bp, $downstream_bp
    );

    return 0 if !$cr_set || $cr_set->size() == 0;

    return $cr_set->total_gc_content();
}

=head2 store

 Title   : store
 Usage   : $cra->store($region);
 Function: Store a single ConservedRegion object in the database.
 Args    : An OPOSSUM::ConservedRegion object
 Returns : True on success, false otherwise.

=cut

sub store
{
    my ($self, $cr) = @_;

    return if !$cr;

    if (!$cr->isa('OPOSSUM::ConservedRegion')) {
        carp "Not an OPOSSUM::ConservedRegion object";
        return;
    }

    my $sql = qq{insert into conserved_regions
        (gene_id, start, end, conservation_level, conservation, gc_content)
        values (?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert conserved regions statement - "
            . $self->errstr;
        return;
    }

    if (!$sth->execute(
            $cr->gene_id, $cr->start, $cr->end, $cr->conservation_level,
            $cr->conservation, $cr->gc_content
        )
    )
    {
        carp "Error inserting OPOSSUM::ConservedRegion - " . $sth->errstr;
        return;
    }

    return 1;
}

=head2 store_list

 Title   : store_list
 Usage   : $cra->store_list($regions, $gene_id, $cons_level);
 Function: Store a list of conserved regions in the database.
 Args    : The list of conserved regions to store (ref to a list of
           either Bio::SeqFeature::FeaturePair objects or a
           OPOSSUM::ConservedRegion objects);
           The ID of the Gene with which this conserved region
           list is associated;
           The conservation level of these conserved regions.
 Returns : True on success, false otherwise.

=cut

sub store_list
{
    my ($self, $regions, $gid, $cons_level) = @_;

    return if !$regions || !$regions->[0];

    if (   !$regions->[0]->isa('Bio::SeqFeature::FeaturePair')
        && !$regions->[0]->isa('OPOSSUM::ConservedRegion'))
    {
        carp "Not a Bio::SeqFeature::FeaturePair or OPOSSUM::ConservedRegion"
            . " list";
        return;
    }

    my $sql = qq{insert into conserved_regions
        (gene_id, start, end, conservation_level, conservation, gc_content)
        values (?,?,?,?,?,?,?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert conserved regions statement - "
            . $self->errstr;
        return;
    }

    my $ok = 1;
    foreach my $cr (@$regions) {
        if ($cr->isa('Bio::SeqFeature::FeaturePair')) {
			my ($gc_content) = $cr->get_tag_values('gc_content');
            if (!$sth->execute(
                $gid, $cr->start, $cr->end, $cons_level, $cr->score, $gc_content
            ))
            {
                carp "Error inserting Bio::SeqFeature::FeaturePair - "
                    . $sth->errstr;
                $ok = 0;
            }
        } elsif ($cr->isa('OPOSSUM::ConservedRegion')) {
            if (!$sth->execute($gid, $cr->start, $cr->end, $cons_level,
                $cr->conservation, $cr->gc_content))
            {
                carp "Error inserting OPOSSUM::ConservedRegion - "
                    . $sth->errstr;
                $ok = 0;
            }
        } else {
            carp "Not a Bio::SeqFeature::FeaturePair or "
                . " OPOSSUM::ConservedRegion object";
            return;
        }
    }

    return $ok;
}

=head2 store_set

 Title   : store_set
 Usage   : $cra->store_set($cr_set);
 Function: Store a set of conserved regions in the database.
 Args    : An OPOSSUM::ConservedRegionSet object
 Returns : True on success, false otherwise.

=cut

sub store_set
{
    my ($self, $cr_set) = @_;

    return if !$cr_set;

    if (!$cr_set->isa('OPOSSUM::ConservedRegionSet')) {
        carp "Not an OPOSSUM::ConservedRegionSet";
        return;
    }

    my $cons_level = $cr_set->param('conservation_level');
    if (!$cons_level) {
        carp
            "OPOSSUM::ConservedRegionSet conservation_level parameter not set";
        return;
    }

    my $gid = $cr_set->param('gene_id');
    if (!$gid) {
        carp "OPOSSUM::ConservedRegionSet gene_id parameter not set";
        return;
    }

    my $sql = qq{insert into conserved_regions
        (gene_id, start, end, conservation_level, conservation, gc_content)
		values (?,?,?,?,?,?,?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert conserved regions statement - "
            . $self->errstr;
        return;
    }

    my $ok = 1;
    foreach my $cr (@{$cr_set->conserved_regions}) {
        if (!$sth->execute(
            $gid, $cr->start, $cr->end, $cons_level,
			$cr->conservation, $cr->gc_content))
        {
            carp "Error inserting OPOSSUM::ConservedRegion - " . $sth->errstr;
            $ok = 0;
        }
    }

    return $ok;
}

#
# Filter conserved regions by promoter search regions. Truncate conserved
# regions at promoter search region boundaries. This routine also takes into
# account that a conserved region could overlap more than one search region
# and therefore be split into two or more new conserved regions.
#
sub _filter_by_promoter_search_regions
{
    my ($self, $gid, $cr_list, $upstream_bp, $downstream_bp) = @_;

    return if !$cr_list;

    #
    # If neither upstream nor downstream bp are provided we don't need to
    # do any filtering.
    #
    return $cr_list if !defined $upstream_bp && !defined $downstream_bp;

    my $ga = $self->db->get_GeneAdaptor;
    if (!$ga) {
        carp "error getting GeneAdaptor";
        return;
    }

    my $gene = $ga->fetch_by_id($gid);
    if (!$gene) {
        carp "error fetching Gene $gid";
        return;
    }

    my $search_regions = $gene->promoter_search_regions(
        $upstream_bp, $downstream_bp
    );

    my @filtered_crs;
    foreach my $cr (@$cr_list) {
        my $cr_start = $cr->start;
        my $cr_end   = $cr->end;

        foreach my $sr (@$search_regions) {
            my $sr_start = $sr->start;
            my $sr_end   = $sr->end;

            my $filt_cr_start = $cr_start;
            my $filt_cr_end   = $cr_end;

            # Only need to process if original CR bounds overlap current SR
            if ($cr_end >= $sr_start && $cr_start <= $sr_end) {

                # CR overlaps SR start; set new CR start to SR start
                if ($cr_start < $sr_start) {
                    $filt_cr_start = $sr_start;
                }

                # CR overlaps SR end; set new CR end to SR end
                if ($cr_end > $sr_end) {
                    $filt_cr_end = $sr_end;
                }

                #
                # This is in inner (SR) loop as a single original CR may
                # overlap multiple SRs resulting in multiple new CRs
                #
                # NOTE: we keep the conservation score of the original CR
                # even though the new CR would not really have the same
                # score if it was truncated at an SR boundary.
                #
                push @filtered_crs, OPOSSUM::ConservedRegion->new(
                    -adaptor            => $self,
                    -gene_id            => $cr->gene_id,
                    -conservation_level => $cr->conservation_level,
                    -conservation       => $cr->conservation,
					-gc_content			=> $cr->gc_content,
                    -start              => $filt_cr_start,
                    -end                => $filt_cr_end
                );
            }
        }
    }

    return @filtered_crs ? \@filtered_crs : undef;
}

1;
