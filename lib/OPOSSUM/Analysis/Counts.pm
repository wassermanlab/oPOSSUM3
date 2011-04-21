=head1 NAME

OPOSSUM::Analysis::Counts - Object to store the count of TFBS hits on the
gene sequences

=head1 SYNOPSIS

 my $aca = $db_adaptor->get_AnalysisCountsAdaptor();

 my $counts = $aca->fetch_counts(
     -conservation_level     => 2,
     -threshold_level        => 3,
     -search_region_level    => 1
 );

=head1 DESCRIPTION

This object stores a count of the number of times each TF profile was
found on each gene. These counts can be retrieved from the database
by the OPOSSUM::DBSQL::Analysis::CountsAdaptor. This object can be passed to
the OPOSSUM::Analysis::Fisher and OPOSSUM::Analysis::Zscore modules.

=head1 MODIFICATIONS

 2010/03/10
 - No longer takes -tf_info_set and -conserved_region_lengths_set arguments.
   These have been decoupled from this module and should be passed to the
   appropriate Fisher/Zscore analysis module separately.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

 Modified by Shannan Ho Sui on Dec 21, 2006 to accommodate schema changes

=head1 METHODS

=cut

package OPOSSUM::Analysis::Counts;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $counts = OPOSSUM::Analysis::Counts->new();
 Function : Create a new OPOSSUM::Analysis::Counts object.
 Returns  : An OPOSSUM::Analysis::Counts object.

=cut

sub new
{
    my ($class, %args) = @_;

    if ($args{-tf_info_set}) {
        carp "Support for TFInfoSet is deprecated.\n";
    }

    if ($args{-conserved_region_length_set}) {
        carp "Support for ConservedRegionLengthSet is deprecated.\n";
    }

    my $gene_ids = $args{-gene_ids};
    my $tf_ids   = $args{-tf_ids};

    my $self = bless {
        -gene_ids           => undef,
        -tf_ids             => undef,
        -gene_exists        => undef,
        -tf_exists          => undef,
        _gene_tfbs_counts   => {},
        _tfbs_gene_exists   => {},
        _params             => {}
    }, ref $class || $class;

    #
    # If gene ID and TF ID list provided, initialize counts with 0's
    #
    if ($gene_ids && $tf_ids) {
        foreach my $gene_id (@$gene_ids) {
            foreach my $tf_id (@$tf_ids) {
                $self->gene_tfbs_count($gene_id, $tf_id, 0);
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

=head2 tf_ids

 Title    : tf_ids
 Usage    : $tfids = $counts->tf_ids()
 Function : Get the list of TF IDs stored in the counts object.
 Returns  : A reference to a list of TF IDs.
 Args     : None.

=cut

sub tf_ids
{
    my $self = shift;

    return $self->{-tf_ids}
}

=head2 get_all_tf_ids

 Title    : get_all_tf_ids
 Usage    : $tfids = $counts->get_all_tf_ids()
 Function : Get the list of TF IDs. This is a synonym for the get
            variant of the tf_ids method.
 Returns  : A reference to a list of TF IDs.
 Args     : Optionally a reference to a list of TF IDs.

=cut

sub get_all_tf_ids
{
    my $self = shift;

    return $self->tf_ids();
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

=head2 tf_info_set

 Title    : tf_info_set
 Usage    : This method is deprecated.

=cut

sub tf_info_set
{
    carp "tf_info_set() is deprecated\n";
}

=head2 get_tf_info

 Title    : get_tf_info
 Usage    : This method is deprecated.

=cut

sub get_tf_info
{
    carp "get_tf_info() is deprecated\n";
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

=head2 num_tfs

 Title    : num_tfs
 Usage    : $num = $counts->num_tfs()
 Function : Get the number of TFs in the counts object
 Returns  : An integer.
 Args     : None.

=cut

sub num_tfs
{
    my $self = shift;

    my $num_tfs = 0;
    if ($self->tf_ids()) {
        $num_tfs = scalar @{$self->tf_ids()};
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

=head2 tf_exists

 Title    : tf_exists
 Usage    : $bool = $counts->tf_exists($id)
 Function : Return whether sites for the TF with the given ID exists in
            the counts object.
 Returns  : Boolean.
 Args     : TF ID.

=cut

sub tf_exists
{
    my ($self, $id) = @_;

    return $self->{-tf_exists}->{$id};
}

=head2 exists

 Title    : exists
 Usage    : $bool = $counts->exists($gene_id, $tf_id)
 Function : Return whether the gene/TF pair with the given
            IDs exist in the counts object.
 Returns  : Boolean.
 Args     : A gene ID and a TF ID.

=cut

sub exists
{
    my ($self, $gene_id, $tf_id) = @_;

    return $self->gene_exists($gene_id) && $self->tf_exists($tf_id);
}

=head2 total_cr_length

 Title    : total_cr_length
 Usage    : This method is deprecated

=cut

sub total_cr_length
{
    carp "total_cr_length() is deprecated\n";
}

=head2 tfbs_width

 Title    : tfbs_width
 Usage    : This method is deprecated

=cut

sub tfbs_width
{
    carp "tfbs_width() is deprecated\n";
}

=head2 gene_tfbs_count

 Title    : gene_tfbs_count
 Usage    : $count = $counts->gene_tfbs_count($gene_id, $tf_id);
            $counts->gene_tfbs_count($gene_id, $tf_id, $count);
 Function : Get/set the count of the number of times sites for the given
            TF were detected for the given gene/sequence.
 Returns  : An integer.
 Args     : A gene/sequence ID,
            A TF ID, 
            Optionally a new count for this gene/TF pair

=cut

sub gene_tfbs_count
{
    my ($self, $gene_id, $tf_id, $count) = @_;

    return if !defined $gene_id || !defined $tf_id;

    if (defined $count) {
        $self->{_gene_tfbs_counts}->{$gene_id}->{$tf_id} = $count;

        if ($count > 0) {
            $self->{_tfbs_gene_exists}->{$tf_id}->{$gene_id} = 1;
        }

        unless ($self->gene_exists($gene_id)) {
            $self->_add_gene($gene_id);
        }

        unless ($self->tf_exists($tf_id)) {
            $self->_add_tf($tf_id);
        }
    }

    if ($self->{_gene_tfbs_counts}->{$gene_id}) {
        return $self->{_gene_tfbs_counts}->{$gene_id}->{$tf_id} || 0;
    }

    return 0;
}

=head2 tfbs_gene_count

 Title    : tfbs_gene_count
 Usage    : $count = $counts->tfbs_gene_count($tf_id)
 Function : Get the count of the number of sequences/genes for which
            sites for the given TF were detected.
 Returns  : An integer.
 Args     : A TF ID. 

=cut

sub tfbs_gene_count
{
    my ($self, $tf_id) = @_;

    return if !$tf_id;

    my $gene_count = 0;

    if ($self->{_tfbs_gene_exists}->{$tf_id}) {
        $gene_count = scalar(keys %{$self->{_tfbs_gene_exists}->{$tf_id}}) || 0;
    }


    return $gene_count;
}

=head2 set_all_gene_tfbs_counts

 Title    : set_all_gene_tfbs_counts
 Usage    : $count = $counts->set_all_gene_tfbs_counts($data);
 Function : Set the count of the number of times sites for the given
            TF were detected for the given gene/sequence.
 Returns  : Nothing.
 Args     : An arrayref of arrayrefs of gene ID, TF ID, count, e.g.
            as returned by DBI->fetchall_arrayref.

=cut

sub set_all_gene_tfbs_counts
{
    my ($self, $data) = @_;

    return if !defined $data;

    foreach my $row (@$data) {
        my $gene_id = $row->[0];
        my $tf_id   = $row->[1];
        my $count   = $row->[2];

        $self->{_gene_tfbs_counts}->{$gene_id}->{$tf_id} = $count;

        if ($count > 0) {
            $self->{_tfbs_gene_exists}->{$tf_id}->{$gene_id} = 1;
        }

        #
        # Removed checks for efficiency
        # DJA 11/04/21
        #
        #unless ($self->gene_exists($gene_id)) {
        #    $self->_add_gene($gene_id);
        #}

        #unless ($self->tf_exists($tf_id)) {
        #    $self->_add_tf($tf_id);
        #}
    }
}

=head2 tfbs_gene_ids

 Title    : tfbs_gene_ids
 Usage    : $ids = $counts->tfbs_gene_ids($tf_id)
 Function : Get the list of gene/sequence IDs for which sites for 
            the given TF were detected.
 Returns  : A ref to a list of gene/sequence IDs.
 Args     : A TF ID. 

=cut

sub tfbs_gene_ids
{
    my ($self, $tf_id) = @_;

    return if !$tf_id;

    my @gene_ids;
    if ($self->{_tfbs_gene_exists}->{$tf_id}) {
        @gene_ids = keys %{$self->{_tfbs_gene_exists}->{$tf_id}};
    }

    return @gene_ids ? \@gene_ids : undef;
}

=head2 tfbs_count

 Title    : tfbs_count
 Usage    : $count = $counts->tfbs_count($tf_id)
 Function : For the given TF, return the total number of TFBS which appear
            for all the genes/sequences in the counts object.
 Returns  : An integer.
 Args     : A TF ID. 

=cut

sub tfbs_count
{
    my ($self, $tf_id) = @_;

    return if !$tf_id;

    my $count = 0;
    if ($self->{_tfbs_gene_exists}->{$tf_id}) {
        my @gene_ids = keys %{$self->{_tfbs_gene_exists}->{$tf_id}};
        foreach my $gene_id (@gene_ids) {
            $count += $self->gene_tfbs_count($gene_id, $tf_id);
        }
    }

    return $count;
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

=head2 missing_tf_ids

 Title    : missing_tf_ids
 Usage    : $ids = $counts->missing_tf_ids()
 Function : Get a list of missing TF IDs. For convenience, the counts
            object allows storage of TFs which may have been entered
            for analysis but could not be found in the database.
 Returns  : A ref to a list of TF IDs.
 Args     : None.

=cut

sub missing_tf_ids
{
    $_[0]->{-missing_tf_ids};
}

=head2 subset

 Title    : subset
 Usage    : $subset = $counts->subset(
				    -gene_ids   => $gene_ids,
				    -tf_ids     => $tf_ids
            );

            OR

            $subset = $counts->subset(
                    -gene_start => $gene_start,
                    -gene_end   => $gene_end,
                    -tf_start   => $tf_start,
                    -tf_end     => $tf_end
            );

            OR some combination of above.

 Function : Get a new counts object from the current one with only a
            subset of the TFs and/or genes. The subset of TFs and 
            genes may be either specified with explicit ID lists or by
            starting and/or ending IDs.
 Returns  : A new OPOSSUM::Analysis::Counts object.
 Args     : gene_ids    - optionally a list ref of gene IDs,
            tf_ids      - optionally a list ref of TF IDs,
            gene_start  - optionally the starting gene ID,
            gene_end    - optionally the ending gene ID,
            tf_start    - optionally the starting TFBS ID,
            tf_end      - optionally the ending TFBS ID.

=cut

sub subset
{
    my ($self, %args) = @_;

    my $gene_start = $args{-gene_start};
    my $gene_end   = $args{-gene_end};
    my $gene_ids   = $args{-gene_ids};
    my $tf_start   = $args{-tf_start};
    my $tf_end     = $args{-tf_end};
    my $tf_ids     = $args{-tf_ids};

    my $all_gene_ids = $self->gene_ids;
    my $all_tf_ids   = $self->tf_ids;

    my $subset_gene_ids;
    my $subset_tf_ids;

    my @missing_gene_ids;
    my @missing_tf_ids;

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

    if (!defined $tf_ids) {
        if (!$tf_start && !$tf_end) {
            $subset_tf_ids = $all_tf_ids;
        } else {
            if (!$tf_start) {
                $tf_start = $all_tf_ids->[0];
            }
            if (!$tf_end) {
                $tf_end = $all_tf_ids->[scalar @$all_tf_ids - 1];
            }
            foreach my $tf_id (@$all_tf_ids) {
                if ($tf_id ge $tf_start && $tf_id le $tf_end) {
                    push @$subset_tf_ids, $tf_id;
                }
            }
        }
    } else {
        foreach my $tf_id (@$tf_ids) {
            if (grep(/^$tf_id$/, @$all_tf_ids)) {
                push @$subset_tf_ids, $tf_id;
            } else {
                carp "warning: TF ID $tf_id not in super set,"
                    . " omitting from subset";
                push @missing_tf_ids, $tf_id;
            }
        }
    }

    my $subset = OPOSSUM::Analysis::Counts->new();

    return if !$subset;

    foreach my $gene_id (@$subset_gene_ids) {
        foreach my $tf_id (@$subset_tf_ids) {
            $subset->gene_tfbs_count(
                $gene_id,
                $tf_id,
                $self->gene_tfbs_count($gene_id, $tf_id)
            );
        }
    }

    $subset->{-missing_gene_ids} =
        @missing_gene_ids ? \@missing_gene_ids : undef;
    $subset->{-missing_tf_ids} = @missing_tf_ids ? \@missing_tf_ids : undef;

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

sub _add_tf
{
    my ($self, $tf_id) = @_;
    
    return if !defined $tf_id;

    unless ($self->{-tf_exists}->{$tf_id}) {
        push @{$self->{-tf_ids}}, $tf_id;
        $self->{-tf_exists}->{$tf_id} = 1;   
    }
}

1;
