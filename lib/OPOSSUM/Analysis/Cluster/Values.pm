=head1 NAME

OPOSSUM::Analysis::Cluster::Values - Object to store some values per TFCluster (mainly values
between predicted TFBS clusters and ChIP-Seq peaks, but could be used for others)

=head1 DESCRIPTION

This object stores values associated with each TFCluster.
This object can be passed to OPOSSUM::Analysis::Cluster::KS module.

=head1 AUTHOR

 Andrew Kwon, adapted from module by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::Values;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $values = OPOSSUM::Analysis::Cluster::Values->new();
 Function : Create a new OPOSSUM::Analysis::Cluster::Values object.
 Returns  : An OPOSSUM::Analysis::Cluster::Values object.

=cut

sub new
{
    my ($class, %args) = @_;

    my $seq_ids = $args{-seq_ids};
    my $cluster_ids   = $args{-cluster_ids};
    #my $values   = $args{-values};

    my $self = bless {
        -seq_ids           		=> undef,
        -cluster_ids            => undef,
        -seq_exists        		=> undef,
        -cluster_exists         => undef,
        _seq_cluster_values   	=> {},
        _cluster_seq_exists   	=> {},
        _params             	=> {}
    }, ref $class || $class;

    if ($seq_ids && $cluster_ids) {
        #
        # If seq ID and TFCluster ID list provided, initialize as arrays
        #
        foreach my $seq_id (@$seq_ids) {
			$self->_add_seq($seq_id);
            foreach my $cluster_id (@$cluster_ids) {
				$self->_add_cluster($cluster_id);
				my $values = []; # initialize with empty array ref
                $self->seq_cluster_values($seq_id, $cluster_id, $values);
            }
        }
    }

    return $self;
}

=head2 param

 Title    : param
 Usage    : $val = $values->param($param)
	        or $values->param($param, $value);
 Function : Get/set a value of a values parameter.
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

=head2 seq_ids

 Title    : seq_ids
 Usage    : $gids = $values->seq_ids()
 Function : Get the list of seq IDs stored in the values object.
 Returns  : A reference to a list of seq IDs.
 Args     : None.

=cut

sub seq_ids
{
    my $self = shift;

    return $self->{-seq_ids};
}

=head2 get_all_seq_ids

 Title    : get_all_seq_ids
 Usage    : $sids = $values->get_all_seq_ids()
 Function : Get the list of seq IDs. This is a synonym for the
            get variant of the seq_ids method.
 Returns  : A reference to a list of seq IDs.
 Args     : None.

=cut

sub get_all_seq_ids
{
    my $self = shift;

    return $self->seq_ids();
}

=head2 cluster_ids

 Title    : cluster_ids
 Usage    : $tfids = $values->cluster_ids()
 Function : Get the list of TF IDs stored in the values object.
 Returns  : A reference to a list of TF IDs.
 Args     : None.

=cut

sub cluster_ids
{
    my $self = shift;

    return $self->{-cluster_ids}
}

=head2 get_all_cluster_ids

 Title    : get_all_cluster_ids
 Usage    : $tfids = $values->get_all_cluster_ids()
 Function : Get the list of TF IDs. This is a synonym for the get
            variant of the cluster_ids method.
 Returns  : A reference to a list of TF IDs.
 Args     : Optionally a reference to a list of TF IDs.

=cut

sub get_all_cluster_ids
{
    my $self = shift;

    return $self->cluster_ids();
}


=head2 num_seqs

 Title    : num_seqs
 Usage    : $num = $values->num_seqs()
 Function : Get the number of seqs in the values object
 Returns  : An integer.
 Args     : None.

=cut

sub num_seqs
{
    my $self = shift;

    my $num_seqs = 0;
    if ($self->seq_ids()) {
        $num_seqs = scalar @{$self->seq_ids()};
    }

    return $num_seqs;
}

=head2 num_clusters

 Title    : num_clusters
 Usage    : $num = $values->num_clusters()
 Function : Get the number of TFs in the values object
 Returns  : An integer.
 Args     : None.

=cut

sub num_clusters
{
    my $self = shift;

    my $num_clusters = 0;
    if ($self->cluster_ids()) {
        $num_clusters = scalar @{$self->cluster_ids()};
    }

    return $num_clusters;
}

=head2 seq_exists

 Title    : seq_exists
 Usage    : $bool = $values->seq_exists($id)
 Function : Return whether the seq with the given ID exists in the
            values object.
 Returns  : Boolean.
 Args     : Gene/promoter ID.

=cut

sub seq_exists
{
    my ($self, $id) = @_;

    return $self->{-seq_exists}->{$id};
}

=head2 cluster_exists

 Title    : cluster_exists
 Usage    : $bool = $values->cluster_exists($id)
 Function : Return whether sites for the TF with the given ID exists in
            the values object.
 Returns  : Boolean.
 Args     : TF ID.

=cut

sub cluster_exists
{
    my ($self, $id) = @_;

    return $self->{-cluster_exists}->{$id};
}

=head2 exists

 Title    : exists
 Usage    : $bool = $values->exists($seq_id, $cluster_id)
 Function : Return whether the seq/TF pair with the given
            IDs exist in the values object.
 Returns  : Boolean.
 Args     : A seq ID and a TF ID.

=cut

sub exists
{
    my ($self, $seq_id, $cluster_id) = @_;

    return $self->seq_exists($seq_id) && $self->cluster_exists($cluster_id);
}


=head2 add_seq_cluster_val

 Title    : add_seq_cluster_val
 Usage    : $values->add_seq_cluster_val($seq_id, $cluster_id, $val);
 Function : Add a single val for a site
 Returns  : An integer.
 Args     : A sequence ID,
            A TFCluster ID, 
            a new val for this seq/TFCluster pair

=cut

sub add_seq_cluster_value
{
    my ($self, $seq_id, $cluster_id, $val) = @_;

    return if !defined $seq_id || !defined $cluster_id || !defined $val;

	$self->_add_seq($seq_id);
	$self->_add_cluster($cluster_id);

	push @{$self->{_seq_cluster_values}->{$seq_id}->{$cluster_id}}, $val;
	$self->{_cluster_seq_exists}->{$cluster_id}->{$seq_id} = 1;

    return scalar @{$self->{_seq_cluster_values}->{$seq_id}->{$cluster_id}};
}

=head2 seq_cluster_values

 Title    : seq_cluster_values
 Usage    : $values->seq_cluster_val($seq_id, $cluster_id, $vals);
			or $vals = $values->seq_cluster_val($seq_id, $cluster_id);
 Function : Get/set the listref of values for the given TFBS sites in the
			given sequence
 Returns  : A listref of integers.
 Args     : A sequence ID,
            A TF ID,
			Optionally a listref of integers

=cut

sub seq_cluster_values
{
    my ($self, $seq_id, $cluster_id, $vals) = @_;

    return if !defined $seq_id || !defined $cluster_id;

	if (defined($vals)) {
		$self->_add_seq($seq_id);
		$self->_add_cluster($cluster_id);
		
		$self->{_seq_cluster_values}->{$seq_id}->{$cluster_id} = $vals;
		$self->{_cluster_seq_exists}->{$cluster_id}->{$seq_id} = 1;	
	}
	
    if ($self->{_seq_cluster_values}->{$seq_id}) {
        return $self->{_seq_cluster_values}->{$seq_id}->{$cluster_id};
    }

    return;
}


=head2 cluster_seq_ids

 Title    : cluster_seq_ids
 Usage    : $ids = $values->cluster_seq_ids($cluster_id)
 Function : Get the list of sequence IDs for which sites for 
            the given TFCluster were detected.
 Returns  : A ref to a list of sequence IDs.
 Args     : A TFCluster ID. 

=cut

sub cluster_seq_ids
{
    my ($self, $cluster_id) = @_;

    return if !$cluster_id;

    my @seq_ids;
    if ($self->{_cluster_seq_exists}->{$cluster_id}) {
        @seq_ids = keys %{$self->{_cluster_seq_exists}->{$cluster_id}};
    }

    return @seq_ids ? \@seq_ids : undef;
}

=head2 all_cluster_values

 Title    : all_cluster_values
 Usage    : $val = $values->all_cluster_values($cluster_id)
 Function : For the given TFCluster, return the listref of all values 
            for all the sequences in the values object.
 Returns  : An listref of integers.
 Args     : A TFCluster ID. 

=cut

sub all_cluster_values
{
    my ($self, $cluster_id) = @_;

    return if !$cluster_id;

    my @vals;
    if ($self->{_cluster_seq_exists}->{$cluster_id}) {
        my @seq_ids = keys %{$self->{_cluster_seq_exists}->{$cluster_id}};
        foreach my $seq_id (@seq_ids) {
            push @vals, @{$self->seq_cluster_values($seq_id, $cluster_id)};
        }
    }

    return \@vals;
}

=head2 missing_seq_ids

 Title    : missing_seq_ids
 Usage    : $ids = $values->missing_seq_ids()
 Function : Get a list of missing seq/sequence IDs. For convenience,
            the values object allows storage of seqs which may have
            been entered for analysis but could not be found in the
            database.
 Returns  : A ref to a list of seq/sequence IDs.
 Args     : None.

=cut

sub missing_seq_ids
{
    $_[0]->{-missing_seq_ids};
}

=head2 missing_cluster_ids

 Title    : missing_cluster_ids
 Usage    : $ids = $values->missing_cluster_ids()
 Function : Get a list of missing TF IDs. For convenience, the values
            object allows storage of TFs which may have been entered
            for analysis but could not be found in the database.
 Returns  : A ref to a list of TF IDs.
 Args     : None.

=cut

sub missing_cluster_ids
{
    $_[0]->{-missing_cluster_ids};
}

=head2 subset

 Title    : subset
 Usage    : $subset = $values->subset(
				    -seq_ids   => $seq_ids,
				    -cluster_ids     => $cluster_ids
            );

            OR

            $subset = $values->subset(
                    -seq_start => $seq_start,
                    -seq_end   => $seq_end,
                    -cluster_start   => $cluster_start,
                    -cluster_end     => $cluster_end
            );

            OR some combination of above.

seoul Function : Get a new values object from the current one with only a
            subset of the TFs and/or seqs. The subset of TFs and 
            seqs may be either specified with explicit ID lists or by
            starting and/or ending IDs.
 Returns  : A new OPOSSUM::Analysis::Cluster::Values object.
 Args     : seq_ids    - optionally a list ref of seq IDs,
            cluster_ids      - optionally a list ref of TF IDs,
            seq_start  - optionally the starting seq ID,
            seq_end    - optionally the ending seq ID,
            cluster_start    - optionally the starting TFBS ID,
            cluster_end      - optionally the ending TFBS ID.

=cut

sub subset
{
    my ($self, %args) = @_;

    my $seq_start = $args{-seq_start};
    my $seq_end   = $args{-seq_end};
    my $seq_ids   = $args{-seq_ids};
    my $cluster_start   = $args{-cluster_start};
    my $cluster_end     = $args{-cluster_end};
    my $cluster_ids     = $args{-cluster_ids};

    my $all_seq_ids = $self->seq_ids;
    my $all_cluster_ids   = $self->cluster_ids;

    my $subset_seq_ids;
    my $subset_cluster_ids;

    my @missing_seq_ids;
    my @missing_cluster_ids;

    if (!defined $seq_ids) {
        if (!$seq_start && !$seq_end) {
            $subset_seq_ids = $all_seq_ids;
        } else {
            if (!$seq_start) {
                $seq_start = $all_seq_ids->[0];
            }
            if (!$seq_end) {
                $seq_end = $all_seq_ids->[scalar @$all_seq_ids - 1];
            }
            foreach my $seq_id (@$all_seq_ids) {
                if ($seq_id ge $seq_start && $seq_id le $seq_end) {
                    push @$subset_seq_ids, $seq_id;
                }
            }
        }
    } else {
        foreach my $seq_id (@$seq_ids) {
            if (grep(/^$seq_id$/, @$all_seq_ids)) {
                push @$subset_seq_ids, $seq_id;
            } else {
                carp "warning: seq ID $seq_id not in super set,"
                    . " omitting from subset";
                push @missing_seq_ids, $seq_id;
            }
        }
    }

    if (!defined $cluster_ids) {
        if (!$cluster_start && !$cluster_end) {
            $subset_cluster_ids = $all_cluster_ids;
        } else {
            if (!$cluster_start) {
                $cluster_start = $all_cluster_ids->[0];
            }
            if (!$cluster_end) {
                $cluster_end = $all_cluster_ids->[scalar @$all_cluster_ids - 1];
            }
            foreach my $cluster_id (@$all_cluster_ids) {
                if ($cluster_id ge $cluster_start && $cluster_id le $cluster_end) {
                    push @$subset_cluster_ids, $cluster_id;
                }
            }
        }
    } else {
        foreach my $cluster_id (@$cluster_ids) {
            if (grep(/^$cluster_id$/, @$all_cluster_ids)) {
                push @$subset_cluster_ids, $cluster_id;
            } else {
                carp "warning: TFCluster ID $cluster_id not in super set,"
                    . " omitting from subset";
                push @missing_cluster_ids, $cluster_id;
            }
        }
    }

    my $subset = OPOSSUM::Analysis::Cluster::Values->new();

    return if !$subset;

    foreach my $seq_id (@$subset_seq_ids) {
        foreach my $cluster_id (@$subset_cluster_ids) {
            $subset->seq_cluster_values(
                $seq_id,
                $cluster_id,
                $self->seq_cluster_values($seq_id, $cluster_id)
            );
        }
    }

    $subset->{-missing_seq_ids} =
        @missing_seq_ids ? \@missing_seq_ids : undef;
    $subset->{-missing_cluster_ids} = @missing_cluster_ids ? \@missing_cluster_ids : undef;

    return $subset;
}

sub _add_seq
{
    my ($self, $seq_id) = @_;
    
    return if !defined $seq_id;

    unless ($self->{-seq_exists}->{$seq_id}) {
        push @{$self->{-seq_ids}}, $seq_id;
        $self->{-seq_exists}->{$seq_id} = 1;   
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
