=head1 NAME

OPOSSUM::Analysis::Values - Object to store some values per TF (mainly values
between predicted TFBSs and ChIP-Seq peaks, but could be used for others)

=head1 DESCRIPTION

This object stores values associated with each TF.
This object can be passed to OPOSSUM::Analysis::KS module.

=head1 AUTHOR

 Andrew Kwon, adapted from module by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Values;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $values = OPOSSUM::Analysis::Values->new();
 Function : Create a new OPOSSUM::Analysis::Values object.
 Returns  : An OPOSSUM::Analysis::Values object.

=cut

sub new
{
    my ($class, %args) = @_;

    my $seq_ids = $args{-seq_ids};
    my $tf_ids   = $args{-tf_ids};
    #my $values   = $args{-values};

    my $self = bless {
        -seq_ids           => undef,
        -tf_ids             => undef,
        -seq_exists        => undef,
        -tf_exists          => undef,
        _seq_tfbs_values   => {},
        _tfbs_seq_exists   => {},
        _params             => {}
    }, ref $class || $class;

    if ($seq_ids && $tf_ids) {
        #
        # If seq ID and TF ID list provided, initialize as arrays
        #
        foreach my $seq_id (@$seq_ids) {
			$self->_add_seq($seq_id);
            foreach my $tf_id (@$tf_ids) {
				$self->_add_tf($tf_id);
				my $values = []; # initialize with empty array ref
                $self->seq_tfbs_values($seq_id, $tf_id, $values);
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

=head2 tf_ids

 Title    : tf_ids
 Usage    : $tfids = $values->tf_ids()
 Function : Get the list of TF IDs stored in the values object.
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
 Usage    : $tfids = $values->get_all_tf_ids()
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

=head2 num_tfs

 Title    : num_tfs
 Usage    : $num = $values->num_tfs()
 Function : Get the number of TFs in the values object
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

=head2 tf_exists

 Title    : tf_exists
 Usage    : $bool = $values->tf_exists($id)
 Function : Return whether sites for the TF with the given ID exists in
            the values object.
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
 Usage    : $bool = $values->exists($seq_id, $tf_id)
 Function : Return whether the seq/TF pair with the given
            IDs exist in the values object.
 Returns  : Boolean.
 Args     : A seq ID and a TF ID.

=cut

sub exists
{
    my ($self, $seq_id, $tf_id) = @_;

    return $self->seq_exists($seq_id) && $self->tf_exists($tf_id);
}


=head2 add_seq_tfbs_val

 Title    : add_seq_tfbs_val
 Usage    : $values->add_seq_tfbs_val($seq_id, $tf_id, $val);
 Function : Add a single val for a site
 Returns  : An integer.
 Args     : A sequence ID,
            A TF ID, 
            a new val for this seq/TF pair

=cut

sub add_seq_tfbs_value
{
    my ($self, $seq_id, $tf_id, $val) = @_;

    return if !defined $seq_id || !defined $tf_id || !defined $val;

	$self->_add_seq($seq_id);
	$self->_add_tf($tf_id);

	push @{$self->{_seq_tfbs_values}->{$seq_id}->{$tf_id}}, $val;
	$self->{_tfbs_seq_exists}->{$tf_id}->{$seq_id} = 1;

    return scalar @{$self->{_seq_tfbs_values}->{$seq_id}->{$tf_id}};
}

=head2 seq_tfbs_values

 Title    : seq_tfbs_values
 Usage    : $values->seq_tfbs_val($seq_id, $tf_id, $vals);
			or $vals = $values->seq_tfbs_val($seq_id, $tf_id);
 Function : Get/set the listref of values for the given TFBS sites in the
			given sequence
 Returns  : A listref of integers.
 Args     : A sequence ID,
            A TF ID,
			Optionally a listref of integers

=cut

sub seq_tfbs_values
{
    my ($self, $seq_id, $tf_id, $vals) = @_;

    return if !defined $seq_id || !defined $tf_id;

	if (defined($vals)) {
		$self->_add_seq($seq_id);
		$self->_add_tf($tf_id);
		
		$self->{_seq_tfbs_values}->{$seq_id}->{$tf_id} = $vals;
		$self->{_tfbs_seq_exists}->{$tf_id}->{$seq_id} = 1;	
	}
	
    if ($self->{_seq_tfbs_values}->{$seq_id}) {
        return $self->{_seq_tfbs_values}->{$seq_id}->{$tf_id};
    }

    return;
}


=head2 tfbs_seq_ids

 Title    : tfbs_seq_ids
 Usage    : $ids = $values->tfbs_seq_ids($tf_id)
 Function : Get the list of sequence IDs for which sites for 
            the given TF were detected.
 Returns  : A ref to a list of sequence IDs.
 Args     : A TF ID. 

=cut

sub tfbs_seq_ids
{
    my ($self, $tf_id) = @_;

    return if !$tf_id;

    my @seq_ids;
    if ($self->{_tfbs_seq_exists}->{$tf_id}) {
        @seq_ids = keys %{$self->{_tfbs_seq_exists}->{$tf_id}};
    }

    return @seq_ids ? \@seq_ids : undef;
}

=head2 all_tfbs_values

 Title    : all_tfbs_values
 Usage    : $val = $values->all_tfbs_values($tf_id)
 Function : For the given TF, return the listref of all values 
            for all the sequences in the values object.
 Returns  : An listref of integers.
 Args     : A TF ID. 

=cut

sub all_tfbs_values
{
    my ($self, $tf_id) = @_;

    return if !$tf_id;

    my @vals;
    if ($self->{_tfbs_seq_exists}->{$tf_id}) {
        my @seq_ids = keys %{$self->{_tfbs_seq_exists}->{$tf_id}};
        foreach my $seq_id (@seq_ids) {
            push @vals, @{$self->seq_tfbs_values($seq_id, $tf_id)};
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

=head2 missing_tf_ids

 Title    : missing_tf_ids
 Usage    : $ids = $values->missing_tf_ids()
 Function : Get a list of missing TF IDs. For convenience, the values
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
 Usage    : $subset = $values->subset(
				    -seq_ids   => $seq_ids,
				    -tf_ids     => $tf_ids
            );

            OR

            $subset = $values->subset(
                    -seq_start => $seq_start,
                    -seq_end   => $seq_end,
                    -tf_start   => $tf_start,
                    -tf_end     => $tf_end
            );

            OR some combination of above.

seoul Function : Get a new values object from the current one with only a
            subset of the TFs and/or seqs. The subset of TFs and 
            seqs may be either specified with explicit ID lists or by
            starting and/or ending IDs.
 Returns  : A new OPOSSUM::Analysis::Values object.
 Args     : seq_ids    - optionally a list ref of seq IDs,
            tf_ids      - optionally a list ref of TF IDs,
            seq_start  - optionally the starting seq ID,
            seq_end    - optionally the ending seq ID,
            tf_start    - optionally the starting TFBS ID,
            tf_end      - optionally the ending TFBS ID.

=cut

sub subset
{
    my ($self, %args) = @_;

    my $seq_start = $args{-seq_start};
    my $seq_end   = $args{-seq_end};
    my $seq_ids   = $args{-seq_ids};
    my $tf_start   = $args{-tf_start};
    my $tf_end     = $args{-tf_end};
    my $tf_ids     = $args{-tf_ids};

    my $all_seq_ids = $self->seq_ids;
    my $all_tf_ids   = $self->tf_ids;

    my $subset_seq_ids;
    my $subset_tf_ids;

    my @missing_seq_ids;
    my @missing_tf_ids;

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

    my $subset = OPOSSUM::Analysis::Values->new();

    return if !$subset;

    foreach my $seq_id (@$subset_seq_ids) {
        foreach my $tf_id (@$subset_tf_ids) {
            $subset->seq_tfbs_values(
                $seq_id,
                $tf_id,
                $self->seq_tfbs_values($seq_id, $tf_id)
            );
        }
    }

    $subset->{-missing_seq_ids} =
        @missing_seq_ids ? \@missing_seq_ids : undef;
    $subset->{-missing_tf_ids} = @missing_tf_ids ? \@missing_tf_ids : undef;

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
