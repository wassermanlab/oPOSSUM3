=head1 NAME

OPOSSUM::ConservedTFBSDimer - ConservedTFBSDimer object (conserved_tfbss
DB record)

=head1 DESCRIPTION

A ConservedTFBSDimer object models a record retrieved from the
conserved_tfbss table of the oPOSSUM DB with extra fields for the two species
repeat types.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::ConservedTFBSDimer;

use strict;

use Carp;


=head2 new

 Title   : new
 Usage   : $tfbs = OPOSSUM::ConservedTFBSDimer->new(
			    -tf_id		    => 1,
			    -gene_pair_id	=> 1,
			    -start1		    => 646823928,
			    -end1		    => 646823937,
			    -rel_start1		=> 21,
			    -rel_end1		=> 30,
			    -strand1		=> 1,
			    -score1		    => 1.897,
			    -rel_score1		=> 0.765,
			    -seq1		    => 'CCAAGGATAG',
			    -start2		    => 78459304,
			    -end2		    => 78459313,
			    -rel_start1		=> 237,
			    -rel_end1		=> 246,
			    -strand2		=> -1,
			    -score2		    => 2.345,
			    -rel_score2		=> 0.854,
			    -seq2		    => 'CCAAGGAGAG',
			    -conservation_level	=> 1,
			    -conservation	=> 0.832);

 Function: Construct a new ConservedTFBSDimer object
 Returns : a new OPOSSUM::ConservedTFBSDimer object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
		    %args
		}, ref $class || $class;

    return $self;
}

=head2 param

 Title    : param
 Usage    : $value = $ctfs->param($param)
    or $ctfs->param($param, $value);
 Function : Get/set the value of a parameter
 Returns  : Value of the named parameter
 Args     : [1] name of a parameter
            [2] on set, the value of the parameter

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

=head2 id

 Title   : id
 Usage   : $id = $ctfs->id() or $ctfs->id($id);

 Function: Get/set the ID of the TFBS profile (matrix) associated with
           this conserved TFBS.
 Returns : The TFBS profile ID.
 Args    : None or a new ID.

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
        $self->{-id} = $id;
    }
    return $self->{-id};
}

=head2 tf_id

 Title   : tf_id
 Usage   : $id = $ctfs->tf_id() or $ctfs->tf_id($id);

 Function: Synonymous to the 'id' method.

=cut

sub tf_id
{
    my ($self, $id) = @_;

    return $self->id($id);
}

=head2 gene_pair_id

 Title   : gene_pair_id
 Usage   : $gene_pair_id = $ctfs->gene_pair_id()
           or $ctfs->gene_pair_id($gene_pair_id);

 Function: Get/set the ID of the GenePair object associated with this
 	   conserved TFBS.
 Returns : A string.
 Args    : None or a new GenePair ID.

=cut

sub gene_pair_id
{
    my ($self, $gene_pair_id) = @_;

    if (defined $gene_pair_id) {
	$self->{-gene_pair_id} = $gene_pair_id;
    }
    return $self->{-gene_pair_id};
}

=head2 start1

 Title   : start1
 Usage   : $start = $ctfs->start1() or $ctfs->start1($start);

 Function: Get/set the start position of this conserved TF site on the
	   species 1 promoter pair sequence.
 Returns : An integer.
 Args    : None or a new start position.

=cut

sub start1
{
    my ($self, $start) = @_;

    if ($start) {
	$self->{-start1} = $start;
    }
    return $self->{-start1};
}

=head2 start2

 Title   : start2
 Usage   : $start = $ctfs->start2() or $ctfs->start2($start);

 Function: Get/set the start position of this conserved TF site on the
	   species 2 promoter pair sequence.
 Returns : An integer.
 Args    : None or a new start position.

=cut

sub start2
{
    my ($self, $start) = @_;

    if ($start) {
	$self->{-start2} = $start;
    }
    return $self->{-start2};
}

=head2 end1

 Title   : end1
 Usage   : $end = $ctfs->end1() or $ctfs->end1($end);

 Function: Get/set the end position of this conserved TF site on the
	   species 1 promoter pair sequence.
 Returns : An integer.
 Args    : None or a new end position.

=cut

sub end1
{
    my ($self, $end) = @_;

    if ($end) {
	$self->{-end1} = $end;
    }
    return $self->{-end1};
}

=head2 end2

 Title   : end2
 Usage   : $end = $ctfs->end2() or $ctfs->end2($end);

 Function: Get/set the end position of this conserved TF site on the
	   species 2 promoter pair sequence.
 Returns : An integer.
 Args    : None or a new end position.

=cut

sub end2
{
    my ($self, $end) = @_;

    if ($end) {
	$self->{-end2} = $end;
    }
    return $self->{-end2};
}

=head2 rel_start1

 Title   : rel_start1
 Usage   : $rel_start = $ctfs->rel_start1()
	   or $ctfs->rel_start1($rel_start);

 Function: Get/set the relative start position of this conserved TFBS
           on the species 1 sequence.
 Returns : An integer.
 Args    : None or a new relative start position.

=cut

sub rel_start1
{
    my ($self, $rel_start) = @_;

    if ($rel_start) {
	$self->{-rel_start1} = $rel_start;
    }
    return $self->{-rel_start1};
}

=head2 rel_start2

 Title   : rel_start2
 Usage   : $rel_start = $ctfs->rel_start2()
	   or $ctfs->rel_start2($rel_start);

 Function: Get/set the relative start position of this conserved TFBS
	   on the species 2 sequence.
 Returns : An integer.
 Args    : None or a new ative start position.

=cut

sub rel_start2
{
    my ($self, $rel_start) = @_;

    if ($rel_start) {
	$self->{-rel_start2} = $rel_start;
    }
    return $self->{-rel_start2};
}

=head2 rel_end1

 Title   : rel_end1
 Usage   : $rel_end = $ctfs->rel_end1() or $ctfs->rel_end1($rel_end);

 Function: Get/set the relative end position of this conserved TFBS
	   on the species 1 sequence.
 Returns : An integer.
 Args    : None or a new relative end position.

=cut

sub rel_end1
{
    my ($self, $rel_end) = @_;

    if ($rel_end) {
	$self->{-rel_end1} = $rel_end;
    }
    return $self->{-rel_end1};
}

=head2 rel_end2

 Title   : rel_end2
 Usage   : $rel_end = $ctfs->rel_end2() or $ctfs->rel_end2($rel_end);

 Function: Get/set the relative end position of this conserved TFBS
	   on the species 2 sequence.
 Returns : An integer.
 Args    : None or a new relative end position.

=cut

sub rel_end2
{
    my ($self, $rel_end) = @_;

    if ($rel_end) {
	$self->{-rel_end2} = $rel_end;
    }
    return $self->{-rel_end2};
}

=head2 strand1

 Title   : strand1
 Usage   : $strand = $ctfs->strand1() or $ctfs->strand1($strand);

 Function: Get/set the strand of this conserved TFBS on the
	   species 1 promoter pair sequence.
 Returns : 1 or -1.
 Args    : None or a new strand.

=cut

sub strand1
{
    my ($self, $strand) = @_;

    if ($strand) {
	$self->{-strand1} = $strand;
    }
    return $self->{-strand1};
}

=head2 strand2

 Title   : strand2
 Usage   : $strand = $ctfs->strand2() or $ctfs->strand2($strand);

 Function: Get/set the strand of this conserved TFBS on the
	   species 2 promoter pair sequence.
 Returns : 1 or -1.
 Args    : None or a new strand.

=cut

sub strand2
{
    my ($self, $strand) = @_;

    if ($strand) {
	$self->{-strand2} = $strand;
    }
    return $self->{-strand2};
}

=head2 seq1

 Title   : seq1
 Usage   : $seq = $ctfs->seq1() or $ctfs->seq1($seq);

 Function: Get/set the sequence of this conserved TF site on the
	   species 1 promoter pair sequence.
 Returns : A string.
 Args    : None or a new sequence.

=cut

sub seq1
{
    my ($self, $seq) = @_;

    if ($seq) {
	$self->{-seq1} = $seq;
    }
    return $self->{-seq1};
}

=head2 seq2

 Title   : seq2
 Usage   : $seq = $ctfs->seq2() or $ctfs->seq2($seq);

 Function: Get/set the sequence of this conserved TF site on the
	   species 2 promoter pair sequence.
 Returns : A string.
 Args    : None or a new sequence.

=cut

sub seq2
{
    my ($self, $seq) = @_;

    if ($seq) {
	$self->{-seq2} = $seq;
    }
    return $self->{-seq2};
}

=head2 score1

 Title   : score1
 Usage   : $score = $ctfs->score1() or $ctfs->score1($score);

 Function: Get/set the matrix score of this conserved TFBS on the
	   species 1 sequence.
 Returns : A real number.
 Args    : None or a new score.

=cut

sub score1
{
    my ($self, $score) = @_;

    if ($score) {
        $self->{-score1} = $score;
    }
    return $self->{-score1};
}

=head2 score2

 Title   : score2
 Usage   : $score = $ctfs->score2() or $ctfs->score2($score);

 Function: Get/set the matrix score of this conserved TFBS on the
	   species 2 sequence.
 Returns : A real number.
 Args    : None or a new score.

=cut

sub score2
{
    my ($self, $score) = @_;

    if ($score) {
        $self->{-score2} = $score;
    }
    return $self->{-score2};
}

=head2 rel_score1

 Title   : rel_score1
 Usage   : $score = $ctfs->rel_score1() or $ctfs->rel_score1($score);

 Function: Get/set the matrix relative score of this conserved TFBS on the
	   species 1 sequence.
 Returns : A real number.
 Args    : None or a new relative score.

=cut

sub rel_score1
{
    my ($self, $score) = @_;

    if ($score) {
        $self->{-rel_score1} = $score;
    }
    return $self->{-rel_score1};
}

=head2 rel_score2

 Title   : rel_score2
 Usage   : $score = $ctfs->rel_score2() or $ctfs->rel_score2($score);

 Function: Get/set the matrix relative score of this conserved TFBS on the
	   species 2 sequence.
 Returns : A real number.
 Args    : None or a new relative score.

=cut

sub rel_score2
{
    my ($self, $score) = @_;

    if ($score) {
        $self->{-rel_score2} = $score;
    }
    return $self->{-rel_score2};
}

=head2 conservation_level

 Title   : conservation_level
 Usage   : $level = $ctfs->conservation_level()
	   or $ctfs->conservation_level($level);

 Function: Get/set the conservation level of the region that this conserved
	   TF site falls within.
 Returns : An integer.
 Args    : None or a new conservation level.

=cut

sub conservation_level
{
    my ($self, $conservation_level) = @_;

    if (defined $conservation_level) {
	$self->{-conservation_level} = $conservation_level;
    }
    return $self->{-conservation_level};
}

=head2 conservation

 Title   : conservation
 Usage   : $cons = $ctfs->conservation() or $ctfs->conservation($cons);

 Function: Get/set the conservation score of the region that this conserved
	   TF site falls within.
 Returns : A real number.
 Args    : None or a new conservation score.

=cut

sub conservation
{
    my ($self, $conservation) = @_;

    if (defined $conservation) {
	$self->{-conservation} = $conservation;
    }
    return $self->{-conservation};
}

=head2 repeat_type1

 Title   : repeat_type1
 Usage   : $rt = $ctfs->repeat_type1() or $ctfs->repeat_type1($rt);

 Function: Get/set the repeat type of the TFBS dimer.
 Returns : A real number.
 Args    : None or a new repeat type.

=cut

sub repeat_type1
{
    my ($self, $rt) = @_;

    if (defined $rt) {
        $self->{-repeat_type1} = $rt;
    }
    return $self->{-repeat_type1};
}

=head2 repeat_type2

 Title   : repeat_type2
 Usage   : $rt = $ctfs->repeat_type2() or $ctfs->repeat_type2($rt);

 Function: Get/set the repeat type of the TFBS dimer.
 Returns : A real number.
 Args    : None or a new repeat type.

=cut

sub repeat_type2
{
    my ($self, $rt) = @_;

    if (defined $rt) {
        $self->{-repeat_type2} = $rt;
    }
    return $self->{-repeat_type2};
}

=head2 tf_info

 Title   : tf_info
 Usage   : $tf_info = $ctfs->tf_info() or $ctfs->tf_info($tf_info);

 Function: Get/set the TFInfo object associated with this conserved TFBS
 Returns : An OPOSSUM::TFInfo object.
 Args    : None or a new OPOSSUM::TFInfo object.

=cut

sub tf_info
{
    my ($self, $info) = @_;

    if (defined $info) {
	if ($info->isa("OPOSSUM::TFInfo")) {
	    $self->{-info} = $info;
	} else {
	    carp "not an OPOSSUM::TFInfo";
	    return undef;
	}
    }
    return $self->{-info};
}

=head2 gene_pair

 Title   : gene_pair
 Usage   : $gp = $ctfs->gene_pair() or $ctfs->gene_pair($gp);

 Function: Get/set the GenePair object associated with this conserved
	   TFBS.
 Returns : An OPOSSUM::GenePair object.
 Args    : None or a new OPOSSUM::GenePair object.

=cut

sub gene_pair
{
    my ($self, $gp) = @_;

    if ($gp) {
	if ($gp->isa("OPOSSUM::GenePair")) {
	    $self->{-gene_pair} = $gp;
	} else {
	    carp "not an OPOSSUM::GenePair";
	    return undef;
	}
    }
    return $self->{-gene_pair};
}


1;
