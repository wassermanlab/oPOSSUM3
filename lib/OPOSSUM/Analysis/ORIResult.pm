=head1 NAME

OPOSSUM::ORIResult.pm - module to hold the result of an ORI analysis for
a single TFBS

=head1 AUTHOR

 Shannan Ho Sui
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: shosui@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::ORIResult;

use strict;

use Carp;

=head2 new

 Title    : new
 Usage    : $or = OPOSSUM::Analysis::ORIResult->new(
					-id		=> $id,
					-t_rate		=> $t_rate,
					-bg_rate	=> $bg_rate,
					-t_gene_prop	=> $t_gene_prop,
					-bg_gene_prop	=> $bg_gene_prop,
					-ori     	=> $ori);
 Function : Create a new OPOSSUM::Analysis::ORIResult object.
 Returns  : An OPOSSUM::Analysis::ORIResult object.
 Args     : id		- TFBS ID
 	    t_hits	- number of times this TFBS was detected in the
	    		  target set of genes
 	    bg_hits	- number of times this TFBS was detected in the
	    		  background set of genes
 	    t_rate	- rate that this TFBS was detected in the target 
                          set of genes
 	    bg_rate	- rate that this TFBS was detected in the 
                          background set of genes
 	    t_gene_prop	- proportion of genes in the target set for which this
			  TFBS was detected
 	    bg_gene_prop- proportion of genes in the background set for which
	    		  this TFBS was detected
	    ori 	- measure of over-representation of the TFBS in the target
                          set versus the background rate based on TFBS density
                          and the proportion of genes containing the TFBS

=cut

sub new
{
    my ($class, %args) = @_;

    my $id = $args{-id};
    if (!$id) {
        carp "must provide ID";
	return;
    }
    my $ori = $args{-ori};
    my $t_rate = $args{-t_rate};
    my $bg_rate = $args{-bg_rate};
    my $t_gene_prop = $args{-t_gene_prop};
    my $bg_gene_prop = $args{-bg_gene_prop};

    my $self = bless {
    			-id		=> $id,
			-ori	        => $ori,
			-t_rate		=> $t_rate,
			-bg_rate	=> $bg_rate,
			-t_gene_prop	=> $t_gene_prop,
			-bg_gene_prop	=> $bg_gene_prop
		    }, ref $class || $class;

    return $self;
}

=head2 id

 Title    : id
 Usage    : $id = $or->id() or $or->id($id);
 Function : Get/set the TFBS ID of this result.
 Returns  : TFBS ID string.
 Args     : Optional TFBS ID string.

=cut

sub id
{
    my ($self, $id) = @_;

    if (defined $id) {
    	$self->{-id} = $id;
    }

    return $self->{-id};
}

=head2 ori

 Title    : ori
 Usage    : $ori = $or->ori() or $or->ori($ori);
 Function : Get/set the ori of this result.
 Returns  : real ori.
 Args     : Optional real ori.

=cut

sub ori
{
    my ($self, $ori) = @_;

    if (defined $ori) {
    	$self->{-ori} = $ori;
    }

    return $self->{-ori};
}

=head2 t_rate

 Title    : t_rate
 Usage    : $t_rate = $or->t_rate() or $or->t_rate($t_rate);
 Function : Get/set the target rate of this result.
 Returns  : real target rate.
 Args     : Optional real target rate.

=cut

sub t_rate
{
    my ($self, $t_rate) = @_;

    if (defined $t_rate) {
    	$self->{-t_rate} = $t_rate;
    }

    return $self->{-t_rate};
}

=head2 bg_rate

 Title    : bg_rate
 Usage    : $bg_rate = $or->bg_rate() or $or->bg_rate($bg_rate);
 Function : Get/set the background rate of this result.
 Returns  : real background rate.
 Args     : Optional real background rate.

=cut

sub bg_rate
{
    my ($self, $bg_rate) = @_;

    if (defined $bg_rate) {
    	$self->{-bg_rate} = $bg_rate;
    }

    return $self->{-bg_rate};
}

=head2 t_gene_prop

 Title    : t_gene_prop
 Usage    : $t_gene_prop = $or->t_gene_prop()
	    or $or->t_gene_prop($t_gene_prop);
 Function : Get/set the target gene proportion of this result.
 Returns  : integer target gene hits.
 Args     : Optional integer target gene hits.

=cut

sub t_gene_prop
{
    my ($self, $t_gene_prop) = @_;

    if (defined $t_gene_prop) {
    	$self->{-t_gene_prop} = $t_gene_prop;
    }

    return $self->{-t_gene_prop};
}

=head2 bg_gene_prop

 Title    : bg_gene_hits
 Usage    : $bg_gene_prop = $or->bg_gene_prop()
	    or $or->bg_gene_prop($bg_gene_prop);
 Function : Get/set the background gene hits of this result.
 Returns  : integer background gene hits.
 Args     : Optional integer background gene hits.

=cut

sub bg_gene_prop
{
    my ($self, $bg_gene_prop) = @_;

    if (defined $bg_gene_prop) {
    	$self->{-bg_gene_prop} = $bg_gene_prop;
    }

    return $self->{-bg_gene_prop};
}

1;
