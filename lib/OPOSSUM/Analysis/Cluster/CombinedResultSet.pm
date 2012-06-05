=head1 NAME

OPOSSUM::Analysis::Cluster::CombinedResultSet.pm - module to hold the results
    of a Combined analysis

=head1 DESCRIPTION

Implements a set of CombinedResult objects

=head1 MODIFICATIONS

 Andrew Kwon, Nov. 16, 2011
 - added provisions for K-S test results

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::CombinedResultSet;

use strict;
use bignum;

use Carp;

use OPOSSUM::Analysis::Cluster::CombinedResult;

=head2 new

 Title    : new
 Usage    : $crs = OPOSSUM::Analysis::Cluster::CombinedResultSet->new();
 Function : Create a new OPOSSUM::Analysis::Cluster::CombinedResultSet object.
 Returns  : An OPOSSUM::Analysis::Cluster::CombinedResultSet object.
 Args     : None.

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        _result_hash => {},
    }, ref $class || $class;

    my $frs = $args{-fisher_result_set};
    my $zrs = $args{-zscore_result_set};
    my $ksrs = $args{-ks_result_set};

    if (   $frs && $frs->isa('OPOSSUM::Analysis::Cluster::FisherResultSet')
        && $zrs && $zrs->isa('OPOSSUM::Analysis::Cluster::ZscoreResultSet')
    )
    {
        $self->_combine_results($frs, $zrs, $ksrs);
    }

    return $self;
}

=head2 size

 Title    : size
 Usage    : $size = $crs->size();
 Function : Return the size of the set
 Returns  : An integer
 Args     : None

=cut

sub size
{
    return scalar keys %{$_[0]->{_result_hash}} || 0;
}

=head2 num_results

 Title    : num_results
 Usage    : Synonym for size

=cut

sub num_results
{
    return $_[0]->size();
}

=head2 ids

 Title    : ids
 Usage    : $ids = $crs->ids();
 Function : Return IDs of all the results in the set
 Returns  : A list or listref of result (TFCluster) IDs
 Args     : None

=cut

sub ids
{
    my @ids = keys %{$_[0]->{_result_hash}};

    if (wantarray()) {
        return @ids;
    } else {
        return @ids ? \@ids : undef;
    }
}

=head2 add_result

 Title    : add_result
 Usage    : $crs->add_result($result);
 Function : Add a new result to the set
 Returns  : The OPOSSUM::Analysis::Cluster::CombinedResult object just added
 Args     : An OPOSSUM::Analysis::Cluster::CombinedResult object

=cut

sub add_result
{
    my ($self, $result) = @_;

    return if !$result;

    unless ($result->isa("OPOSSUM::Analysis::Cluster::CombinedResult")) {
        carp "Argument is not an OPOSSUM::Analysis::Cluster::CombinedResult\n";
        return;
    }

    my $id = $result->id();
    if ($self->{_result_hash}->{$id}) {
        carp "Result with ID $id already exists in set\n";
        return;
    }

    $self->{_result_hash}->{$id} = $result;
}

=head2 get_result

 Title    : get_result
 Usage    : $result = $crs->get_result($id);
 Function : Return a result from the set by it's (TFCluster) ID
 Returns  : An OPOSSUM::Analysis::Cluster::CombinedResult object
 Args     : ID of the result (TFCluster ID)

=cut

sub get_result
{
    my ($self, $id) = @_;

    return if !defined $id;

    my $result = $self->{_result_hash}->{$id};
    unless ($result) {
        carp "No Combined result with ID $id; make sure you provided an ID"
            . " and not an index value\n";
        return;
    }

    return $result;
}

=head2 get_result_by_id

 Title    : get_result_by_id
 Usage    : Synonym for get_result()

=cut

sub get_result_by_id
{
    my ($self, $id) = @_;

    return $self->get_result($id);
}

=head2 get_list

 Title    : get_list
 Usage    : $cr_list = $crs->get_list(%params);
 Function : Get a list of results
 Returns  : A listref of OPOSSUM::Analysis::Cluster::CombinedResult objects
 Args     : Parameters specifying which results to return:
                -num_results    = The top number of results to return.
                -zscore_cutoff  = Return only results whose z-score is at
                                  least this value.
                -fisher_cutoff  = Return only results whose Fisher p-valie
                                  is at most this value.
                -ks_cutoff      = Return only results whose KS p-value
                                  is at most this value.
                -sort_by        = A field to sort on, e.g. 'id', 'zscore',
                                  'fisher_p_value'.
                -reverse        = Boolean which if true, causes list to
                                  be sorted in reverse alphanumeric order.

=cut

sub get_list
{
    my ($self, %params) = @_;

    #
    # Generally, either zscore/fisher cutoff should be specified or sorted
    # number of results.
    #
    my $nresults = $params{-num_results};
    my $zcut     = $params{-zscore_cutoff};
    my $fcut     = $params{-fisher_cutoff};
    my $kscut    = $params{-ks_cutoff};
    my $sort_by  = $params{-sort_by};
    my $reverse  = $params{-reverse};

    my @list;
    my $ids = $self->ids();
    foreach my $id (@$ids) {
        my $result = $self->get_result($id);
        my $zscore = $result->zscore();
        my $fisher_p_value = $result->fisher_p_value();
        my $ks_p_value = $result->ks_p_value();

        my $include = 1;

        if (defined $zcut && $zscore < $zcut) {
            $include = 0;
        }

        if (defined $fcut && $fisher_p_value < $fcut) {
            $include = 0;
        }
        
        if (defined $kscut && $ks_p_value < $kscut) {
            $include = 0;
        }

        if ($include) {
            push @list, $result;
        }
    }

    return if !@list;

    #
    # Sort list. Default sort on ID if nothing specified
    #
    $sort_by = 'id' if !$sort_by;

    if ($list[0]->can($sort_by)) {
        if ($sort_by eq 'id') {
            #
            # String sort
            #
            if ($reverse) {
                @list = sort {$b->id() cmp $a->id()} @list;
            } else {
                @list = sort {$a->id() cmp $b->id()} @list;
            }
        } else {
            #
            # Numeric sort. If any values are undef, sort them to the bottom,
            # i.e. treat them as +/-infinity
            #
            if ($reverse) {
                @list = sort {
                        (defined $b->$sort_by ? $b->$sort_by : -&inf)
                    <=> (defined $a->$sort_by ? $a->$sort_by : -&inf)
                } @list;
            } else {
                @list = sort {
                        (defined $a->$sort_by ? $a->$sort_by : &inf)
                    <=> (defined $b->$sort_by ? $b->$sort_by : &inf)
                } @list;
            }
        }
    } else {
        carp "Invalid sort field $sort_by\n";
    }

    if (defined $nresults && $nresults !~ /^all$/i && $nresults < scalar @list)
    {
        return [@list[0..$nresults-1]] || undef;
    } else {
        return \@list || undef;
    }
}

=head2 param

 Title    : param
 Usage    : $value = $crs->param($param);
            $crs->param($param, $value);
 Function : Get/set a parameter value
 Returns  : A parameter value
 Args     : [1] A parameter name,
            [2] On set, a parameter value

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

#
# Combine Fisher and Z-score results
#
sub _combine_results
{
    my ($self, $frs, $zrs, $ksrs) = @_;

    return if !$frs || !$zrs;

    my $frs_size = $frs->size() || 0;
    my $zrs_size = $zrs->size() || 0;

    unless ($frs_size == $zrs_size) {
        carp "Warning: Fisher and Z-score result set sizes differ\n";
    }

    my @fids = $frs->ids();
    my @zids = $zrs->ids();
    my %id_processed;
    foreach my $id (@fids, @zids) {
        next if $id_processed{$id};

        $id_processed{$id} = 1;

        my $fresult = $frs->get_result($id);
        if (!$fresult) {
            carp "No Fisher result for ID $id\n";
        }

        my $zresult = $zrs->get_result($id);
        if (!$zresult) {
            carp "No Z-score result for ID $id\n";
        }
        
        my $ksresult;
        if ($ksrs)
        {
            my $ksrs_size = $ksrs->size() || 0;
            $ksresult = $ksrs->get_result($id);
            if (!$ksresult) {
                carp "No KS-test result for ID $id\n";
            }
        }

        #if ($zresult->bg_gene_hits() == 0) {
        #    carp "No background gene hits for TF $id\n";
        #}

        my $cresult = OPOSSUM::Analysis::Cluster::CombinedResult->new(-id => $id);

        if ($fresult) {
            $cresult->t_gene_hits($fresult->t_hits());
            $cresult->bg_gene_hits($fresult->bg_hits());
            $cresult->t_gene_no_hits($fresult->t_no_hits());
            $cresult->bg_gene_no_hits($fresult->bg_no_hits());
            $cresult->fisher_p_value($fresult->p_value());
        }

        if ($zresult) {
            $cresult->t_cluster_hits($zresult->t_hits());
            $cresult->bg_cluster_hits($zresult->bg_hits());
            $cresult->t_cluster_rate($zresult->t_rate());
            $cresult->bg_cluster_rate($zresult->bg_rate());
            $cresult->zscore($zresult->z_score());
            $cresult->zscore_p_value($zresult->p_value());
        }

        if ($ksresult) {
            $cresult->ks_p_value($ksresult->p_value());
        }

        $self->add_result($cresult);
    }
}

1;
