
=head1 NAME

OPOSSUM::Analysis::Cluster::FisherResultsIO - Object for the I/O of
OPOSSUM::Analysis::Cluster::FisherResultSet objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::FisherResultsIO;

use strict;

use Carp;
use OPOSSUM::Analysis::Cluster::FisherResult;
use OPOSSUM::Analysis::Cluster::FisherResultSet;

=head2 new

 Title    : new
 Usage    : $frIO = OPOSSUM::Analysis::Cluster::FisherResultsIO->new(
                -file	=> $in_file
            );
 Function : Create a new OPOSSUM::Analysis::Cluster::FisherResultsIO object.
 Returns  : An OPOSSUM::Analysis::Cluster::FisherResultsIO object.
 Args     : file - name of a file for input/output
            fh   - a filehandle for input/output

=cut

sub new
{
    my ($class, %args) = @_;

    my $file = $args{-file};
    my $fh   = $args{-fh};

    if (!$file && !$fh) {
        carp "must provide either a file name or a file handle";
        return;
    }

    if ($file && $fh) {
        carp "must provide either a file name or a file handle, not both";
        return;
    }

    if ($file) {
        open($fh, $file);
        if (!$fh) {
            carp "error opening $file - $!";
            return;
        }
    }

    my $self = bless {
        -file => $file,
        -fh   => $fh
    }, ref $class || $class;

    return $self;
}

sub DESTROY
{
    my $self = shift;

   #
   # Only close if it was opened in this module (i.e. a file name was provided
   # rather than a file handle
   #
    if ($self->file) {
        $self->close;
    }
}

=head2 fh

 Title    : fh
 Usage    : $fh = $frIO->fh() or $frIO->fh(\*FH);
 Function : Get/set the filehandle
 Returns  : A filehandle
 Args     : Optional filehandle

=cut

sub fh
{
    my ($self, $fh) = @_;

    if ($fh) {
        $self->{-fh} = $fh;
    }
    return $self->{-fh};
}

=head2 file

 Title    : file
 Usage    : $file = $frIO->file() or $frIO->file($file);
 Function : Get/set the file name
 Returns  : A file name
 Args     : Optional file name

=cut

sub file
{
    my ($self, $file) = @_;

    if ($file) {
        $self->{-file} = $file;
    }
    return $self->{-file};
}

=head2 close

 Title    : close
 Usage    : $frIO->close();
 Function : Close the filehandle if it is open.
 Returns  : Nothing
 Args     : None

=cut

sub close
{
    my $self = shift;

    if ($self->fh) {
        close($self->fh);
        $self->{-fh} = undef;
    }
}

=head2 read_results

 Title    : read_results
 Usage    : $results = $frIO->read_results();
 Function : Read results from the open filehandle.
 Returns  : An OPOSSUM::Analysis::Cluster::FisherResultSet object
 Args     : None

=cut

sub read_results
{
    my ($self) = @_;

    my $fh = $self->fh;
    if (!$fh) {
        carp "file handle is no longer valid";
        return;
    }

    if ($self->file =~ /^>{1,2}/) {
        carp "file is not open for reading";
        return;
    }

    return _read_results($fh);
}

=head2 write_results

 Title    : write_results
 Usage    : $frIO->write_results($results, $exclude_single_hits);
 Function : Write results to the open filehandle. If $exclude_single_hits
            is true, does not write results in which a TFBS was only
            found for a single gene.
 Returns  : Nothing
 Args     : An OPOSSUM::Analysis::Cluster::FisherResultSet object
            Optionally a flag which, if true, excludes results in which
            a TFBS was only found for a single gene

=cut

sub write_results
{
    my ($self, $results, $exclude_single_hits) = @_;

    if (!$results || !$results->isa("OPOSSUM::Analysis::Cluster::FisherResultSet")) {
        carp "no results provided or results is not a"
            . " OPOSSUM::Analysis::Cluster::FisherResults object";
        return;
    }

    my $fh = $self->fh;
    if (!$fh) {
        carp "file handle is no longer valid";
        return;
    }

    if ($self->file !~ /^>{1,2}/) {
        carp "file is not open for writing";
        return;
    }

    return _write_results($fh, $results, $exclude_single_hits);
}

sub _read_results
{
    my ($fh) = @_;

    my @headers;
    my $line = <$fh>;
    if ($line) {
        chomp($line);
        @headers = split /\s+/, $line;
    }
    if (!@headers || scalar @headers != 6) {
        carp "error reading Fisher results file column headers";
        return;
    }

    my $results = OPOSSUM::Analysis::Cluster::FisherResultSet->new;
    if (!$results) {
        carp "error creating Fisher result set";
        return;
    }

    #
    # The Fisher module ouputs the columns in one order (with target
    # hits/nohits before background hits/nohits), but this module outputs
    # them in the opposite order (to be compatible with the web results
    # output) so check the column headers to see which order to read them in.
    #
    if (   $headers[0] eq 'TFCluster'
        && $headers[1] eq 'tHits'
        && $headers[2] eq 'tNoHits'
        && $headers[3] eq 'bHits'
        && $headers[4] eq 'bNoHits'
        && $headers[5] eq 'p-value')
    {
        while (my $line = <$fh>) {
            chomp $line;
            my @vals = split /\s+/, $line;
            {
                $results->add_result(
                    OPOSSUM::Analysis::Cluster::FisherResult->new(
                        -id         => $vals[0],
                        -t_hits     => $vals[1],
                        -t_no_hits  => $vals[2],
                        -bg_hits    => $vals[3],
                        -bg_no_hits => $vals[4],
                        -p_value    => $vals[5]
                    )
                );
            }
        }
    } elsif ($headers[0] eq 'TFCluster'
        && $headers[3] eq 'bHits'
        && $headers[4] eq 'bNoHits'
        && $headers[1] eq 'tHits'
        && $headers[2] eq 'tNoHits'
        && $headers[5] eq 'p-value')
    {
        while (my $line = <$fh>) {
            chomp $line;
            my @vals = split /\s+/, $line;
            {
                $results->add_result(
                    OPOSSUM::Analysis::Cluster::FisherResult->new(
                        -id         => $vals[0],
                        -bg_hits    => $vals[1],
                        -bg_no_hits => $vals[2],
                        -t_hits     => $vals[3],
                        -t_no_hits  => $vals[4],
                        -p_value    => $vals[5]
                    )
                );
            }
        }
    } else {
        carp "incorrect Fisher results file column headers";
        return;
    }

    return $results;
}

sub _write_results
{
    my ($fh, $results, $exclude_single_hits) = @_;

    return if !$results;

    print $fh "TFCluster\tbHits\tbNoHits\ttHits\ttNoHits\tp-value\n";

    my $num_results = $results->num_results;
    my $idx         = 0;
    while ($idx < $num_results) {
        my $result = $results->get_result($idx);
        if ($result->t_hits > 0) {
            if (!$exclude_single_hits || $result->t_hits > 1) {
                printf $fh "%-15s %d\t%d\t%d\t%d\t%0.4g\n",
                    $result->id,
                    $result->bg_hits,
                    $result->bg_no_hits,
                    $result->t_hits,
                    $result->t_no_hits,
                    $result->p_value;
            }
        }
        $idx++;
    }

    return 1;
}

1;
