
=head1 NAME

OPOSSUM::Analysis::ZscoreResultsIO - Object for the I/O of
OPOSSUM::Analysis::ZscoreResultSet objects

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::ZscoreResultsIO;

use strict;

use Carp;
use OPOSSUM::Analysis::ZscoreResult;
use OPOSSUM::Analysis::ZscoreResultSet;

=head2 new

 Title    : new
 Usage    : $zrIO = OPOSSUM::Analysis::ZscoreResultsIO->new(
                -file   => $in_file
            );
 Function : Create a new OPOSSUM::Analysis::ZscoreResultsIO object.
 Returns  : An OPOSSUM::Analysis::ZscoreResultsIO object.
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
        },
        ref $class || $class;

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
 Usage    : $fh = $zrIO->fh() or $zrIO->fh(\*FH);
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
 Usage    : $file = $zrIO->file() or $zrIO->file($file);
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
 Usage    : $zrIO->close();
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
 Usage    : $results = $zrIO->read_results();
 Function : Read results from the open filehandle.
 Returns  : An OPOSSUM::Analysis::ZscoreResultSet object
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

    my $format = $self->format;
    if (!$format) {
        carp "no file format provided";
        return;
    }

    return _read_results($fh);
}

=head2 write_results

 Title    : write_results
 Usage    : $zrIO->write_results($results, $exclude_single_hits);
 Function : Write results to the open filehandle. If $exclude_single_hits
            is set to a true value, does not write results in which a
            given TFBS was only found for one gene. 
 Returns  : Nothing
 Args     : An OPOSSUM::Analysis::ZscoreResultSet object
            Optionally a flag indication that results in which a given
            TFBS was only found for one gene should not be written

=cut

sub write_results
{
    my ($self, $results, $exclude_single_hits) = @_;

    if (!$results || !$results->isa("OPOSSUM::Analysis::ZscoreResultSet")) {
        carp "no results provided or results is not a"
            . " OPOSSUM::Analysis::ZscoreResultSet object";
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
        carp "error reading Zscore results file column headers";
        return;
    }
    if (   $headers[0] ne 'TF'
        || $headers[1] ne 'bHits'
        || $headers[2] ne 'tHits'
        || $headers[1] ne 'bRate'
        || $headers[2] ne 'tRate'
        || $headers[3] ne 'Z-score')
    {
        carp "incorrect Zscore results file column headers";
        return;
    }

    my $results = OPOSSUM::Analysis::ZscoreResultSet->new;
    if (!$results) {
        carp "error creating Zscore result set";
        return;
    }
    while ($line = <$fh>) {
        chomp($line);
#
# Won't work with scientific notation
#if ($line
#    =~ /^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d|\d*\.\d+)\s+(\d|\d*\.\d+)\s+(\d|\d*\.\d+)\s+(\d|\d*\.\d+)\s*$/)
        my @vals = split /\s+/, $line;
        {
            $results->add_result(
                OPOSSUM::Analysis::ZscoreResult->new(
                    -id      => $vals[0],
                    -bg_hits => $vals[1],
                    -t_hits  => $vals[2],
                    -bg_rate => $vals[3],
                    -t_rate  => $vals[4],
                    -z_score => $vals[5]
                )
            );
            #-p_value    => $vals[6]));
        }
    }

    return $results;
}

sub _write_results
{
    my ($fh, $results, $exclude_single_hits) = @_;

    return if !$results;

    print $fh "TF\t\tbHits\ttHits\tbRate\ttRate\tZ-score\n";

    my $num_results = $results->num_results;
    my $idx         = 0;
    while ($idx < $num_results) {
        my $result = $results->get_result($idx);
        if ($result->t_gene_hits > 0) {
            if (!$exclude_single_hits || $result->t_gene_hits > 1) {
                printf $fh "%-15s %d\t%d\t%0.5f\t%0.5f\t%0.4f\n",
                    $result->id,
                    $result->bg_hits, $result->t_hits,
                    $result->bg_rate, $result->t_rate,
                    $result->z_score;
            }
        }
        $idx++;
    }

    return 1;
}

1;
