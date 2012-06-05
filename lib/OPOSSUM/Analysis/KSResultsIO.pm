
=head1 NAME

OPOSSUM::Analysis::KSResultsIO - Object for the I/O of
OPOSSUM::Analysis::KSResultSet objects

=head1 AUTHOR

 Andrew Kwon, based on module by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::KSResultsIO;

use strict;

use Carp;
use OPOSSUM::Analysis::KSResult;
use OPOSSUM::Analysis::KSResultSet;

=head2 new

 Title    : new
 Usage    : $frIO = OPOSSUM::Analysis::KSResultsIO->new(
                -file	=> $in_file
            );
 Function : Create a new OPOSSUM::Analysis::KSResultsIO object.
 Returns  : An OPOSSUM::Analysis::KSResultsIO object.
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
 Returns  : An OPOSSUM::Analysis::KSResultSet object
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
 Args     : An OPOSSUM::Analysis::KSResultSet object
            Optionally a flag which, if true, excludes results in which
            a TFBS was only found for a single gene

=cut

sub write_results
{
    my ($self, $results) = @_;

    if (!$results || !$results->isa("OPOSSUM::Analysis::KSResultSet")) {
        carp "no results provided or results is not a"
            . " OPOSSUM::Analysis::KSResults object";
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

    return _write_results($fh, $results);
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
    if (!@headers || scalar @headers != 3) {
        carp "error reading KS results file column headers";
        return;
    }

    my $results = OPOSSUM::Analysis::KSResultSet->new;
    if (!$results) {
        carp "error creating KS result set";
        return;
    }

    #
    # The KS module ouputs the columns in one order (with target
    # hits/nohits before background hits/nohits), but this module outputs
    # them in the opposite order (to be compatible with the web results
    # output) so check the column headers to see which order to read them in.
    #
    if (   $headers[0] eq 'TF'
        && $headers[1] eq 'p-value'
        && $headers[2] eq 'BG_distribution')
    {
        while (my $line = <$fh>) {
            chomp $line;
            my @vals = split /\s+/, $line;
            {
                $results->add_result(
                    OPOSSUM::Analysis::KSResult->new(
                        -id         => $vals[0],
                        -p_value    => $vals[1],
                        -bg_distribution => $vals[2]
                    )
                );
            }
        }
    } else {
        carp "incorrect KS results file column headers";
        return;
    }

    return $results;
}

sub _write_results
{
    my ($fh, $results) = @_;

    return if !$results;

    print $fh "TF\tp-value\tBG_distribution\n";

    my $num_results = $results->num_results;
    my $idx         = 0;
    while ($idx < $num_results) {
        my $result = $results->get_result($idx);
        if ($result->t_hits > 0) {
            printf $fh "%s\t%0.4g\t%s\n",
                $result->id,
                $result->p_value,
                $result->bg_distribution;
        }
        $idx++;
    }

    return 1;
}

1;
