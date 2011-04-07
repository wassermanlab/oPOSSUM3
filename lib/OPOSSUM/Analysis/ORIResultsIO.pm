=head1 NAME

OPOSSUM::Analysis::ORIResultsIO - Object for the I/O of
OPOSSUM::Analysis::ORIResultSet objects

=head1 AUTHOR

 Shannan Ho Sui
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: shosui@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::ORIResultsIO;

use strict;

use Carp;
use OPOSSUM::Analysis::ORIResult;
use OPOSSUM::Analysis::ORIResultSet;

=head2 new

 Title    : new
 Usage    : $orIO = OPOSSUM::Analysis::ORIResultsIO->new(
					    -file	=> $in_file);
 Function : Create a new OPOSSUM::Analysis::ORIResultsIO object.
 Returns  : An OPOSSUM::Analysis::ORIResultsIO object.
 Args     : file - name of a file for input/output
 	    fh	 - a filehandle for input/output

=cut

sub new
{
    my ($class, %args) = @_;

    my $file = $args{-file};
    my $fh = $args{-fh};

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
			-file	=> $file,
			-fh	=> $fh
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
 Usage    : $fh = $orIO->fh() or $orIO->fh(\*FH);
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
 Usage    : $file = $orIO->file() or $orIO->file($file);
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
 Usage    : $orIO->close();
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
 Usage    : $results = $orIO->read_results();
 Function : Read results from the open filehandle.
 Returns  : An OPOSSUM::Analysis::ORIResultSet object
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
 Usage    : $orIO->write_results($results, $exclude_single_hits);
 Function : Write results to the open filehandle. If $exclude_single_hits
 	    is set to a true value, does not write results in which a
	    given TFBS was only found for one gene. 
 Returns  : Nothing
 Args     : An OPOSSUM::Analysis::ORIResultSet object
 	    Optionally a flag indication that results in which a given
	    TFBS was only found for one gene should not be written

=cut

sub write_results
{
    my ($self, $results, $exclude_single_hits) = @_;

    if (!$results || !$results->isa("OPOSSUM::Analysis::ORIResultSet")) {
	carp "no results provided or results is not a"
		. " OPOSSUM::Analysis::ORIResultSet object";
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
	carp "error reading ORI results file column headers";
	return;
    }
    if ($headers[0] ne 'TF'
	    || $headers[1] ne 'bRate'
	    || $headers[2] ne 'tRate'
	    || $headers[1] ne 'bProp'
	    || $headers[2] ne 'tProp'
	    || $headers[3] ne 'ORI')
    {
	carp "incorrect ORI results file column headers";
	return;
    }

    my $results = OPOSSUM::Analysis::ORIResultSet->new;
    if (!$results) {
	carp "error creating ORI result set";
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
	    $results->add_result(OPOSSUM::Analysis::ORIResult->new(
						    -id		=> $vals[0],
						    -bg_rate	=> $vals[1],
						    -t_rate	=> $vals[2],
						    -bg_prop	=> $vals[3],
						    -t_prop	=> $vals[4],
						    -ori	=> $vals[5]));
	}
    }

    return $results;
}

sub _write_results
{
    my ($fh, $results, $exclude_single_hits) = @_;

    return if !$results;

    print $fh "TF\t\tbRate\ttRate\tbProp\ttProp\tORI\n";

    my $num_results = $results->num_results;
    my $idx = 0;
    while ($idx < $num_results) {
    	my $result = $results->get_result($idx);
	if ($result->t_gene_hits > 0) {
	    if (!$exclude_single_hits || $result->t_gene_hits > 1) {
		printf $fh "%-15s %0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n",
				$result->id,
				$result->bg_rate, $result->t_rate,
				$result->bg_gene_prop, $result->t_gene_prop,
				$result->ori;
	    }
	}
	$idx++;
    }

    return 1;
}

1;
