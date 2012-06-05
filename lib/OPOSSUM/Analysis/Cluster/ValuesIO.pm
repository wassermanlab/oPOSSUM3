
=head1 NAME

OPOSSUM::Analysis::Cluster::ValuesIO - Object for the I/O of OPOSSUM::Analysis::Cluster::Values
objects

=head1 AUTHOR

 Andrew Kwon, based on module by David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::ValuesIO;

use strict;

use Carp;
use OPOSSUM::Analysis::Cluster::Values;

=head2 new

 Title    : new
 Usage    : $valsIO = OPOSSUM::Analysis::Cluster::ValuesIO->new(
                -file   => $in_file,
                -format => 'KS'
            );

 Function : Create a new OPOSSUM::Analysis::Cluster::ValuesIO object.
 Returns  : An OPOSSUM::Analysis::Cluster::ValuesIO object.
 Args     : file    - name of a file for input/output
            fh      - a filehandle for input/output
            format  - format of the file: KS
                      or 'detail'

=cut

sub new
{
    my ($class, %args) = @_;

    my $file   = $args{-file};
    my $fh     = $args{-fh};
    my $format = $args{-format};

    if (!$file && !$fh) {
        carp "must provide either a file name or a file handle";
        return;
    }

    if ($file && $fh) {
        carp "must provide either a file name or a file handle, not both";
        return;
    }

    if (!$format) {
        carp "must provide a file format";
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
        -file   => $file,
        -fh     => $fh,
        -format => $format
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
 Usage    : $fh = $valsIO->fh() or $valsIO->fh(\*FH);
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
 Usage    : $file = $valsIO->file() or $valsIO->file($file);
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

=head2 format

 Title    : format
 Usage    : $format = $valsIO->format() or $valsIO->format($format);
 Function : Get/set the file format: 'ks' for now
 Returns  : A file format
 Args     : Optional file format

=cut

sub format
{
    my ($self, $format) = @_;

    if ($format) {
        $self->{-format} = $format;
    }
    return $self->{-format};
}

=head2 close

 Title    : close
 Usage    : $valsIO->close();
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

=head2 write_vals

 Title    : write_vals
 Usage    : $valsIO->write_vals($vals);
 Function : Write vals to the open filehandle
 Returns  : Nothing
 Args     : An OPOSSUM::Analysis::Cluster::Values object

=cut

sub write_vals
{
    my ($self, $vals, $cluster_id) = @_;

    if (!$vals || !$vals->isa("OPOSSUM::Analysis::Cluster::Values")) {
        carp
            "no vals provided or vals is not an OPOSSUM::Analysis::Cluster::Values"
            . " object";
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

    my $format = $self->format;
    if (!$format) {
        carp "no file format provided";
        return;
    }

    $format = lc $format;

    # forget the formatting for now. not used.
    
    if ($cluster_id) {
        my $val_list = $vals->all_cluster_values($cluster_id);
        foreach my $val (@$val_list) {
            print $fh "$cluster_id\t$val\n";
        }
    } else {
        my $cluster_ids = $vals->cluster_ids();
        if (!$cluster_ids) {
            carp "no TFBS IDs in vals\n";
            return 0;
        }
        
        foreach my $cluster_id (@$cluster_ids) {
            my $val_list = $vals->cluster_seq_val($cluster_id);
            foreach my $val (@$val_list) {
                print $fh "$cluster_id\t$val\n";
            }
        }
    }
    
    return 1;
}



1;
