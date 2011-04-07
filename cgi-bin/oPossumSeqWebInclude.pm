#
# This module should be included in all the oPossumSeq*Web.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to the sequence-based oPossum variants.
#

use oPossumWebOpt;
use oPossumSeqWebOpt;

use File::Temp qw/ tempdir /;

use TFBS::MatrixSet;
use TFBS::Matrix::PFM;
use TFBS::Matrix::PWM;

use strict;

sub get_seq_file
{
    my ($self, $seq_input_method, $tempdir) = @_;

    my $q = $self->query;

    my $filename;
    if ($seq_input_method eq "paste") {
        my $sl = $q->param('seq_list');
        if (!$sl) {
            $self->_error("Sequence input not specified");
            return;
        }

        $filename = "$tempdir/seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target sequences file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($seq_input_method eq "upload") {
        my $file = $q->param('seq_file');
        my $fh   = $q->upload('seq_file');

        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("File $file is empty\n");
            return;
        }

        $filename = "$tempdir/seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target sequences file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown sequence input method");
        return;
    }

    return $filename;
}

sub get_back_seq_file
{
    my ($self, $seq_input_method, $tempdir) = @_;

    my $q = $self->query;

    my $filename;
    if ($seq_input_method eq "default") {
        my $bg_seq_set_files = BG_SEQ_SET_FILES;
        my $bg_seq_set_key = $q->param('bg_seq_set_key');
        $filename = ABS_HTDOCS_DATA_PATH . "/"
            . $bg_seq_set_files->{$bg_seq_set_key};
        #printf STDERR "background seq file: $filename\n";
    } elsif ($seq_input_method eq "paste") {
        my $sl = $q->param('bg_seq_list');

        if (!$sl) {
            $self->_error("Background sequence input not specified");
            return;
        }

        $filename = "$tempdir/back_seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error(
                "Unable to create background sequences file $filename\n"
            );
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($seq_input_method eq "upload") {
        my $file = $q->param('bg_seq_file');
        my $fh   = $q->upload('bg_seq_file');
        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("File $file is empty\n");
            return;
        }

        $filename = "$tempdir/back_seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error(
                "Unable to create background sequences file $filename\n"
            );
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown background sequence input method");
        return;
    }

    return $filename;
}

sub write_matrix_file
{
    my ($self, $matrix_spec, $tempdir) = @_;

    my $filename = "$tempdir/matrices.txt";

    unless (open(FH, ">$filename")) {
        $self->_error("Unable to create matrix file $filename\n");
        return;
    }
    
    if ($matrix_spec->isa('TFBS::MatrixSet')) {
        my $iter = $matrix_spec->Iterator();
        while (my $matrix = $iter->next()) {
            printf FH ">%s\n", $matrix->name();
            printf FH "%s\n", $matrix->prettyprint();
        }
    } elsif ($matrix_spec->isa('TFBS::Matrix')) {
        printf FH ">%s\n", $matrix_spec->name();
        printf FH "%s\n", $matrix_spec->prettyprint();
    } else {
        $self->_error(
            "Unknown matrix specification - neither a Matrix nor a MatrixSet");
        return;
    }

    close(FH);

    return $filename;
}

sub upload_matrix_file
{
    my ($self, $fh) = @_;

    return if !$fh;

    my $text;
    while (my $line = <$fh>) {
        $text .= $line;
    }
    close($fh);

    return $text;
}

sub parse_matrix_text
{
    my ($self, $text) = @_;

    return if !$text;

    my $matrix_set    = TFBS::MatrixSet->new();
    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    my $matrix_count  = 0;
    my @lines         = split /\n/, $text;
    foreach my $line (@lines) {
        chomp $line;
        next if !$line;

        if ($line =~ /^>\s*(\S+)/) {
            $name = $1;
        } else {
            if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            } elsif ($line =~ /^\s*\d+/) {
                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            } else {
                next;
            }
            $line_count++;

            if ($line_count == 4) {
                my $id = sprintf "matrix%d", $matrix_count + 1;

                unless ($name) {
                    $name = $id;
                }

                #
                # Simplistic determination of whether matrix looks more like
                # a PWM than a PFM.
                #
                my $matrix_type = 'PFM';
                if ($matrix_string =~ /\d*\.\d+/) {
                    $matrix_type = 'PWM';
                }

                my $matrix;
                if ($matrix_type eq 'PWM') {
                    $matrix = TFBS::Matrix::PWM->new(
                        -ID           => $id,
                        -name         => $name,
                        -matrixstring => $matrix_string
                    );
                } else {
                    $matrix = TFBS::Matrix::PFM->new(
                        -ID           => $id,
                        -name         => $name,
                        -matrixstring => $matrix_string
                    );
                }


                $matrix_set->add_Matrix($matrix);

                $line_count    = 0;
                $name          = '';
                $matrix_string = '';

                $matrix_count++;
            }
        }
    }

    return $matrix_count ? $matrix_set : undef;
}

1;
