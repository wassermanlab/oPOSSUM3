#
# This module should be included in all the oPossumSeq*Web.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to the sequence-based oPossum variants.
#

use OPOSSUM::Web::Opt::BaseOpt;
use OPOSSUM::Web::Opt::SeqOpt;

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

        unless ($self->check_seq_text($sl)) {
            return;
        }

        $filename = "$tempdir/seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target sequence file $filename");
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
            $self->_error("Target sequence upload file is empty");
            return;
        }

        unless ($self->check_seq_text($sl)) {
            return;
        }

        $filename = "$tempdir/seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target sequence file $filename");
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

        unless ($self->check_seq_text($sl)) {
            return;
        }

        $filename = "$tempdir/back_seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error(
                "Unable to create background sequence file $filename"
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
            $self->_error("Background sequence upload file is empty");
            return;
        }

        unless ($self->check_seq_text($sl)) {
            return;
        }

        $filename = "$tempdir/back_seqs.fa";
        unless (open(FH, ">$filename")) {
            $self->_error(
                "Unable to create background sequence file $filename"
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

sub get_peak_pos_file
{
    my ($self, $peak_pos_input_method, $tempdir) = @_;

    my $q = $self->query;

    my $filename;
    if ($peak_pos_input_method eq "paste") {
        my $sl = $q->param('peak_pos_list');
        if (!$sl) {
            $self->_error("Sequence input not specified");
            return;
        }

        $filename = "$tempdir/peak_pos.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target sequence peak position file $filename");
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($peak_pos_input_method eq "upload") {
        my $file = $q->param('peak_pos_file');
        my $fh   = $q->upload('peak_pos_file');

        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("Target peak position file is empty");
            return;
        }

        $filename = "$tempdir/peak_pos.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target sequence peak position file $filename");
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown peak position input method");
        return;
    }

    return $filename;
}

sub get_bg_peak_pos_file
{
    my ($self, $peak_pos_input_method, $tempdir) = @_;

    my $q = $self->query;

    my $filename;
    if ($peak_pos_input_method eq "paste") {
        my $sl = $q->param('bg_peak_pos_list');
        if (!$sl) {
            $self->_error("Background sequence input not specified");
            return;
        }

        $filename = "$tempdir/bg_peak_pos.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create background sequence peak position file $filename");
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($peak_pos_input_method eq "upload") {
        my $file = $q->param('bg_peak_pos_file');
        my $fh   = $q->upload('bg_peak_pos_file');

        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("Background peak position upload file is empty");
            return;
        }

        $filename = "$tempdir/bg_peak_pos.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create background sequence peak position file $filename");
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown peak position input method");
        return;
    }

    return $filename;
}


sub write_matrix_file
{
    my ($self, $matrix_spec, $tempdir, $fname) = @_;

    my $filename;
    
    if ($fname) {
        $filename = "$tempdir/$fname";
    } else {
        $filename = "$tempdir/matrices.txt";
    }

    unless (open(FH, ">$filename")) {
        $self->_error("Unable to create matrix file $filename");
        return;
    }
    
    if ($matrix_spec->isa('TFBS::MatrixSet')) {
        my $iter = $matrix_spec->Iterator();
        while (my $matrix = $iter->next()) {
            printf FH ">%s %s\n", $matrix->ID(), $matrix->name();
            printf FH "%s\n", $matrix->prettyprint();
        }
    } elsif ($matrix_spec->isa('TFBS::Matrix')) {
        printf FH ">%s %s\n", $matrix_spec->ID(), $matrix_spec->name();
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
    my $id;
    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    my $matrix_count  = 0;
    my @lines         = split /\n/, $text;
    foreach my $line (@lines) {
        chomp $line;
        next if !$line;

        if ($line =~ /^>\s*(.*)/) {
            my $desc = $1;

            if ($desc) {
                ($id, $name) = split /\s+/, $desc;
            }

            unless ($id) {
                $id = sprintf "matrix%d", $matrix_count + 1;
            }

            unless ($name) {
                $name = $id;
            }
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

sub check_seq_text
{
    my ($self, $seq_text) = @_;

    my $seq_text_len = length $seq_text;
    #carp "Sequence text length is $seq_text_len\n";
    if ($seq_text_len > MAX_SEQ_FILE_SIZE) {
        $self->_error(
            "Total sequence size exceeds maximum allowed size of "
            .  MAX_SEQ_FILE_SIZE . " bytes"
        );
        return;
    }

    my @lines = split /\n/, $seq_text;
    foreach my $line (@lines) {
        chomp $line;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        next unless $line;
        unless ($line =~ /^\s*>/) {
            if ($line =~ /[^actgnxACTGNX]/) {
                $self->_error(
                      "Sequences contain invalid characters. Please make sure"
                      . " the sequences are in plain text fasta format"
                );
                return;
            }
        }
    }

    return 1;
}

1;
