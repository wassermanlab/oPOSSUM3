#!/usr/local/bin/perl -w

=head1 NAME

  fetch_tf_info.pl [-f matrix_file] [-d jaspar_db_name] [-h jaspar_db_host]
                   [-c collection] [-s source] [-t tax_group] [-pwm]
                   -o out_file

=head1 SYNOPSIS

  fetch_tf_info.pl

=head1 ARGUMENTS

  -d jaspar_db_name If specified, use this as the JASPAR database from
                    which to select TFs.
  -h jaspar_db_host If specified, use this as the JASPAR database host.
  -s source         This value is entered in the source field in the
                    tf_info table (if not specified and -d option is, use
                    jaspar_db_name).
  -t tax_group      If specified, limit TFs to this taxonomic supergroup
                    (phylum/species).
  -f matrix_file    If specified use this PFM file to create from,
                    instead of the JASPAR database
  -pwm              If -f is specified, input matrix file contains PWMs
                    (default is PFMs)
  -o out_file       Output file cotaining TF data formatted for loading
                    into oPOSSUM tf_info table via mysqlimport

=head1 DESCRIPTION

Create a text file from which to populate the tf_info table of the oPOSSUM
database from the various JASPAR databases or an input PFM (or PWM) file.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use lib '/space/devel/oPOSSUM_2010/lib';

use Getopt::Long;
use Pod::Usage;
use TFBS::DB::JASPAR5;
use OPOSSUM::TFInfo;

use constant TFBS_MATRIXTYPE    => 'PFM';     # for JASPAR matrices
use constant TFBS_MIN_IC        => 8;

use constant SOURCE             => 'JASPAR';

use constant JASPAR_DB_NAME     => 'JASPAR_2010';
use constant JASPAR_DB_HOST     => "vm5.cmmt.ubc.ca";
use constant JASPAR_DB_USER     => "jaspar_r";
use constant JASPAR_DB_PASS     => "";

use constant JASPAR_COLLECTIONS  => (
    'CORE',
    'CNE',
    'PHYLOFACTS',
    'SPLICE',
    'POLII',
    'FAM',
    'PBM',
    'PBM_HOMEO',
    'PBM_HLH'
);

my $jaspar_db_name;
my $jaspar_db_host;
my $collection;
my $tax_group;
my $matrix_file;
my $pwm_flag = 0;
my $source;
my $out_file;
GetOptions(
    'd=s' => \$jaspar_db_name,
    'h=s' => \$jaspar_db_host,
    'c=s' => \$collection,
    't=s' => \$tax_group,
    'f=s' => \$matrix_file,
    'pwm' => \$pwm_flag,
    's=s' => \$source,
    'o=s' => \$out_file
);

if ($matrix_file && !$source) {
    pod2usage(
        -msg        => "If specifying a matrix file, please also specify"
                       . " a TFBS profile source, e.g. 'CUSTOM'",
        -verbose    => 1
    );
}

if (!$out_file) {
    pod2usage(
        -msg        => "No output file specified",
        -verbose    => 1
    );
}

$jaspar_db_name = JASPAR_DB_NAME if !$jaspar_db_name && !$matrix_file;
$jaspar_db_host = JASPAR_DB_HOST if !$jaspar_db_host && !$matrix_file;

$source = $jaspar_db_name if !$source && $jaspar_db_name;
$source = SOURCE if !$source;

my $tf_id = 1;
my @tf_info_list;
if ($matrix_file) {
    my $matrix_set;
    if ($pwm_flag) {
        $matrix_set = read_PWMs($matrix_file);
    } else {
        $matrix_set = read_PFMs($matrix_file);
    }

    die "Error reading matrix set from file $matrix_file\n"
        if !$matrix_set;

    my $iter = $matrix_set->Iterator(-sort_by => 'name');
    while (my $matrix = $iter->next) {
        my $total_ic = 0;
        if (!$pwm_flag) {
            $total_ic = $matrix->to_ICM->total_ic;
        }

        my $tf_info = OPOSSUM::TFInfo->new(
            -id             => $tf_id++,
            -source         => $source,
            -collection     => $collection,
            -external_id    => $matrix->ID,
            -name           => $matrix->name,
            -class          => $matrix->class,
            -family         => $matrix->{tags}{'family'},
            -tax_group      => $matrix->{tags}{'tax_group'},
            -width          => $matrix->length,
            -ic             => $total_ic
        );

        push @tf_info_list, $tf_info;
    }
} else {
    my $jdb = TFBS::DB::JASPAR5->connect(
        "dbi:mysql:$jaspar_db_name:$jaspar_db_host",
        JASPAR_DB_USER,
        JASPAR_DB_PASS
    );

    die "Error connecting to JASPAR database $jaspar_db_name\n" if !$jdb;

    my @collections;
    if ($collection) {
        push @collections, $collection;
    } else {
        @collections = JASPAR_COLLECTIONS;
    }

    my %matrix_args;
    #$matrix_args{-min_ic}       = TFBS_MIN_IC; # don't limit
    $matrix_args{-taxgroup}     = $tax_group if $tax_group;
    $matrix_args{-collection}   = \@collections;

    my $matrix_set = $jdb->get_MatrixSet(%matrix_args);
    die "No matrices fetched from JASPAR DB $jaspar_db_name" if !$matrix_set;

    my $iter = $matrix_set->Iterator(-sort_by => 'ID');
    while (my $matrix = $iter->next) {
        my $tf_info = OPOSSUM::TFInfo->new(
            -id             => $tf_id++,
            -source         => $source,
            -collection     => $matrix->{tags}{'collection'},
            -external_id    => $matrix->ID,
            -name           => $matrix->name,
            -class          => $matrix->class,
            -family         => $matrix->{tags}{'family'},
            -tax_group      => $matrix->{tags}{'tax_group'},
            -width          => $matrix->length,
            -ic             => $matrix->to_ICM->total_ic
        );

        push @tf_info_list, $tf_info;
    }
}

open(OFH, ">$out_file") || die "Error opening output file $out_file - $!\n";

foreach my $tf_info (@tf_info_list) {
    printf OFH "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.3f\n",
        $tf_info->id,
        $tf_info->source,
        $tf_info->collection,
        $tf_info->external_id,
        $tf_info->name,
        $tf_info->class,
        $tf_info->family,
        $tf_info->tax_group,
        $tf_info->width,
        $tf_info->ic;
}

close(OFH);

exit;

sub read_PFMs
{
    my ($file) = @_;

    open(FH, $file) || die("Error opening matrix file $file\n");

    my $matrix_set = TFBS::MatrixSet->new();

    my $name = '';
    my $id   = '';
    $tax_group = 'human' if !$tax_group;
    my $matrix_string = '';
    my $line_count    = 0;
    while (my $line = <FH>) {
        chomp $line;
        next if !$line;
        if ($line =~ /^>\s*(.+)/) {
            my $tagline = $1;
            # very specific code for Jenn Gardy's (hancock) matrix set!!!
            if ($tagline =~ /(\S+)\s+\[JASPAR\s+(\S+)\]/) {
                $name     = $1;
                $id       = $1;
                $tax_group = $2;
                $tax_group = '' if $tax_group eq 'unspecified';
            } elsif ($tagline =~ /(\S+)\s+\|\s+(\S+)/) {
                $name = $1;
                $id   = $1;
            } else {
                #warn "Unrecognized matrix tagline: $tagline\n";
                $name = $tagline;
                $id   = $tagline;
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
                my $pfm = TFBS::Matrix::PFM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $id || '',
                    -class        => '',
                    -tags         => ({tax_group => $tax_group})
                );

                $matrix_set->add_Matrix($pfm);

                $line_count    = 0;
                $name          = '';
                $id            = '';
                $matrix_string = '';
            }
        }
    }
    close(FH);

    return $matrix_set;
}

sub read_PWMs
{
    my ($file) = @_;

    open(FH, $file) || die("Error opening matrix file $file\n");

    my $matrix_set = TFBS::MatrixSet->new();

    my $name = '';
    my $id   = '';
    $tax_group = 'human' if !$tax_group;
    my $matrix_string = '';
    my $line_count    = 0;
    while (my $line = <FH>) {
        chomp $line;
        if ($line =~ /^>\s*(.+)/) {
            my $tagline = $1;
            # very specific code for Jenn Gardy's (hancock) matrix set!!!
            if ($tagline =~ /(\S+)\s+\[JASPAR\s+(\S+)\]/) {
                $name     = $1;
                $id       = $1;
                $tax_group = $2;
                $tax_group = '' if $tax_group eq 'unspecified';
            } elsif ($tagline =~ /(\S+)\s+\|\s+(\S+)/) {
                $name = $1;
                $id   = $1;
            } else {
                #warn "Unrecognized matrix tagline: $tagline\n";
                $name = $tagline;
                $id   = $tagline;
            }
        } else {
            $line_count++;
            if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            } else {
                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            }
            if ($line_count == 4) {
                my $pwm = TFBS::Matrix::PWM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $id || '',
                    -class        => '',
                    -tags         => ({tax_group => $tax_group})
                );

                $matrix_set->add_Matrix($pwm);

                $line_count    = 0;
                $name          = '';
                $id            = '';
                $matrix_string = '';
            }
        }
    }
    close(FH);

    return $matrix_set;
}
