#!/usr/local/bin/perl -w

#
# Foreach profile in the JASPAR database, compute the GC content and write
# out SQL statements to add a tag/value par for insertion into the
# MATRIX_ANNOTATION table.
#

use TFBS::DB::JASPAR5;

use constant JASPAR_DB_NAME => 'JASPAR_2010';
use constant JASPAR_DB_HOST => 'vm5.cmmt.ubc.ca';
use constant JASPAR_DB_USER => 'jaspar_r';
use constant JASPAR_DB_PASS => '';

use strict;

my $jdb = TFBS::DB::JASPAR5->connect(
    "dbi:mysql:" . JASPAR_DB_NAME . ":" . JASPAR_DB_HOST,
    JASPAR_DB_USER,
    JASPAR_DB_PASS
);

if (!$jdb) {
    die "Could not connect to JASPAR database";
}

my $matrix_set = $jdb->get_MatrixSet(
    -all            => 1,
    -matrixtype     => 'PFM'
);

unless ($matrix_set) {
    die "Error fetching matrices from JASPAR";
}

my $iter = $matrix_set->Iterator();

while (my $matrix = $iter->next) {
    my ($base_id, $version) = split (/\./,  $matrix->ID);

    my $internal_id = $jdb->_get_internal_id($base_id, $version);

    my $gc_content = compute_matrix_gc_content($matrix);

    printf("insert into MATRIX_ANNOTATION (ID, TAG, VAL) values (%d, 'gc_content', %.3f);\n", $internal_id, $gc_content);
}

exit;

sub compute_matrix_gc_content
{
    my $pfm = shift;

    my $matrix = $pfm->matrix();

    my $gc_count = 0;
    my $total_count = 0;
    my $row_num = 0;
    foreach my $row (@$matrix) {
        $row_num++;
        foreach my $val (@$row) {
            if ($row_num == 2 || $row_num == 3) {
                $gc_count += $val;
            }

            $total_count += $val;
        }
    }

    return $gc_count / $total_count;
}
