#!/usr/local/bin/perl -w

=head1 NAME

populate_db_info.pl

=head1 SYNOPSIS

  populate_db_info.pl
      -h opossum_db_host -d opossum_db_name -u opossum_db_user
      -p opossum_db_password
      -s species_name -ln latin_name -a assembly -edb ensembl_db_name
      -udb ucsc_db_name -uct cons_table_name -t min_threshold
      -ubp max_upstream_bp -crl min_cr_length -ic min_ic

=head1 ARGUMENTS

  -h opossum_db_host     = Host name of oPOSSUM DB.
  -d opossum_db_name     = Name of oPOSSUM DB.
  -u opossum_db_user     = Name of oPOSSUM DB user (write access).
  -p opossum_db_password = oPOSSUM database user password.
  -s species_name        = Species common name, e.g. 'human'.
  -ln latin_name         = Species latin name, e.g. 'homo sapiens'.
  -a assembly            = NCBI genome assembly, e.g. 'NCBI37'.
  -edb ensembl_db_name   = Name of Ensembl DB on which gene information
                           is based.
  -udb ucsc_db_name      = Name of UCSC DB on which phastCons conservation
                           is based.
  -uct cons_table_name   = Name of UCSC DB table from which phastCons
                           conservation scores are extracted.
  -t min_threshold       = The absolute minimum PWM threshold used to
                           compute TFBSs. This must be less than or equal
                           to the threshold value at threshold level 1 in
                           the threshold_levels table.
  -ubp max_upstream_bp   = The maximum amount of upstream bp to include in
                           the TFBS searching. This must be at least as
                           much as the upstream_bp value at level 1 in
                           the search_region_levels table.
  -crl min_cr_length     = The minimum length of conserved regions to use
                           for TFBS searching.
  -ic min_ic             = The minimum information content of JASPAR TFBS
                           profile matrices to use in TFBS searching.

=head1 DESCRIPTION

Populate the oPOSSUM db_info table with the provided values.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use lib '/apps/oPOSSUM3/lib';


use Getopt::Long;
use Pod::Usage;
use OPOSSUM::DBInfo;
use OPOSSUM::DBSQL::DBAdaptor;


my $opossum_db_host;
my $opossum_db_name;
my $opossum_db_user;
my $opossum_db_pass;
my $species;
my $latin_name;
my $assembly;
my $ens_db_name;
my $ucsc_db_name;
my $cons_table_name;
my $min_threshold;
my $max_upstream_bp;
my $min_cr_length;
my $min_ic;
GetOptions(
    'h=s'   => \$opossum_db_host,
    'd=s'   => \$opossum_db_name,
    'u=s'   => \$opossum_db_user,
    'p=s'   => \$opossum_db_pass,
    's=s'   => \$species,
    'ln=s'  => \$latin_name,
    'a=s'   => \$assembly,
    'edb=s' => \$ens_db_name,
    'udb=s' => \$ucsc_db_name,
    'uct=s' => \$cons_table_name,
    't=s'   => \$min_threshold,
    'ubp=i' => \$max_upstream_bp,
    'crl=i' => \$min_cr_length,
    'ic=i'  => \$min_ic
);

if (!$opossum_db_host) {
    pod2usage(
        -msg        => "Please specify the oPOSSUM DB host",
        -verbose    => 1
    );
}

if (!$opossum_db_name) {
    pod2usage(
        -msg        => "Please specify the oPOSSUM DB name",
        -verbose    => 1
    );
}

if (!$opossum_db_user) {
    pod2usage(
        -msg        => "Please specify the oPOSSUM DB user name",
        -verbose    => 1
    );
}

if (!$opossum_db_pass) {
    pod2usage(
        -msg        => "Please specify the oPOSSUM DB password",
        -verbose    => 1
    );
}

if (!$species) {
    pod2usage(
        -msg        => "Please specify the species name",
        -verbose    => 1
    );
}

if (!$latin_name) {
    pod2usage(
        -msg        => "Please specify the species latin name",
        -verbose    => 1
    );
}

if (!$assembly) {
    pod2usage(
        -msg        => "Please specify the genome assembly",
        -verbose    => 1
    );
}

if (!$ens_db_name) {
    pod2usage(
        -msg        => "Please specify the Ensembl DB name",
        -verbose    => 1
    );
}

if (!$ucsc_db_name) {
    pod2usage(
        -msg        => "Please specify the UCSC DB name",
        -verbose    => 1
    );
}

if (!$cons_table_name) {
    pod2usage(
        -msg        => "Please specify the UCSC phastCons table name",
        -verbose    => 1
    );
}

if (!$min_threshold) {
    pod2usage(
        -msg        => "Please specify the minimum PWM score threshold",
        -verbose    => 1
    );
}

if (!defined $max_upstream_bp) {
    pod2usage(
        -msg        => "Please specify the maximum upstream bp",
        -verbose    => 1
    );
}

if (!$min_cr_length) {
    pod2usage(
        -msg        => "Please specify the minimum conserved region length",
        -verbose    => 1
    );
}

if (!defined $min_ic) {
    pod2usage(
        -msg        => "Please specify the minimum IC",
        -verbose    => 1
    );
}

my $db_info = OPOSSUM::DBInfo->new(
    -species         => $species,
    -latin_name      => $latin_name,
    -assembly        => $assembly,
    -ensembl_db      => $ens_db_name,
    -ucsc_db         => $ucsc_db_name,
    -ucsc_cons_table => $cons_table_name,
    -min_threshold   => $min_threshold,
    -max_upstream_bp => $max_upstream_bp,
    -min_cr_length   => $min_cr_length,
    -min_ic          => $min_ic
);

my $dba = OPOSSUM::DBSQL::DBAdaptor->new(
    -host       => $opossum_db_host,
    -dbname     => $opossum_db_name,
    -user       => $opossum_db_user,
    -password   => $opossum_db_pass
) || die "Error connecting to oPOSSUM DB\n";

my $dbia = $dba->get_DBInfoAdaptor() || die "Error getting DBInfoAdaptor\n";

$dbia->store($db_info) || die "Error storing db_info record\n";
