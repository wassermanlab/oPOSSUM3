#!/usr/local/bin/perl -w

=head1 NAME

update_external_gene_ids.pl

=head1 SYNOPSIS

  update_external_gene_ids.pl
      -s species -u opossum_db_user -p opossum_db_pass -ed ensembl_db_name
      [-eh ensembl_db_host] [-d opossum_db_name] [-h opossum_db_name]
      [-t id_type | -all] [-o out_file] [-l log_file]

=head1 ARGUMENTS

  -s species             = Name of species (common name, e.g. 'human').
  -u opossum_db_user     = Name of oPOSSUM DB user.
  -p opossum_db_pass     = oPOSSUM DB password for above user.
  -ed ensembl_db_name    = Name of Ensembl DB from which to extract gene
                           information.
  -eh ensembl_db_host    = Ensembl DB host name.
                           Default = 'vm2.cmmt.ubc.ca'
  -d opossum_db          = Name of oPOSSUM DB to update. If not provided,
                           DB is assumed to be oPOSSUM_2010_<species>.
  -h opossum_db_host     = oPOSSUM DB host name.
                           Default = 'vm5.cmmt.ubc.ca'
  -t gene_it_type        = Type of external gene ID. This must be a
                           number corresponding to the id_type column in
                           the external_gene_id_types table.
  -all                   = Update external gene IDs for all ID types in
                           the external_gene_id_types table.
  -o out_file            = Name of output file to write for later loading
                           into oPOSSUM DB via mysqlimport.
  -l log_file            = Name of log file to which processing and error
                           messages are written.
                           Default = update_external_gene_ids.log

=head1 DESCRIPTION

For each gene in the oPOSSUM DB, retrieve the external gene ID from Ensembl
for the given external gene ID type (or all types in the external_gene_id_types
table. If UPDATE_DB is true, directly update the external_gene_ids table
in the oPOSSUM database, otherwise write the information to a text file for
later loading into the DB using mysqlimport.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

# Ensembl API lib
use lib '/usr/local/src/ensembl-57/ensembl/modules';
use lib '/space/devel/oPOSSUM_2010/lib';
use lib '/raid2/local/src/ensembl-57/ensembl/modules';
use lib '/home/dave/devel/oPOSSUM_2010/lib';


use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use Log::Dispatch::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use OPOSSUM::DBSQL::DBAdaptor;
#use OPOSSUM::Gene;
#use OPOSSUM::Sequence;
#use OPOSSUM::Promoter;
#use OPOSSUM::Exon;

use constant DEBUG                  => 0;
use constant UPDATE_DB              => 0;

# oPOSSUM database settings
use constant OPOSSUM_DB_HOST    	=> 'vm5.cmmt.ubc.ca';

# Local Ensembl database settings
use constant ENSEMBL_DB_HOST    	=> 'vm2.cmmt.ubc.ca';
use constant ENSEMBL_DB_USER    	=> 'ensembl_r';
use constant ENSEMBL_DB_PASS    	=> '';

my $species;
my $opossum_db_name;
my $opossum_db_host;
my $opossum_db_user;
my $opossum_db_pass;
my $ens_db_name;
my $ens_db_host;
my $id_type;
my $all_id_types;
my $out_file;
my $log_file;
GetOptions(
    's=s'	=> \$species,
    'd=s'	=> \$opossum_db_name,
    'h=s'	=> \$opossum_db_host,
    'u=s'	=> \$opossum_db_user,
    'p=s'	=> \$opossum_db_pass,
    'ed=s'	=> \$ens_db_name,
    'eh=s'	=> \$ens_db_host,
    't=i'	=> \$id_type,
    'all'   => \$all_id_types,
    'o=s'	=> \$out_file,
    'l=s'	=> \$log_file
);

$opossum_db_host = OPOSSUM_DB_HOST if !$opossum_db_host;
$ens_db_host = ENSEMBL_DB_HOST if !$ens_db_host;

if (!$species) {
    pod2usage(
        -msg        => "Please specify the species name",
        -verbose    => 1
    );
}

if (!$opossum_db_name) {
    $opossum_db_name = "oPOSSUM_2010_$species";
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

if (!$ens_db_name) {
    pod2usage(
        -msg        => "Please specify the Ensembl DB name",
        -verbose    => 1
    );
}

unless ($id_type || $all_id_types) {
    pod2usage(
        -msg        => "Please specify a gene ID type to update or -all",
        -verbose    => 1
    );
}

#
# -all option overrides individual gene ID type specification
#
if ($id_type && $all_id_types) {
    $id_type = undef;
}

unless (UPDATE_DB || $out_file) {
    pod2usage(
        -msg        => "Please provide an output file name",
        -verbose    => 1
    );
}

if (!$log_file) {
    $log_file = "update_external_gene_ids_$species";
    $log_file .= "_$id_type" if $id_type;
    $log_file .= "_all" if $all_id_types;
    $log_file .= ".log";
}

#
# Initialize logging
#
my $logger = get_logger();
if (DEBUG) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}
my $appender = Log::Log4perl::Appender->new("Log::Dispatch::File",
                    filename    => $log_file,
                    mode        => "write");

#my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("%p: %m%n");

$appender->layout($layout);
$logger->add_appender($appender);

my $start_time = time;
my $localtime = localtime($start_time);

$logger->info("fetch_gene_info started at $localtime\n");

my $opdb = OPOSSUM::DBSQL::DBAdaptor->new(
    -host	    => $opossum_db_host,
    -dbname	    => $opossum_db_name,
    -user	    => $opossum_db_user,
    -password	=> $opossum_db_pass
);

if (!$opdb) {
    $logger->logdie("connecting to oPOSSUM DB $opossum_db_name");
}

my $opxgita = $opdb->get_ExternalGeneIDTypeAdaptor()
    || $logger->logdie("Could not get oPOSSUM ExternalGeneIDTypeAdaptor"); 

my $xgits = $opxgita->fetch_external_gene_id_type_hash()
    || $logger->logdie("Could not fetch external gene ID types");

if ($id_type) {
    my $xgit = $xgits->{$id_type};

    unless ($xgit) {
        $logger->logdie("Unknown extenal gene ID type $id_type");
    }

    my $name = $xgit->name();

    $logger->info("Updating external gene ID type $id_type - $name");
} else {
    $logger->info("Updating all external gene ID types");
    foreach my $id_type (sort keys %$xgits) {
        my $name = $xgits->{$id_type}->name();
        $logger->info("$id_type - $name");
    }
}

my $opga = $opdb->get_GeneAdaptor()
    || $logger->logdie("Could not get oPOSSUM GeneAdaptor"); 

my $gids = $opga->fetch_gene_ids()
    || $logger->logdie("Could not fetch oPOSSUM gene IDs"); 

#
# Get gene ID to Ensembl ID mapping up front
#
my %gid_ensids;
foreach my $gid (@$gids) {
    my $opgene = $opga->fetch_by_gene_id($gid);

    if (!$opgene) {
        $logger->error("Could not fetch oPOSSUM gene $gid"); 
        next;
    }

    my $ensid = $opgene->ensembl_id();

    $gid_ensids{$gid} = $ensid;
}

my $ensdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host	    => $ens_db_host,
    -dbname	    => $ens_db_name,
    -user	    => ENSEMBL_DB_USER,
    -pass	    => ENSEMBL_DB_PASS,
    -species	=> $species,
    -driver	    => 'mysql'
);

if (!$ensdb) {
    $logger->logdie("connecting to Ensembl $species DB $ens_db_name");
}

my $ensga = $ensdb->get_GeneAdaptor
    || $logger->logdie("Could not get Ensembl GeneAdaptor");

my $sth;
if (UPDATE_DB) {
    my $sql = qq{insert into external_gene_ids (gene_id, id_type, external_id)
        values (?,?,?)};

    $sth = $opdb->dbc->prepare($sql)
        || $logger->logdie("Could not prepare insert into external_gene_ids");
}

if ($out_file) {
    open(OFH, ">$out_file")
        || $logger->logdie("Could not open output file $out_file");
}

foreach my $gid (@$gids) {
    my $ensid = $gid_ensids{$gid};

    $logger->info("Processing gene $gid\t$ensid");

    my $ensgene = $ensga->fetch_by_stable_id($ensid);

    if (!$ensgene) {
        $logger->error("fetching Ensembl $species gene $ensid");
        next;
    }

    my $ensextname = $ensgene->external_name();

    #printf "%7d\t%s\t%s:\n", $gid, $ensid, $ensextname;

    foreach my $id_type (sort keys %$xgits) {
        my $name        = $xgits->{$id_type}->name();
        my $dblink_name = $xgits->{$id_type}->dblink_name();

        my $xids = $ensgene->get_all_DBLinks($dblink_name);

        my %inc;
        foreach my $xid (@$xids) {
            my $xdbname = $xid->dbname();
            my $id      = $xid->primary_id;
            my $symbol  = $xid->display_id;

            my $external_id;
            if ($id_type == 1) {
                $external_id = $symbol;
            } else {
                $external_id = $id;
            }

            unless ($inc{$external_id}) {
                #printf "%23s\t%15s\t%15s\t%15s\n",
                #    $name, $xdbname, $id, $symbol;

                #printf "%23s\t%15s\n",
                #    $name, $external_id;

                $logger->info("inserting $name\t$external_id\n");

                if (UPDATE_DB) {
                    $sth->execute($gid, $id_type, $external_id)
                        || $logger->logdie(
                            "Could not execute insert into external_gene_ids"
                        );
                } else {
                    print OFH "$gid\t$id_type\t$external_id\n";
                }

                $inc{$external_id} = 1;
            }
        }
    }
    #print "\n";
}

if ($sth) {
    $sth->finish();
}

if ($out_file) {
    close(OFH);
}

my $end_time = time;
$localtime = localtime($end_time);
my $elapsed_secs = $end_time - $start_time;

$logger->info("update_external_gene_ids completed at $localtime\n");
$logger->info("Elapsed time (s): $elapsed_secs");

exit;
