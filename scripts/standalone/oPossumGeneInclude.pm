
=head1 NAME

 oPossumGeneInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the gene-based analysis
  scripts.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

use oPossumOpt;

# Form selection values and default settings
use constant DFLT_BG_NUM_RAND_GENES     => 5000;

#
# oPOSSUM DB Access
#
# The species name is dynamically appended to OPOSSUM_DB_NAME for full
# oPOSSUM DB name, e.g. oPOSSUM_2010_human
#
use constant OPOSSUM_DB_HOST    => 'vm5.cmmt.ubc.ca';
use constant OPOSSUM_DB_NAME    => 'oPOSSUM_2010';
use constant OPOSSUM_DB_USER    => 'opossum_r';
use constant OPOSSUM_DB_PASS    => '';

#
# Form selection values and default settings. These may need to be overriden 
# by specific oPOSSUM variants.
#
use constant DFLT_GENE_ID_TYPE          => 0;
use constant DFLT_CONSERVATION_LEVEL    => 3;
use constant DFLT_THRESHOLD_LEVEL       => 2;
use constant DFLT_SEARCH_REGION_LEVEL   => 3;

# Leave as undef or set to 'all' for all biotypes (not tested).
use constant DFLT_BIOTYPE               => undef;


use lib OPOSSUM_LIB_PATH;

use Getopt::Long;
use Pod::Usage;
use File::Temp;
use Carp;
#use CGI::Carp qw(carpout);
#use Template;
use File::Temp qw/ tempfile tempdir /;
use Log::Log4perl qw(get_logger :levels);

use Bio::SeqIO;

use TFBS::DB::JASPAR5;

use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::Analysis::Counts;


sub opossum_db_connect
{
    my ($species) = @_;
    
    my $db_name = sprintf("%s_%s", OPOSSUM_DB_NAME, $species);
    
    my $opdba = OPOSSUM::DBSQL::DBAdaptor->new(
        -host     => OPOSSUM_DB_HOST,
        -dbname   => $db_name,
        -user     => OPOSSUM_DB_USER,
        -password => OPOSSUM_DB_PASS
    );
    
    return $opdba;
}

#
# Connect to JASPAR database
#
sub jaspar_db_connect
{
    my ($tf_db) = @_;
    
    my $jdb = TFBS::DB::JASPAR5->connect(
        "dbi:mysql:" . $tf_db . ":" . JASPAR_DB_HOST,
        JASPAR_DB_USER,
        JASPAR_DB_PASS
    );

    return $jdb;
}


#
# For the provided list of Ensembl / external gene IDs, fetch the associated
# oPOSSUM gene IDs from the DB. Also keep track of which provided gene IDs
# mapped (included) or did not map (missing) to oPOSSUM gene IDs.
#
# NOTE: There can be a 1-to-many mapping of oPOSSUM gene IDs to external gene
# IDs. Make sure all of these are retained in the included external gene IDs
# list.
#
# Andrew's NOTE: This function now takes an additional argument $has_operon
# This determines whether the operon information should be used when fetching
# the genes.
#
# If operons are used, then for each gene which is in an operon, store the
# mapping of the genes ID to the ID of the first gene in the operon.
#
# Also return an additional array which contains the actual gene IDs to use
# in fetching TFBSs, TFBS counts, conserved regions lengths etc. This list
# will contain:
# - for a gene which is NOT in an operon, it's own gene ID
# - for a gene which is in an operon, the gene ID of the first gene in the
#   operon
#
sub fetch_opossum_gene_ids
{
    my ($ga, $oa, $id_type, $ids, $has_operon) = @_;

    my @mapped_ids;
    my @unmapped_ids;
    my %gene_id_ext_ids;

    #
    # This will be set to a mapping of operon gene IDs to the first gene ID
    # in the operon.
    #
    my %operon_first_gene_id;

    #
    # This will be set to a list of gene IDs for both genes which do not
    # belong to an operon and genes which are the first gene in an operon 
    # that one or more of the input genes belongs to.
    #
    my @operon_unique_gene_ids;

    #
    # These are never used
    #
    #my %ext_id_gene_ids;
    #my %gene_id_ensembl_ids; 

    my $gene_ids;
    if (defined $id_type && $id_type == DFLT_GENE_ID_TYPE) {
        #
        # Ensembl IDs provided.
        # 
        $gene_ids = $ga->fetch_gene_ids_by_ensembl_ids(
            $ids,
            -mapped_ensembl_ids     => \@mapped_ids,
            -unmapped_ensembl_ids   => \@unmapped_ids,
            #-gene_id_ensembl_ids    => \%gene_id_ensembl_ids
            #-ensembl_id_gene_ids    => \%ext_id_gene_ids
        );
        #%gene_id_ext_ids = %gene_id_ensembl_ids;
    } else { 
        #
        # Fetch by external gene IDs
        # 
        $gene_ids = $ga->fetch_gene_ids_by_external_ids(
            $ids,
            $id_type,
            -mapped_ext_ids         => \@mapped_ids,
            -unmapped_ext_ids       => \@unmapped_ids,
            -gene_id_ext_ids        => \%gene_id_ext_ids
            #-ext_id_gene_ids        => \%ext_id_gene_ids
            #-gene_id_ensembl_ids    => \%gene_id_ensembl_ids
        );
    } 

    if ($has_operon) {
        #my $oa = $opdba->get_OperonAdaptor();

        my %fgene_id_inc;
        foreach my $gene_id (@$gene_ids) {
            # first check to see if this gene belongs to an operon
            # it is possible to have multiple genes belonging to the same
            # operon to be included in the genes to be fetched.
            # map each gene id in the operon to the first gene id
            my $operon = $oa->fetch_by_gene_id($gene_id);

            if ($operon) {
                #print STDERR "Operon " . $operon->id . " for gene $gene_id\n";
                my $fgene = $operon->fetch_first_gene();

                $operon_first_gene_id{$gene_id} = $fgene->id;

                unless ($fgene_id_inc{$fgene->id}) {
                    push @operon_unique_gene_ids, $fgene->id;
                }
            } else {
                push @operon_unique_gene_ids, $gene_id;
            }
        }
    } 

    return (
        $gene_ids,
        @mapped_ids                 ? \@mapped_ids : undef,
        @unmapped_ids               ? \@unmapped_ids : undef,
        #%gene_id_ensembl_ids        ? \%gene_id_ensembl_ids : undef,
        %gene_id_ext_ids            ? \%gene_id_ext_ids : undef,
        #%ext_id_gene_ids            ? \%ext_id_gene_ids : undef,
        #%operon_first_gene_to_gene_ids
        #                            ? \%operon_first_gene_to_gene_ids : undef
        %operon_first_gene_id       ? \%operon_first_gene_id : undef,
        @operon_unique_gene_ids     ? \@operon_unique_gene_ids : undef
    );
}

sub fetch_all_gene_ids
{
    # AK: added biotype parameter
    my ($ga, $biotype) = @_;
    
    my $where = _biotype_where_clause($biotype);
    
    my $gids = $ga->fetch_gene_ids($where);

    return $gids;
}

sub fetch_random_opossum_gene_ids
{
    my ($ga, $oa, $num_genes, $has_operon, $biotype) = @_;

    return if !$num_genes;

    $biotype = _parse_biotype($biotype) if $biotype;

    my $where = _biotype_where_clause($biotype);

    my $total_num_genes = $ga->fetch_gene_count($where);

    $num_genes = $total_num_genes if $num_genes > $total_num_genes;

    my %fetch_args = (-num_genes => $num_genes);

    if ($biotype) {
        if (ref $biotype eq 'ARRAY') {
            $fetch_args{-biotypes} = $biotype;
        } else {
            $fetch_args{-biotype} = $biotype;
        }
    }

    my $rand_gids = $ga->fetch_random_gene_ids(%fetch_args);

    #if ($has_operon) {
    #    my %incl_operon;
    #    my %incl_gid;
    #    my @passed_gids;

    #    my $oa = $opdba->get_OperonAdaptor();
    #    unless ($oa) {
    #        $self->_error("Error getting OperonAdaptor");
    #        return;
    #    }

    #    # so that we can skip them if they turn up again
    #    foreach my $rand_gid (@$rand_gids) {
    #        $incl_gid{$rand_gid} = 1;
    #    }

    #    my $gid_opid_map = $oa->fetch_gene_ids_to_operon_ids_hash($rand_gids);

    #    foreach my $gid (@$rand_gids) {
    #        my $operon = $$gid_opid_map{$gid};
    #        if (!$operon) {
    #            push @passed_gids, $gid;
    #            next;
    #        }
    #        if (!$incl_operon{$operon->id}) {
    #            $incl_operon{$operon->id} = $gid;
    #            push @passed_gids, $gid;
    #        }
    #    }

    #    # just in case num_genes > scalar(@$rand_gids), don't use num_genes
    #    while (scalar(@passed_gids) < scalar(@$rand_gids)) {
    #        my $new_gid;
    #        if (@biotypes) {
    #            $new_gid = $ga->fetch_random_gene_ids(
    #                -num_genes => 1,
    #                -biotypes   => \@biotypes
    #            );
    #        } else {
    #            $new_gid = $ga->fetch_random_gene_ids(
    #                -num_genes => 1,
    #                -biotype   => $biotype
    #            );
    #        }
    #        next if $incl_gid{$new_gid};
    #        my $op = $oa->fetch_by_gene_id($$new_gid[0]);
    #        next if $incl_operon{$op->id};
    #        push @passed_gids, $$new_gid[0];
    #    }

    #    $rand_gids = \@passed_gids;
    #}

    my %operon_first_gene_id;
    my @operon_unique_gene_ids;

    if ($has_operon) {
        # if dflt_gene_id_type, # genes = 1
        my %fgene_id_inc;
        foreach my $gid (@$rand_gids) {
            # first check to see if this gene belongs to an operon
            # it is possible to have multiple genes belonging to the same operon
            # to be included in the genes to be fetched.
            # thus need to keep an array of gene ids
            # map each non-first gene id in the operon to the first gene id
            my $operon = $oa->fetch_by_gene_id($gid);
            if ($operon) {
                my $fgene = $operon->fetch_first_gene();
                #if ($$rand_gids[$i] != $fgene->id) {
                    $operon_first_gene_id{$gid} = $fgene->id;
                #}

                unless ($fgene_id_inc{$fgene->id}) {
                    push @operon_unique_gene_ids, $fgene->id;
                }
            } else {
                push @operon_unique_gene_ids, $gid;
            }
        }
    }

    return (
        $rand_gids,
        %operon_first_gene_id   ? \%operon_first_gene_id : undef,
        @operon_unique_gene_ids ? \@operon_unique_gene_ids : undef
    );
} 


sub read_ids_file
{
    my ($file, $logger) = @_;

    open(FH, $file) || $logger->logdie("Could not open IDs file $file - $!");

    my @ids;
    while (my $line = <FH>) {
        chomp $line;

        push @ids, $line;
    }

    close(FH);

    return @ids ? \@ids : undef;
}

sub fetch_cr_gc_content
{
    my ($ga, $cra, $gids, $clevel, $upstream_bp, $downstream_bp, $biotype) = @_;
    
    if (!$gids or scalar(@$gids) == 0) {
        # all gene ids
        my $where = _biotype_where_clause($biotype);

        $gids = $ga->fetch_gene_ids($where);
    }
        
    my $sum_gc = 0;
    my $sum_length = 0;
    foreach my $gid (@$gids) {
        my $gc_content = $cra->fetch_gc_content_by_upstream_downstream(
            $gid, $clevel, $upstream_bp, $downstream_bp
        );

        my $cr_length = $cra->fetch_length_by_upstream_downstream(
            $gid, $clevel, $upstream_bp, $downstream_bp
        );
        
        $sum_gc += $gc_content * $cr_length;
        $sum_length += $cr_length;
    }
    
    #print STDERR "fetch_cr_gc_content: # gids = " . scalar(@$gids) . "\n";
    #print STDERR "fetch_cr_gc_content: sum(gc_content) = $sum_gc_content\n";
    my $avg_gc_content = $sum_gc / $sum_length;
    $avg_gc_content = sprintf("%.2f", $avg_gc_content);
    
    return $avg_gc_content;
}


#
# Biotype can be a single value, a comma separated string or an array ref
#
sub _biotype_where_clause
{
    my ($biotype) = @_;

    return unless $biotype;

    $biotype = _parse_biotype($biotype);

    my $where;
    if ($biotype) {
        if (ref $biotype eq 'ARRAY') {
            $where = "where biotype in ('"
                . join("','", @$biotype)
                . "')";
        } elsif ($biotype !~ /^all$/i) {
            # don't set if biotype is 'all'
            $where = "where biotype = '$biotype'";
        }
    }

    return $where;
}

#
# Parse biotype argument. Biotype may be a single value, a multi-value string
# separated by spaces or commas or an array ref. If it is a multi-value string, # convert it to an array and return an array ref. If it is a single value and
# it is equal to 'all' return as undef. otherwise return passed argument
# as is.
#
sub _parse_biotype
{
    my $biotype = shift;

    return unless $biotype;

    return $biotype if ref $biotype eq 'ARRAY';

    $biotype =~ s/^\s+//;   # remove leading spaces
    $biotype =~ s/\s+$//;   # remove trailing spaces

    if ($biotype =~ /,/) {
        my @biotypes = split /\s*,\s*/, $biotype;

        return \@biotypes;
    } elsif ($biotype =~ /\s+/) {
        my @biotypes = split /\s+/, $biotype;

        return \@biotypes;
    } elsif ($biotype =~ /^all$/i) {
        return undef;
    }

    return $biotype;
}

1;