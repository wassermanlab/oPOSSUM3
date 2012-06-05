#
# This module should be included in all the oPossumGene*Web.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to all the oPossum gene-based variants.
#

use OPOSSUM::Web::Opt::BaseOpt;
use OPOSSUM::Opt::GeneOpt;

use lib OPOSSUM_LIB;

use OPOSSUM::DBSQL::DBAdaptor;

use Data::Dumper;    # for debugging only

use strict;

sub opossum_db_connect
{
    my ($self) = @_;

    my $species = $self->state->species();

    unless ($species) {
        $self->_error("Species not set");
        return;
    }

    my $db_name = sprintf("%s_%s", OPOSSUM_DB_NAME, $species);

    my $dba = OPOSSUM::DBSQL::DBAdaptor->new(
        -host     => OPOSSUM_DB_HOST,
        -dbname   => $db_name,
        -user     => OPOSSUM_DB_USER,
        -password => OPOSSUM_DB_PASS
    );

    unless ($dba) {
        $self->_error("Could not connect to oPOSSUM database $db_name");
        return;
    }

    $self->opdba($dba);
}

sub opdba
{
    my $self = shift;

    if (@_) {
        $self->{-opdba} = shift;
    }

    return $self->{-opdba};
}

sub fetch_all_gene_ids
{
    # AK: added biotype parameter
    my ($self, $biotype) = @_;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $ga = $opdba->get_GeneAdaptor();
    if (!$ga) {
        $self->_error("Could not get GeneAdaptor");
    }
    
    # let's use this to skip over pseudogenes
    #my $where = _biotype_where_clause($biotype);
    my $where = "biotype != 'pseudogene'";
    
    my $gids = $ga->fetch_gene_ids($where);
    if (!$gids) {
        $self->_error("Could not fetch oPOSSUM gene IDs");
    }

    #$self->state->db_info($db_info);

    return $gids;
}

sub fetch_gene_count
{
    # AK: added where parameter
    my ($self, $where) = @_;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $ga = $opdba->get_GeneAdaptor();
    if (!$ga) {
        $self->_error("Could not get GeneAdaptor");
    }

    my $count = $ga->fetch_gene_count($where);
    if (!$count) {
        $self->_error("Could not fetch oPOSSUM gene count");
    }

    #$self->state->db_info($db_info);

    return $count;
}

sub fetch_db_info
{
    my $self = shift;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $dbia = $opdba->get_DBInfoAdaptor();
    if (!$dbia) {
        $self->_error("Could not get DBInfoAdaptor");
    }

    my $db_info = $dbia->fetch_db_info();
    if (!$db_info) {
        $self->_error("Could not fetch DB info");
    }

    #$self->state->db_info($db_info);

    return $db_info;
}

sub fetch_external_gene_id_types
{
    my $self = shift;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $xgita = $opdba->get_ExternalGeneIDTypeAdaptor();
    if (!$xgita) {
        $self->_error("Could not get ExternalGeneIDTypeAdaptor");
    }

    my $xgits = $xgita->fetch_where();
    if (!$xgits) {
        $self->_error("Could not fetch external gene ID types");
    }

    #$self->state->external_gene_id_types($xgits);

    return $xgits;
}

sub fetch_conservation_levels
{
    my $self = shift;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $cla = $opdba->get_ConservationLevelAdaptor();
    if (!$cla) {
        $self->_error("Could not get ConservationLevelAdaptor");
    }

    my $clh = $cla->fetch_conservation_level_hash();
    if (!$clh) {
        $self->_error("Could not fetch conservation levels");
    }

    #$self->state->conservation_level_hash($clh);

    return $clh;
}

sub fetch_cr_gc_content
{
    my ($self, $gids, $clevel, $upstream_bp, $downstream_bp, $biotype) = @_;
    
    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }
    
    if (!$gids or scalar(@$gids) == 0) {
        # all gene ids
        my $ga = $opdba->get_GeneAdaptor();

        my $where = _biotype_where_clause($biotype);

        $gids = $ga->fetch_gene_ids($where);
    }
    
    my $cra = $opdba->get_ConservedRegionAdaptor();
    if (!$cra) {
        $self->error("Could not get ConservedRegionAdaptor");
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

sub fetch_threshold_levels
{
    my $self = shift;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $tla = $opdba->get_ThresholdLevelAdaptor();
    if (!$tla) {
        $self->_error("Could not get ThresholdLevelAdaptor");
    }

    my $tlh = $tla->fetch_threshold_level_hash();
    if (!$tlh) {
        $self->_error("Could not fetch threshold levels");
    }

    #$self->state->threshold_level_hash($tlh);

    return $tlh;
}

sub fetch_search_region_levels
{
    my $self = shift;

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $srla = $opdba->get_SearchRegionLevelAdaptor();
    if (!$srla) {
        $self->_error("Could not get SearchRegionLevelAdaptor");
    }

    my $srlh = $srla->fetch_search_region_level_hash();
    if (!$srlh) {
        $self->_error("Could not fetch search region levels");
    }

    #$self->state->search_region_level_hash($srlh);

    return $srlh;
}


sub parse_gene_id_text
{
    my ($self, $text) = @_;

    #
    # Strip anything out that is NOT a part of the gene ID or a valid
    # separator
    #
    $text =~ s/[^\w\.\/_\-,;:\s\n]+//g;

    #
    # Strip out leading and trailing separators
    #
    $text =~ s/^[,;:\s\n]+//g;
    $text =~ s/[,;:\s\n]+$//g;

    #print LOG "processed gene text:\n"
    #   . Data::Dumper::Dumper($gene_text) . "\n";

    my @raw_gene_list = split /[,;:\n\s]+/, $text;
    #print LOG "raw gene list:\n"
    #   . Data::Dumper::Dumper(@raw_gene_list) . "\n";

    my %gene_included;
    my @unique_gene_list;
    if (@raw_gene_list) {
        foreach my $gene (@raw_gene_list) {
            unless ($gene_included{$gene}) {
                push @unique_gene_list, $gene;
                $gene_included{$gene} = 1;
            }
        }
    }

    return @unique_gene_list ? \@unique_gene_list : undef;
}

sub get_t_gene_file
{
    my ($self, $id_input_method, $tempdir) = @_;

    my $q = $self->query;

    my $filename;
    if ($id_input_method eq "paste") {
        my $sl = $q->param('t_gene_text');
        if (!$sl) {
            $self->_error("Target gene ID input not specified");
            return;
        }

        $filename = "$tempdir/t_gene_ids.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target gene IDs file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($id_input_method eq "upload") {
        my $file = $q->param('t_gene_file');
        my $fh   = $q->upload('t_gene_file');

        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("File $file is empty\n");
            return;
        }

        $filename = "$tempdir/t_gene_ids.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create target gene IDs file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown target gene ID input method");
        return;
    }

    return $filename;
}

sub get_bg_gene_file
{
    my ($self, $id_input_method, $tempdir) = @_;

    my $q = $self->query;

    my $filename;
    if ($id_input_method eq "paste") {
        my $sl = $q->param('bg_gene_text');
        if (!$sl) {
            $self->_error("Background gene ID input not specified");
            return;
        }

        $filename = "$tempdir/bg_gene_ids.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create background gene IDs file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($id_input_method eq "upload") {
        my $file = $q->param('bg_gene_file');
        my $fh   = $q->upload('bg_gene_file');

        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("File $file is empty\n");
            return;
        }

        $filename = "$tempdir/bg_gene_ids.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create background gene IDs file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown background gene ID input method");
        return;
    }

    return $filename;
}

1;
