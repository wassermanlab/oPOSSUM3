#
# This module should be included in all the oPossumGene*Web.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to all the oPossum gene-based variants.
#

use oPossumWebOpt;
use oPossumGeneWebOpt;

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
    
    my $where = _biotype_where_clause($biotype);
    
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

sub fetch_analysis_counts
{
    my ($self, $analysis_type, %args) = @_;

    #printf STDERR "fetch_analysis_counts ($analysis_type) args:\n"
    #    . Data::Dumper::Dumper(%args) . "\n\n";

    my $opdba = $self->opdba();
    if (!$opdba) {
        $opdba = $self->opossum_db_connect();
    }

    my $aca = $opdba->get_AnalysisCountsAdaptor();
    if (!$aca) {
        $self->_error("Could not get AnalysisCountsAdaptor");
    }

    my $counts;
    if ($analysis_type eq 'default') {
        #printf STDERR "\nfetch_analysis_counts: fetching pre-computed counts\n";
        $counts = $aca->fetch_counts(%args);
    } elsif ($analysis_type eq 'custom') {
        #printf STDERR "\nfetch_analysis_counts: fetching custom counts\n";
        $counts = $aca->fetch_custom_counts(%args);
    } else {
        $self->_error("fetch_analysis_counts: unknown analysis type");
    }

    if (!$counts) {
        $self->_error("Could not fetch analysis counts");
    }

    #$self->state->search_region_level_hash($srlh);

    return $counts;
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

#
# For the provided list of Ensembl / external gene IDs, fetch the associated
# oPOSSUM gene IDs from the DB. Also keep track of which provided gene IDs
# mapped (included) or did not map (missing) to oPOSSUM gene IDs.
#
# NOTE: There can be a 1-to-many mapping of oPOSSUM gene IDs to external gene
# IDs. Make sure all of these are retained in the included external gene IDs
# list.
#
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
    my ($self, $id_type, $ids, $has_operon) = @_;

    my $state = $self->state();

    my $opdba = $self->opdba();

    my $ga = $opdba->get_GeneAdaptor();
    unless ($ga) {
        $self->_error("Error getting GeneAdaptor");
        return;
    }

    my @mapped_ids;
    my @unmapped_ids;
    my %gene_id_ext_ids;

    #
    # This will be set to a mapping of operon gene IDs to the first gene ID
    # in the operon.
    #
    my %operon_gene_id_to_first_gene_id;

    #
    # This will be set to a list of gene IDs for both genes which do not
    # belong to an operon and genes which are the first gene in an operon 
    # that one or more of the input genes belongs to.
    #
    my @operon_unique_gene_ids;

    #
    # These are not ever used
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
            -gene_id_ensembl_ids    => \%gene_id_ext_ids
            #-ensembl_id_gene_ids    => \%ext_id_gene_ids
        );
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
            #-ext_id_gene_ids        => \%ext_id_gene_ids,
            #-gene_id_ensembl_ids    => \%gene_id_ensembl_ids
        );
    } 

    if ($has_operon) {
        my $oa = $opdba->get_OperonAdaptor();
        unless ($oa) {
            $self->_error("Error getting OperonAdaptor");
            return;
        }

        my %fgene_id_inc;
        foreach my $gene_id (@$gene_ids) {
            # first check to see if this gene belongs to an operon
            # it is possible to have multiple genes belonging to the same operon
            # to be included in the genes to be fetched.
            # thus need to keep an array of gene ids
            # map each gene id in the operon to the first gene id
            my $operon = $oa->fetch_by_gene_id($gene_id);
            if ($operon) {
                my $fgene = $operon->fetch_first_gene();

                $operon_gene_id_to_first_gene_id{$gene_id} = $fgene->id;

                unless ($fgene_id_inc{$fgene->id}) {
                    push @operon_unique_gene_ids, $fgene->id;
                    $fgene_id_inc{$fgene->id} = 1;
                }
            } else {
                push @operon_unique_gene_ids, $gene_id;
            }
        }
    } 

    #print STDERR "fetch_opossum_gene_ids: # operon genes = " . scalar(keys %operon_gids) . "\n"; 

    return (
        $gene_ids,
        @mapped_ids             ? \@mapped_ids : undef,
        @unmapped_ids           ? \@unmapped_ids : undef,
        #%gene_id_ensembl_ids    ? \%gene_id_ensembl_ids : undef,
        %gene_id_ext_ids        ? \%gene_id_ext_ids : undef,
        #%ext_id_gene_ids        ? \%ext_id_gene_ids : undef,
        #%operon_first_gene_to_gene_ids
        #                        ? \%operon_first_gene_to_gene_ids : undef
        %operon_gene_id_to_first_gene_id
                                ? \%operon_gene_id_to_first_gene_id : undef,
        @operon_unique_gene_ids ? \@operon_unique_gene_ids : undef
    );
}

sub fetch_random_opossum_gene_ids
{
    my ($self, $num_genes, $has_operon, $biotype) = @_;

    return if !$num_genes;

    my $opdba = $self->opdba();

    my $ga = $opdba->get_GeneAdaptor();
    unless ($ga) {
        $self->_error("Error getting GeneAdaptor");
        return;
    }

    my %operon_gene_id_to_first_gene_id;
    my @operon_unique_gene_ids;
    
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

    if ($has_operon) {
        my $oa = $opdba->get_OperonAdaptor();
        unless ($oa) {
            $self->_error("Error getting OperonAdaptor");
            return;
        }

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
                    $operon_gene_id_to_first_gene_id{$gid} = $fgene->id;
                #}

                unless ($fgene_id_inc{$fgene->id}) {
                    push @operon_unique_gene_ids, $fgene->id;
                    $fgene_id_inc{$fgene->id} = 1;
                }
            } else {
                push @operon_unique_gene_ids, $gid;
            }
        }
    }
    
    return (
        $rand_gids,
        %operon_gene_id_to_first_gene_id
                                ? \%operon_gene_id_to_first_gene_id : undef,
        @operon_unique_gene_ids ? \@operon_unique_gene_ids : undef
    );
} 

sub write_results
{
    my ($self, $filename, $results, $tf_set) = @_;

    return unless $results && $results->[0];

    my $text = "TF\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tf = $tf_set->get_tf($result->id());

        my $total_ic;
        if ($tf->isa("TFBS::Matrix::PFM")) {
            $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
        } else {
            $total_ic = 'N/A';
        }

        $text .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n",
            $tf->name,
            $tf->ID,
            $tf->class() || 'N/A',
            $tf->tag('family') || 'N/A',
            $tf->tag('tax_group') || 'N/A',
            $total_ic,
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_tfbs_hits() || 0,
            defined $result->t_tfbs_rate()
                ? sprintf("%.3f", $result->t_tfbs_rate()) : 'N/A',
            $result->bg_tfbs_hits() || 0,
            defined $result->bg_tfbs_rate()
                ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'N/A',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'N/A',
            defined $result->fisher_p_value()
                ? sprintf("%.3g", $result->fisher_p_value()) : 'N/A';
    }

    unless (open(FH, ">$filename")) {
        $self->_error("Unable to create results text file $filename - $!");
        return;
    }

    print FH $text;

    close(FH);
}

sub write_tfbs_details
{
    my ($self, $filename, $tf, $genes, $gid_tfbss, $gene_id_type,
        $gid_gene_ids, $has_operon, $operon_genes) = @_;

    my $total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $total_ic = 'N/A';
    }

    my $text = sprintf("%s\n\n", $tf->name());

    $text .= sprintf("JASPAR ID:\t%s\n", $tf->ID());
    $text .= sprintf("Class:\t%s\n", $tf->class() || 'N/A');
    $text .= sprintf("Family:\t%s\n", $tf->tag('family') || 'N/A');
    $text .= sprintf("Tax group:\t%s\n", $tf->tag('tax_group') || 'N/A');
    $text .= sprintf("Information content:\t%s\n", $total_ic);

    $text .= sprintf("\n\nConserved %s Binding Sites:\n\n", $tf->name());

    unless (defined $gene_id_type && $gene_id_type == DFLT_GENE_ID_TYPE) {
        $text .= "Gene ID(s)\t";
    }

    $text .= "Ensembl ID(s)";

    if ($has_operon) {
        $text .= "\tOperon ID";
    }

    $text .= qq{\tChr\tStart\tEnd\tStrand\tNearest TSS\tTFBS Start\tTFBS End\tTFBS Rel. Start\tTFBS Rel. End\tTFBS Strand\tAbs. Score\tRel. Score\tTFBS Sequence\n};

    #
    # Notes for operon genes
    # $genes -> only contain first operon genes. However, not all first genes
    # may have been included in the input gene set. These can be checked with
    # $operon_genes -> contain all operon genes in the input set
    # So, whenever a gene in $genes maps to operon_genes, must get the
    # corresponding gene array from operon_genes and list them
    #
    foreach my $gene (@$genes) {
        my $gid = $gene->id;

        my $ensembl_id  = $gene->ensembl_id();
        my $chr         = $gene->chr();
        my $gene_start  = $gene->start();
        my $gene_end    = $gene->end();
        my $strand      = $gene->strand();
        my $tss         = $gene->tss();
        my $promoters   = $gene->fetch_promoters();

        my $operon;
        if ($has_operon) {
            $operon = $gene->operon();
        }

        my $prom_start;
        my $prom_end;
        if ($operon) {
            my $operon_start    = $operon->start();
            my $operon_end      = $operon->end();

            if ($strand == 1) {
                $prom_start = $tss;
                $prom_end   = $operon_end;
            } else {
                $prom_start = $operon_start;
                $prom_end   = $tss;
            }
        } else {
            if ($strand == 1) {
                $prom_start = $tss;
                $prom_end   = $gene_end;
            } else {
                $prom_start = $gene_start;
                $prom_end   = $tss;
            }
        }

        if ($has_operon && defined $operon_genes->{$gid}) {
            unless (defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE)
            {
                my $count = 0;
                foreach my $in_gene (@{$operon_genes->{$gid}}){
                    $text .= ',' if $count > 0;
                    $text .= join(',', @{$gid_gene_ids->{$in_gene->id}});
                    $count++;
                }
                $text .= "\t";
            }
            
            my @ensembl_ids;
            foreach my $in_gene (@{$operon_genes->{$gid}}){
                push @ensembl_ids, $in_gene->ensembl_id();
            }
            $text .= join(',', @ensembl_ids) . "\t";
        } else {
            unless (defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE)
            {
                $text .= join(',', @{$gid_gene_ids->{$gid}});
            }
            $text .= "\t";

            $text .= $gene->ensembl_id() . "\t";
        }
            
        if ($has_operon) {
            my $symbol = "-";
            if (defined $operon) {
                $symbol = $operon->symbol();
            }
            $text .= sprintf("\t%s", $symbol);
        }

        $text .= sprintf("\t%s\t%d\t%d\t%s",
            $chr,
            $prom_start,
            $prom_end,
            $strand == 1 ? '+' : '-'
        );

        # gene_start/gene_end should refer to the first gene here
        # what about tss? which one should I use? show both?
        # actually, no. closest tss refers to the actual promoter tss
        # since the first gene tss is used, should be kept.
        
        my $first = 1;
        foreach my $tfbs (@{$gid_tfbss->{$gid}}) {
            my $site_start      = $gene_start + $tfbs->start() - 1;
            my $site_end        = $gene_start + $tfbs->end() - 1;
            my $site_seq        = $tfbs->seq();
            my $site_strand     = $tfbs->strand();
            my $site_score      = $tfbs->score();
            my $site_rel_score  = $tfbs->rel_score();

            my $closest_tss;
            my $min_tss_dist = 999999;
            foreach my $promoter (@$promoters) {
                my $tss = $promoter->tss();

                my $start_tss_dist = abs($site_start - $tss);
                my $end_tss_dist   = abs($site_end - $tss);

                if ($start_tss_dist < $min_tss_dist) {
                    $min_tss_dist = $start_tss_dist;
                    $closest_tss = $tss;
                }

                if ($end_tss_dist < $min_tss_dist) {
                    $min_tss_dist = $end_tss_dist;
                    $closest_tss = $tss;
                }
            }
            
            my ($rel_start, $rel_end);
            if ($strand == 1) {
                $rel_start = $site_start - $closest_tss;
                if ($site_start >= $closest_tss) {
                    $rel_start++;
                }

                $rel_end = $site_end - $closest_tss;
                if ($site_end >= $closest_tss) {
                    $rel_end++;
                }
            } else {
                $rel_start = $closest_tss - $site_start;
                if ($site_start <= $closest_tss) {
                    $rel_start++;
                }

                $rel_end = $closest_tss - $site_end;
                if ($site_end <= $closest_tss) {
                    $rel_end++;
                }

                # swap coords so start is more upstream than end
                ($rel_start, $rel_end) = ($rel_end, $rel_start);
            }

            unless ($first) {
                $text .= "\t\t\t\t\t";
                if ($has_operon) {
                    $text .= "\t";
                }

                unless (   defined $gene_id_type
                        && $gene_id_type == DFLT_GENE_ID_TYPE
                ) {
                    $text .= "\t";
                }
            }

            $text .= sprintf("\t%d\t%d\t%d\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\n",
                $closest_tss,
                $site_start,
                $site_end,
                $rel_start,
                $rel_end,
                $site_strand == 1 ? '+' : '-',
                $site_score,
                $site_rel_score * 100,
                $site_seq
            );

            $first = 0;
        } # end foreach tfbs
    } # end foreach gene

    unless (open(FH, ">$filename")) {
        $self->_error("Unable to create TFBS details file $filename - $!");
        return;
    }

    print FH $text;
    close(FH);
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
