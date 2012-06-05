=head1 NAME

 OPOSSUM::Include::GeneInclude.pm

=head1 SYNOPSIS

=head1 DESCRIPTION

  Contains all options and routines that are common to all the gene-based
  analysis scripts.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

#use oPossumWebOpt;
use OPOSSUM::Opt::GeneOpt;
use OPOSSUM::Include::BaseInclude;

#use Data::Dumper;    # for debugging only

use Bio::SeqIO;

#use TFBS::DB::JASPAR5;


sub read_gene_ids_from_file
{
	my ($file, %job_args) = @_;

	return read_ids_from_file($file, %job_args);
}


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

sub fetch_all_gene_ids
{
    # AK: added biotype parameter
    my ($ga, $biotype) = @_;
    
    # let's skip over 'pseudogenes';
    #my $where = _biotype_where_clause($biotype);
    my $where = "biotype != 'pseudogene'";
    
    my $gids = $ga->fetch_gene_ids($where);

    return $gids;
}

sub fetch_gene_count
{
    # AK: added where parameter
    my ($ga, $where) = @_;

    my $count = $ga->fetch_gene_count($where);
    #if (!$count) {
    #    $self->_error("Could not fetch oPOSSUM gene count");
    #}

    return $count;
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


sub fetch_gene_data_for_opossum_analysis
{
    my ($ga, $oa, $cla, $has_operon, $biotype,
        $t_gene_id_type, $t_gene_ids,
        $bg_gene_id_type, $bg_gene_ids, $bg_num_rand_genes,
        %fatal_args
        ) = @_;
    
    my $job_id = $fatal_args{-job_id};
    my $heading = $fatal_args{-heading};
    my $email = $fatal_args{-email};
    my $logger = $fatal_args{-logger};

    $logger->info("Fetching target oPOSSUM gene IDs");
    my (
        $t_gids,
        $t_included_gene_ids,
        $t_missing_gene_ids,
        #$t_gid_ens_ids,
        $t_gid_gene_ids,
        #t_gene_id_gids,
        $t_operon_first_gids,
        $t_operon_unique_gids
    ) = fetch_opossum_gene_ids(
        $ga, $oa, $t_gene_id_type, $t_gene_ids, $has_operon
    );
    
    if (!$t_gids || !$t_gids->[0]) {
        fatal(
              "Error fetching target oPOSSUM gene IDs. Please make sure"
            . " you specified the correct gene ID type for the IDs entered.",
            %fatal_args
        );
    }
    
    my $bg_gids;
    my $bg_included_gene_ids;
    my $bg_missing_gene_ids;
    my $bg_gid_gene_ids;
    my $bg_operon_first_gids;
    my $bg_operon_unique_gids;
    if ($bg_gene_ids) {
        $logger->info("Fetching background oPOSSUM gene IDs");
        ($bg_gids,
         $bg_included_gene_ids,
         $bg_missing_gene_ids,
         $bg_gid_gene_ids,
         $bg_operon_first_gids,
         $bg_operon_unique_gids
        ) = fetch_opossum_gene_ids(
            $ga, $oa, $bg_gene_id_type, $bg_gene_ids, $has_operon
        );
    
        if (!$bg_gids || !$bg_gids->[0]) {
            fatal(
                  "Error fetching background oPOSSUM gene IDs. Please make"
                . " sure you specified the correct gene ID type for the gene"
                . " IDs entered.",
                %fatal_args
            );
        }
    } elsif ($bg_num_rand_genes) {
        $logger->info("Fetching random background oPOSSUM gene IDs");
        ($bg_gids,
         $bg_operon_first_gids,
         $bg_operon_unique_gids
        ) = fetch_random_opossum_gene_ids(
            $ga, $oa, $bg_num_rand_genes, $has_operon, $biotype
        );
    
        if (!$bg_gids || !$bg_gids->[0]) {
            fatal(
                "Error fetching random background oPOSSUM gene IDs", %fatal_args
            );
        }
    } else {
        $logger->info("Fetching all background oPOSSUM gene IDs");
        $bg_gids = $ga->fetch_gene_ids();
    }
    
    return (
        $t_gids,
        $t_included_gene_ids,
        $t_missing_gene_ids,
        $t_gid_gene_ids,
        $t_operon_first_gids,
        $t_operon_unique_gids,
        $bg_gids,
        $bg_included_gene_ids,
        $bg_missing_gene_ids,
        $bg_gid_gene_ids,
        $bg_operon_first_gids,
        $bg_operon_unique_gids
    );
}

#
# Set default or custom parameter settings
# 
# XXX
# Previously this routine only returned threshold level and search region level
# but did not return the actual threshold and upstream/downstream bp. Fixed
# so it does this. Also fixed some other bugs.
# DJA 2012/03/05
# XXX
#
sub fetch_analysis_parameters
{
    my ($thla, $srla, $species, $threshold_level, $search_region_level,
        $threshold, $upstream_bp, $downstream_bp) = @_;
    
    my $thlh = $thla->fetch_threshold_level_hash();

    my $analysis_type = 'default';
    if ($threshold_level) {
        # threshold level takes precedence over threshold
        $threshold = $thlh->{$threshold_level}->threshold;
    } elsif (defined $threshold) {
		if ($threshold =~ /(.+)%$/) {
			$threshold = $1 / 100;
		}

        my $found = 0;
        foreach my $level (keys %$thlh) {
            if ($thlh->{$level}->threshold == $threshold) {
                $threshold_level = $level;
                $found = 1;
                last;
            }
        }

        $analysis_type = 'custom' if !$found;
    } else {
        $threshold_level = DFLT_THRESHOLD_LEVEL;
        $threshold = $thlh->{$threshold_level}->threshold;
    }
    
    my $srlh = $srla->fetch_search_region_level_hash();

    if ($search_region_level) {
        # search region level takes precedence
        $upstream_bp   = $srlh->{$search_region_level}->upstream_bp();
        $downstream_bp = $srlh->{$search_region_level}->downstream_bp();
    } else {
        if ($species eq 'yeast') {
            if (defined $upstream_bp) {
                my $found = 0;
                foreach my $level (keys %$srlh) {
                    my $sr_up = $srlh->{$level}->upstream_bp();
                    
                    if ($sr_up == $upstream_bp) {
                        $search_region_level = $level;
                        $found = 1;
                        last;
                    }
                }
                $analysis_type = 'custom' if !$found;
            } else {
                $search_region_level = DFLT_SEARCH_REGION_LEVEL;
                $upstream_bp = $srlh->{$search_region_level}->upstream_bp();
            }
        } else {
            if (defined $upstream_bp && defined $downstream_bp) {
                my $found = 0;
                foreach my $level (keys %$srlh) {
                    my $sr_up   = $srlh->{$level}->upstream_bp();
                    my $sr_down = $srlh->{$level}->downstream_bp();
                    
                    if ($sr_up == $upstream_bp && $sr_down == $downstream_bp) {
                        $search_region_level = $level;
                        $found = 1;
                        last;
                    }
                }
                $analysis_type = 'custom' if !$found;
            } else {
                $search_region_level = DFLT_SEARCH_REGION_LEVEL;
                $upstream_bp   = $srlh->{$search_region_level}->upstream_bp();
                $downstream_bp = $srlh->{$search_region_level}->downstream_bp();
            }
        }
    }

    return ($analysis_type, $threshold_level, $search_region_level, $threshold,
            $upstream_bp, $downstream_bp);
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
    my %operon_gene_id_to_first_gene_id;

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

                $operon_gene_id_to_first_gene_id{$gene_id} = $fgene->id;

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
        %operon_gene_id_to_first_gene_id
                                    ? \%operon_gene_id_to_first_gene_id : undef,
        @operon_unique_gene_ids     ? \@operon_unique_gene_ids : undef
    );
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

    my %operon_gene_id_to_first_gene_id;
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

#
# helper function to be used by write_tfbs_details functions
#
sub _fetch_operon_gene_id
{
    my ($t_gid, $t_operon_first_gids) = @_;
    
    my $gid;

    my $fgid = $t_operon_first_gids->{$t_gid};
    if ($fgid) {
        $gid = $fgid;
    } else {
        $gid = $t_gid;
    }
    
    return $gid;
}

sub _fetch_gene_information_text
{
    my ($gene, $gid_gene_ids, $gene_id_type, $has_operon, $operon_genes) = @_;
    
    my $text;
    
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
        $text .= join(',', @ensembl_ids);
    } else {
        unless (defined $gene_id_type
                && $gene_id_type == DFLT_GENE_ID_TYPE)
        {
            $text .= join(',', @{$gid_gene_ids->{$gid}});
        }
        $text .= "\t";

        $text .= $gene->ensembl_id();
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

    return ($text, $gene_start, $tss, $promoters, $strand);
}

#
# to be used by both gene ssa and gene tca methods
#
sub _fetch_site_information_text
{
    my ($gid, $gid_sites, $gene_start, $tss, $promoters, $strand, 
        $gene_id_type, $has_operon) = @_;
    
    # gene_start/gene_end should refer to the first gene here
    # what about tss? which one should I use? show both?
    # actually, no. closest tss refers to the actual promoter tss
    # since the first gene tss is used, should be kept.
    
    my $text;
    
    my $first = 1;
    foreach my $site (@{$gid_sites->{$gid}}) {
        my $site_start      = $gene_start + $site->start() - 1;
        my $site_end        = $gene_start + $site->end() - 1;
        my $site_seq        = $site->seq();
        my $site_strand     = $site->strand();
        my $site_score      = $site->score();
        my $site_rel_score  = $site->rel_score();

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
            $text .= "\t\t\t\t";
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
    } # end foreach site
    
    return $text;
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
