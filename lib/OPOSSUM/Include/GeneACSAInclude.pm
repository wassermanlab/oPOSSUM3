use strict;

use OPOSSUM::Include::GeneInclude;
use OPOSSUM::Include::ACSAInclude;

use constant BG_COLOR_CLASS => 'bgc_gene_acsa';

#
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
	my ($filename, $results, $tf_set, %job_args) = @_;

    return unless $results && $results->[0];

	my $logger = $job_args{-logger};

    open(FH, ">$filename")
        || fatal("Could not create analysis results file $filename", %job_args);

    $logger->info("Writing analysis results to $filename\n");

    printf FH "TF Name\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tGC Content\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tf = $tf_set->get_tf($result->id());

        my $gc_content = sprintf("%.3f", $tf->tag('gc_content'));

        printf FH 
            "%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",
            $tf->name(),
            $tf->ID(),
            $tf->class() || 'NA',
            $tf->tag('family') || 'NA',
            $tf->tag('tax_group') || 'NA',
            $tf->to_ICM->total_ic(),
            $gc_content,
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_tfbs_hits() || 0,
            defined $result->t_tfbs_rate()
                ? sprintf("%.3g", $result->t_tfbs_rate()) : 'NA',
            $result->bg_tfbs_hits() || 0,
            defined $result->bg_tfbs_rate()
                ? sprintf("%.3g", $result->bg_tfbs_rate()) : 'NA',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'NA',
            defined $result->fisher_p_value()
                ? sprintf("%.3f", $result->fisher_p_value()) : 'NA';
    }
    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/gene. Create a
# text file for each TF.
#
=head3
sub write_tfbs_details
{
	my ($results, $t_counts, $t_operon_first_gids, $has_operon,
		$anchor_matrix, $max_site_dist, $ga, $ctfbsa,
		$conservation_level, $threshold, $upstream_bp, $downstream_bp,
		$abs_results_dir, $rel_results_dir, $web, %job_args) = @_;

	my $species = $job_args{-species};

	# if results are truncated, show the relevant tf ids only
    my $tf_ids;
    my $t_gids = $t_counts->gene_ids;
	foreach my $result (@$results) {
		push @$tf_ids, $result->id;
	}

	my $anchor_tf_id = $anchor_matrix->ID;
	
    foreach my $tf_id (@$tf_ids) 
	{
        my $tf = $tf_set->get_tf($tf_id);

        my $tf_name = $tf->name();

        my $text_filename = "$abs_results_dir/$tf_id.txt";
        my $html_filename = "$abs_results_dir/$tf_id.html";

        # XXX
        # Fetch TFBSs for this TF and all genes and pass to routines
        # below.
        #
        # NOTE: for operon genes only fetch for first gene in operon.
        #
        my @tf_genes;
        my %gid_sitepairs;
        my %gid_inc;
        my %operon_genes;
        foreach my $t_gid (@$t_gids) {
            next unless $t_counts->gene_tfbs_count($t_gid, $tf_id);

            my $gid;
            if ($has_operon) {
                my $fgid = $t_operon_first_gids->{$t_gid};
                if ($fgid) {
                    $gid = $fgid;
                } else {
                    $gid = $t_gid;
                }

                my $gene = $ga->fetch_by_id($t_gid);
                push @{$operon_genes{$gid}}, $gene;
            } else {
                $gid = $t_gid;
            }

            #
            # For operon genes, this gene may have already been included via
            # another gene in the same operon
            #
            next if $gid_inc{$gid};

            $gid_inc{$gid} = 1;

            my $sitepairs = $ctfbsa->fetch_anchored_by_tf_id(
                -tf_id              => $tf_id,
                -anchor_tf_id       => $anchor_tf_id,
                -gene_id            => $gid,
                -distance           => $max_site_dist,
                -conservation_level => $conservation_level,
                -threshold          => $threshold,
                -upstream_bp        => $upstream_bp,
                -downstream_bp      => $downstream_bp
            );

            if ($sitepairs) {
                #
                # This should always be true, since we previously checked
                # the count was not zero.
                #
                my $gene = $ga->fetch_by_id($gid);
                $gene->fetch_promoters();

                push @tf_genes, $gene;

                $gid_sitepairs{$gid} = $sitepairs;
            }
        }

        write_tfbs_details_text(
            $text_filename, $tf, \@tf_genes, $anchor_matrix,
			\%gid_sitepairs, \%operon_genes,
			%job_args
        );

        write_tfbs_details_html(
            $html_filename, $tf, $anchor_matrix,
			$rel_results_dir, $species,
			$tf, \@tf_genes, \%gid_sitepairs, \%operon_genes, 
			%job_args
        ) if $web;
    }
}
=cut

#
# Write the details of the putative TFBSs for each TF/gene. Create a
# text file for each TF.
#
sub write_tfbs_details_text
{
    my ($filename, $tf, $anchor_matrix,
	   	$genes, $gid_sitepairs, 
		$gid_gene_ids, $gene_id_type,
		$has_operon, $operon_genes, %job_args) = @_;

	my $logger = $job_args{-logger};

    my $tf_name = $tf->name();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details text file $filename"
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    #
    # It's not really necessary to check this for gene based oPOSSUM as
    # all profiles come from JASPAR and are retrieved as PFMs.
    #
    my $tf_total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $tf_total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $tf_total_ic = 'NA';
    }

    my $anchor_total_ic;
    if ($anchor_matrix->isa("TFBS::Matrix::PFM")) {
        $anchor_total_ic = sprintf("%.3f", $anchor_matrix->to_ICM->total_ic());
    } else {
        $anchor_total_ic = 'NA';
    }

    my $text = sprintf("%s\n\n", $anchor_matrix->name());

    $text .= sprintf("JASPAR ID:\t%s\n", $anchor_matrix->ID());
    $text .= sprintf("Class:    \t%s\n", $anchor_matrix->class() || 'NA');
    $text .= sprintf("Family:   \t%s\n", $anchor_matrix->tag('family') || 'NA');
    $text .= sprintf("Tax group:\t%s\n",
        $anchor_matrix->tag('tax_group') || 'NA');
    $text .= sprintf("IC:       \t%s\n", $tf_total_ic);

    $text .= "\n\n";

    $text .= sprintf("%s\n\n", $tf->name());

    $text .= sprintf("JASPAR ID:\t%s\n", $tf->ID());
    $text .= sprintf("Class:    \t%s\n", $tf->class() || 'NA');
    $text .= sprintf("Family:   \t%s\n", $tf->tag('family') || 'NA');
    $text .= sprintf("Tax group:\t%s\n", $tf->tag('tax_group') || 'NA');
    $text .= sprintf("IC:       \t%s\n", $tf_total_ic);

    $text .= sprintf("\n\n%s - %s Conserved Binding Site Combinations\n\n",
        $anchor_matrix->name(), $tf->name());

    unless (defined $gene_id_type
            && $gene_id_type == DFLT_GENE_ID_TYPE)
    {
        $text .= "Gene ID(s)\t";
    }

    $text .= "Ensembl ID";

    if ($has_operon) {
        $text .= "\tOperon ID";
    }

    $text .= qq{\tChr\tStart\tEnd\tStrand\tNearest TSS\tAnchoring TF\tStart\tEnd\tRel. Start\tRel. End\tStrand\tScore\t%Score\tSequence\tAnchored TF\tStart\tEnd\tRel. Start\tRel. End\tStrand\tScore\t%Score\tSequence\tDistance\n};

    foreach my $gene (@$genes) {
        my $gid = $gene->id();

        my $sitepairs = $gid_sitepairs->{$gid};

        my $ensembl_id  = $gene->ensembl_id();
        my $gene_start  = $gene->start();
        my $gene_end    = $gene->end();
        my $strand      = $gene->strand();
        my $chr         = $gene->chr();
        my $tss         = $gene->tss();
        my $promoters   = $gene->promoters();

        my $prom_start;
        my $prom_end;
        if ($strand == 1) {
            $prom_start = $tss;
            $prom_end   = $gene_end;
        } else {
            $prom_start = $gene_start;
            $prom_end   = $tss;
        }

        if ($has_operon && defined $operon_genes->{$gid}) {
            unless (defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE
            ) {
                my $count = 0;
                foreach my $in_gene (@{$operon_genes->{$gid}}) {
                    $text .= ',' if $count > 0;
                    $text .= join(',', @{$gid_gene_ids->{$in_gene->id()}});
                    $count++;
                }
                $text .= "\t";
            }
            
            my @ensembl_ids;
            foreach my $in_gene (@{$operon_genes->{$gid}}) {
                #push @ensembl_ids, $t_gid_ens_ids->{$in_gene->id()};
                push @ensembl_ids, $in_gene->ensembl_id();
            }
            $text .= join(',', @ensembl_ids);
        } else {
            unless (defined $gene_id_type
                    && $gene_id_type == DFLT_GENE_ID_TYPE)
            {
                $text .= join(',', @{$gid_gene_ids->{$gid}}) . "\t";
            }

            $ensembl_id  = $gene->ensembl_id();
            $strand      = $gene->strand();
            $chr         = $gene->chr();
            $tss         = $gene->tss(); # this tss is not used.
                                         # overwritten later.
            $text .= "$ensembl_id";
        }
        
        if ($has_operon) {
            my $symbol = "NA";
            if (defined $gene->operon()) {
                $symbol = $gene->operon->symbol();
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
        foreach my $sitepair (@$sitepairs) {
            my $tfbs        = $sitepair->{tf_site};
            my $anchor      = $sitepair->{anchor_site};
            my $distance    = $sitepair->{distance};

            my $site_start          = $gene_start + $tfbs->start() - 1;
            my $site_end            = $gene_start + $tfbs->end() - 1;
            my $site_seq            = $tfbs->seq();
            my $site_strand         = $tfbs->strand();
            my $site_score          = $tfbs->score();
            my $site_rel_score      = $tfbs->rel_score();

            my $anchor_start        = $gene_start + $anchor->start() - 1;
            my $anchor_end          = $gene_start + $anchor->end() - 1;
            my $anchor_seq          = $anchor->seq();
            my $anchor_strand       = $anchor->strand();
            my $anchor_score        = $anchor->score();
            my $anchor_rel_score    = $anchor->rel_score();

            my $anchor_closest_tss;
            my $site_closest_tss;
            my $min_anchor_tss_dist = 999999;
            my $min_site_tss_dist = 999999;
            foreach my $promoter (@$promoters) {
                my $tss = $promoter->tss();

                my $anchor_start_tss_dist = abs($anchor_start - $tss);
                my $anchor_end_tss_dist   = abs($anchor_end - $tss);

                if ($anchor_start_tss_dist < $min_anchor_tss_dist) {
                    $min_anchor_tss_dist = $anchor_start_tss_dist;
                    $anchor_closest_tss = $tss;
                }

                if ($anchor_end_tss_dist < $min_anchor_tss_dist) {
                    $min_anchor_tss_dist = $anchor_end_tss_dist;
                    $anchor_closest_tss = $tss;
                }

            }
            
            my ($site_rel_start, $site_rel_end);
            my ($anchor_rel_start, $anchor_rel_end);
            if ($strand == 1) {
                $anchor_rel_start = $anchor_start - $anchor_closest_tss;
                if ($anchor_start >= $anchor_closest_tss) {
                    $anchor_rel_start++;
                }

                $anchor_rel_end = $anchor_end - $anchor_closest_tss;
                if ($anchor_end >= $anchor_closest_tss) {
                    $anchor_rel_end++;
                }

                $site_rel_start = $site_start - $anchor_closest_tss;
                if ($site_start >= $anchor_closest_tss) {
                    $site_rel_start++;
                }

                $site_rel_end = $site_end - $anchor_closest_tss;
                if ($site_end >= $anchor_closest_tss) {
                    $site_rel_end++;
                }
            } else {
                $anchor_rel_start = $anchor_closest_tss - $anchor_start;
                if ($anchor_start <= $anchor_closest_tss) {
                    $anchor_rel_start++;
                }

                $anchor_rel_end = $anchor_closest_tss - $anchor_end;
                if ($anchor_end <= $anchor_closest_tss) {
                    $anchor_rel_end++;
                }

                # swap coords so start is more upstream than end
                ($anchor_rel_start, $anchor_rel_end)
                    = ($anchor_rel_end, $anchor_rel_start);

                $site_rel_start = $anchor_closest_tss - $site_start;
                if ($site_start <= $anchor_closest_tss) {
                    $site_rel_start++;
                }

                $site_rel_end = $anchor_closest_tss - $site_end;
                if ($site_end <= $anchor_closest_tss) {
                    $site_rel_end++;
                }

                # swap coords so start is more upstream than end
                ($site_rel_start, $site_rel_end)
                    = ($site_rel_end, $site_rel_start);
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

            $text .= sprintf("\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\t%d\n",
                $anchor_closest_tss,
                $anchor_matrix->name(),
                $anchor_start,
                $anchor_end,
                $anchor_rel_start,
                $anchor_rel_end,
                $anchor_strand == 1 ? '+' : '-',
                $anchor_score,
                $anchor_rel_score * 100,
                $anchor_seq,
                $tf->name(),
                $site_start,
                $site_end,
                $site_rel_start,
                $site_rel_end,
                $site_strand == 1 ? '+' : '-',
                $site_score,
                $site_rel_score * 100,
                $site_seq,
                $distance
            );

            $first = 0;
        } # end foreach tfbs
    } # end foreach gid

    unless (open(FH, ">$filename")) {
        $logger->error("Unable to create TFBS details file $filename - $!");
        return;
    }

    print FH $text;

    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/gene. Create an
# HTML file for each TF.
#
sub write_tfbs_details_html
{
    my ($filename, $rel_results_dir, $species, $tf_db, $tf, $anchor_matrix, 
		$genes, $gid_sitepairs,
		$gid_gene_ids, $gene_id_type,
		$has_operon, $operon_genes, %job_args) = @_;

	my $logger = $job_args{-logger};
	my $heading = $job_args{-heading};

    my $tf_name = $tf->name();
    my $tf_id   = $tf->ID();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details html file $filename"
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

	#my $heading = sprintf "%s Anchored Combination Site Analysis",
	#    ucfirst $species;

    my $title = "oPOSSUM $heading";

    my $section = sprintf(
        "%s - %s Conserved Binding Site Combinations",
        $anchor_matrix->name(), $tf_name
    );

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        bg_color_class      => BG_COLOR_CLASS,
        title               => $title,
        heading             => $heading,
        section             => $section,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        low_matrix_ic       => LOW_MATRIX_IC,
        high_matrix_ic      => HIGH_MATRIX_IC,
        low_matrix_gc       => LOW_MATRIX_GC,
        high_matrix_gc      => HIGH_MATRIX_GC,
        #low_seq_gc          => LOW_SEQ_GC,
        #high_seq_gc         => HIGH_SEQ_GC,

        jaspar_url          => JASPAR_URL,

        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'NA'
                               },

        tf_db               => $tf_db,
        anchor_matrix       => $anchor_matrix,
        tf                  => $tf,
        gene_id_type        => $gene_id_type,
        dflt_gene_id_type   => DFLT_GENE_ID_TYPE,
        genes               => $genes,
        gid_gene_ids        => $gid_gene_ids,
        gid_sitepairs       => $gid_sitepairs,
        has_operon          => $has_operon,
        operon_genes        => $operon_genes,
        rel_results_dir     => $rel_results_dir,
        tfbs_details_file   => "$tf_id.txt",
        var_template        => "tfbs_details_gene_acsa.html"
    };

    my $output = process_template('master.html', $vars);

    print FH $output;

    close(FH);
}

1;

