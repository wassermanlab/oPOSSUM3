=head1 NAME

 OPOSSUM::Include::GeneSSAInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the gene-based SSA
  scripts and web modules.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

#use oPossumWebOpt;
#use oPossumGeneWebOpt;

use OPOSSUM::Include::GeneInclude;

use constant BG_COLOR_CLASS => 'bgc_gene_ssa';

sub write_results_text
{
    my ($filename, $results, $tf_set, %job_args) = @_;

    return unless $results && $results->[0];

    my $text = "TF\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tGC Content\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tf = $tf_set->get_tf($result->id());

        my $total_ic;
        if ($tf->isa("TFBS::Matrix::PFM")) {
            $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
        } else {
            $total_ic = 'NA';
        }

        my $gc_content = sprintf("%.3f", $tf->tag('gc_content'));

        $text .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n",
            $tf->name,
            $tf->ID,
            $tf->class() || 'NA',
            $tf->tag('family') || 'NA',
            $tf->tag('tax_group') || 'NA',
            $total_ic,
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
    
    unless (open(FH, ">$filename")) {
        fatal("Unable to create results text file $filename", %job_args);
        return;
    }

    print FH $text;

    close(FH);
    
    return $filename;
}


#
# Write the details of the putative TFBSs for each TF/seq. Create a
# text file for each TF.
#
sub write_tfbs_details_text
{
    my ($filename, $tf, $genes, $gid_tfbss,
        $gid_gene_ids, $gene_id_type,
        $has_operon, $operon_genes, %job_args) = @_;

    my $total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $total_ic = 'NA';
    }

    my $text = sprintf("%s\n\n", $tf->name());

    $text .= sprintf("JASPAR ID:\t%s\n", $tf->ID());
    $text .= sprintf("Class:\t%s\n", $tf->class() || 'NA');
    $text .= sprintf("Family:\t%s\n", $tf->tag('family') || 'NA');
    $text .= sprintf("Tax group:\t%s\n", $tf->tag('tax_group') || 'NA');
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
    foreach my $gene (@$genes)
    {
        # fetch gene and promoter information
        my ($gene_text, $gene_start, $tss, $promoters, $strand)
            = _fetch_gene_information_text(
                $gene, $gid_gene_ids, $gene_id_type, $has_operon, $operon_genes
            );
        
        $text .= $gene_text;
        
        my $tfbs_text = _fetch_site_information_text(
            $gene->id, $gid_tfbss, $gene_start, $tss, $promoters, $strand,
            $gene_id_type, $has_operon
        );
        
        $text .= $tfbs_text;
        
    } # end foreach gene    
    
    unless (open(FH, ">$filename")) {
        fatal("Unable to create TFBS details file $filename", %job_args);
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
    my ($filename, $rel_results_dir, $species, $tf, $genes, $gid_tfbss,
        $tf_db, $gene_id_type, $gid_gene_ids, $has_operon, $operon_genes,
        %job_args) = @_;

    my $tf_name = $tf->name();
    my $tf_id   = $tf->ID();

    my $job_id = $job_args{-job_id};
    my $heading = $job_args{-heading};
    my $email = $job_args{-email};
    my $logger = $job_args{-logger};
    
    open(FH, ">$filename") || fatal(
        "Could not create TFBS details html file $filename", %job_args
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    my $title = "oPOSSUM $heading";
    my $section = sprintf("%s Conserved Binding Sites", $tf_name);

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
        tf                  => $tf,
        gene_id_type        => $gene_id_type,
        dflt_gene_id_type   => DFLT_GENE_ID_TYPE,
        genes               => $genes,
        gid_gene_ids        => $gid_gene_ids,
        gid_tfbss           => $gid_tfbss,
        has_operon          => $has_operon,
        operon_genes        => $operon_genes,
        rel_results_dir     => $rel_results_dir,
        tfbs_details_file   => "$tf_id.txt",
        var_template        => "tfbs_details_gene_ssa.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}

1;
