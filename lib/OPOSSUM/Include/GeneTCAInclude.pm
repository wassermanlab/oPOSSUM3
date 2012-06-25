=head1 NAME

 OPOSSUM::Include::GeneTCAInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the gene-based TCA
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
use OPOSSUM::Include::TCAInclude;

use constant BG_COLOR_CLASS => 'bgc_gene_tca';

sub write_results_text
{
    my ($filename, $results, $tf_cluster_set, %job_args) = @_;

    return unless $results && $results->[0];

    my $text = "TFBS Cluster ID\tClass\tFamily";
    $text .= "\tTarget gene hits\tTarget gene non-hits";
    $text .= "\tBackground gene hits\tBackground gene non-hits";
    $text .= "\tTarget cluster hits\tTarget cluster nucleotide rate";
    $text .= "\tBackground cluster hits\tBackground cluster nucleotide rate";
    $text .= "\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $cl = $tf_cluster_set->get_tf_cluster($result->id());

        $text .= sprintf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",
            'c' . $cl->id(),
            $cl->class() || 'NA',
            $cl->family() || 'NA',
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_cluster_hits() || 0,
            defined $result->t_cluster_rate()
                ? sprintf("%.3g", $result->t_cluster_rate()) : 'NA',
            $result->bg_cluster_hits() || 0,
            defined $result->bg_cluster_rate()
                ? sprintf("%.3g", $result->bg_cluster_rate()) : 'NA',
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

sub write_tfbs_cluster_details_text
{
    my ($filename, $tfcl, $genes, $gid_cluster_sites,
        $gid_gene_ids, $gene_id_type,
        $has_operon, $operon_genes, %job_args) = @_;

    my $text = sprintf("C%s\n\n", $tfcl->id());
    $text .= sprintf("Class:\t%s\n", $tfcl->class || 'NA');
    $text .= sprintf("Family:\t%s\n", $tfcl->family || 'NA');

    $text .= sprintf("\n\nConserved C%s Binding Sites:\n\n", $tfcl->id());

    unless (defined $gene_id_type && $gene_id_type == DFLT_GENE_ID_TYPE) {
        $text .= "Gene ID(s)\t";
    }

    $text .= "Ensembl ID(s)";

    if ($has_operon) {
        $text .= "\tOperon ID";
    }

    $text .= qq{\tChr\tProm.Start\tProm.End\tProm.Strand\tNearest TSS\t};
    $text .= qq{Site Start\tSite End\tRel. Start\tRel. End\tStrand\t%%Score\tSequence\n};

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
        
        my $sites_text = _fetch_site_information_text(
            $gene->id, $gid_cluster_sites,
            $gene_start, $tss, $promoters, $strand,
            $gene_id_type, $has_operon
        );
        
        $text .= $sites_text;
        
    } # end foreach gene    
    
    unless (open(FH, ">$filename")) {
        fatal("Unable to create TFBS cluster details file $filename", %job_args);
        return;
    }

    print FH $text;

    close(FH);
}

#
# Write the details of the putative TFBS clusters for each cluster/gene. Create an
# HTML file for each TFBS cluster.
#
sub write_tfbs_cluster_details_html
{
    my ($filename, $rel_results_dir, $species,
        $tfcl, $genes, $gid_cluster_sites,
        $gene_id_type, $gid_gene_ids,
        $has_operon, $operon_genes,
        %job_args) = @_;

    my $job_id = $job_args{-job_id};
    my $heading = $job_args{-heading};
    my $email = $job_args{-email};
    my $logger = $job_args{-logger};
    
    my $tfcl_id   = $tfcl->id();
    
    open(FH, ">$filename") || fatal(
        "Could not create TFBS cluster details html file $filename",
        %job_args
    );

    #$logger->info("Writing 'c$tfcl_id' TFBS cluster details to $filename");

    my $title = "oPOSSUM $heading";
    my $section = sprintf("%s Conserved Binding Sites", "c$tfcl_id");
    
    my %gid_cluster_sites;
    foreach my $gid (keys %$gid_cluster_sites) {
        my $sites = $gid_cluster_sites->{$gid};
        next if !$sites || scalar @$sites == 0;
        $gid_cluster_sites{$gid} = $sites;
    }

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

        tf_cluster          => $tfcl,
        gene_id_type        => $gene_id_type,
        dflt_gene_id_type   => DFLT_GENE_ID_TYPE,
        genes               => $genes,
        gid_gene_ids        => $gid_gene_ids,
        gid_cluster_tfbss   => \%gid_cluster_sites,
        has_operon          => $has_operon,
        operon_genes        => $operon_genes,
        rel_results_dir     => $rel_results_dir,
        tfbs_cluster_details_file   => "c$tfcl_id.txt",
        var_template        => "tfbs_details_gene_tca.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}

1;
