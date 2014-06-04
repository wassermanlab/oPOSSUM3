use strict;

use OPOSSUM::Include::SeqInclude;
use OPOSSUM::Include::ACSAInclude;

use constant BG_COLOR_CLASS => 'bgc_seq_acsa';

#
# Search seqs with the anchoring TF and all the TFs in the TF set.
#
# Return two hashrefs:
#
# The first hashref refers to a hash of hashes of a list of sites keyed on
# TF ID and sequence ID. In the case where there were no TFBSs for a TF/seq
# combo, the value of the hash is undef.
#
# The second hasref refers to a hash of hashes of a list of site pairs.
#
# i.e.:
# $hashref1->{$tf_id}->{$seq_id} = $sites
# $hashref2->{$tf_id}->{$seq_id} = $site_pairs
#
sub anchored_tf_set_search_seqs
{
    my ($tf_set, $anchor_matrix, $max_site_dist, 
		$seq_id_seqs, $threshold, %job_args) = @_;

	my $logger = $job_args{-logger};

    my $tf_ids = $tf_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

    my $anchor_pwm;
    if ($anchor_matrix->isa("TFBS::Matrix::PFM")) {
        $anchor_pwm = $anchor_matrix->to_PWM();
    } else {
        $anchor_pwm = $anchor_matrix;
    }

    #
    # If threshold is specified as a decimal, convert it to a
    # percentage, otherwise the TFBS::Matrix::PWM::search_seq method
    # treats the number as an absolute matrix score which is not what
    # we intended. DJA 2012/06/07
    #
    unless ($threshold =~ /(.+)%$/ || $threshold > 1) {
        $threshold *= 100;
        $threshold .= '%';
    }

    my %tf_seq_sites;
    my %tf_seq_sitepairs;
    foreach my $seq_id (@seq_ids) {
        $logger->debug("\nSequence: $seq_id\n");

        my $seq = $seq_id_seqs->{$seq_id};

        my $anchor_siteset = $anchor_pwm->search_seq(
            -seqobj     => $seq,
            -threshold  => $threshold
        );

        #logprint_siteset("Anchor siteset", $anchor_siteset);

        next if !$anchor_siteset || $anchor_siteset->size() == 0;

        my $filtered_anchor_siteset = filter_overlapping_sites($anchor_siteset);

        #logprint_siteset("Filtered anchor siteset", $filtered_anchor_siteset);

        next if !$filtered_anchor_siteset
            || $filtered_anchor_siteset->size() == 0;

        foreach my $tf_id (@$tf_ids) {
            #next if $tf_id eq $anchor_matrix->ID;

            my $matrix = $tf_set->get_matrix($tf_id);

            my $pwm;
            if ($matrix->isa("TFBS::Matrix::PFM")) {
                $pwm = $matrix->to_PWM();
            } else {
                $pwm = $matrix;
            }

            my $siteset = $pwm->search_seq(
                -seqobj     => $seq,
                -threshold  => $threshold
            );

            #logprint_siteset("TF siteset", $siteset);

            next if !$siteset || $siteset->size == 0;

            my $filtered_siteset = filter_overlapping_sites($siteset);

            #logprint_siteset("Filtered TF siteset", $filtered_siteset);

            next if !$filtered_siteset || $filtered_siteset->size == 0;

            my ($prox_sites, $sitepairs) = proximal_sites(
                $filtered_anchor_siteset, $filtered_siteset, $max_site_dist
            );

            #logprint_siteset("Proximal siteset", $prox_siteset);
            #logprint_sitepairs($sitepairs);

            if ($prox_sites) {
                $tf_seq_sites{$tf_id}->{$seq_id} = $prox_sites;
                $tf_seq_sitepairs{$tf_id}->{$seq_id} = $sitepairs;
            }
        }
    }

    my $retval1 = %tf_seq_sites ? \%tf_seq_sites : undef;
    my $retval2 = %tf_seq_sitepairs ? \%tf_seq_sitepairs : undef;

    return ($retval1, $retval2);
}

#
# Same idea as above except it takes a list of anchoring matrices instead of
# just a single anchor matrix.
#
# The returned hash refs have an extra level keyed on the anchor TF ID.
# i.e.:
# $hashref1->{$anchor_tf_id}->{$tf_id}->{$seq_id} = $sites
# $hashref2->{$anchor_tf_id}->{$tf_id}->{$seq_id} = $site_pairs
#
sub multi_anchored_tf_set_search_seqs
{
    my ($tf_set, $anchor_matrices, $max_site_dist, 
		$seq_id_seqs, $threshold, %job_args) = @_;

	my $logger = $job_args{-logger};

    my $tf_ids = $tf_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

    my @anchor_pwms;
    foreach my $anchor_matrix (@$anchor_matrices) {
        if ($anchor_matrix->isa("TFBS::Matrix::PFM")) {
            push @anchor_pwms, $anchor_matrix->to_PWM();
        } else {
            push @anchor_pwms, $anchor_matrix;
        }
    }

    #
    # If threshold is specified as a decimal, convert it to a
    # percentage, otherwise the TFBS::Matrix::PWM::search_seq method
    # treats the number as an absolute matrix score which is not what
    # we intended. DJA 2012/06/07
    #
    unless ($threshold =~ /(.+)%$/ || $threshold > 1) {
        $threshold *= 100;
        $threshold .= '%';
    }

    my %tf_seq_sites;
    my %tf_seq_sitepairs;
    foreach my $seq_id (@seq_ids) {
        $logger->debug("\nSequence: $seq_id\n");

        my $seq = $seq_id_seqs->{$seq_id};

        my %anchor_siteset;
        foreach my $anchor_pwm (@anchor_pwms) {
            my $siteset = $anchor_pwm->search_seq(
                -seqobj     => $seq,
                -threshold  => $threshold
            );

            next if !$siteset || $siteset->size() == 0;

            my $filtered_anchor_siteset = filter_overlapping_sites(
                $siteset
            );

            next if !$filtered_anchor_siteset
                || $filtered_anchor_siteset->size() == 0;

            $anchor_siteset{$anchor_pwm->ID} = $siteset;
        }

        foreach my $tf_id (@$tf_ids) {
            #next if $tf_id eq $anchor_matrix->ID;

            my $matrix = $tf_set->get_matrix($tf_id);

            my $pwm;
            if ($matrix->isa("TFBS::Matrix::PFM")) {
                $pwm = $matrix->to_PWM();
            } else {
                $pwm = $matrix;
            }

            my $siteset = $pwm->search_seq(
                -seqobj     => $seq,
                -threshold  => $threshold
            );

            next if !$siteset || $siteset->size == 0;

            my $filtered_siteset = filter_overlapping_sites($siteset);

            next if !$filtered_siteset || $filtered_siteset->size == 0;

            foreach my $anchor_tf_id (keys %anchor_siteset) {
                my ($prox_sites, $sitepairs) = proximal_sites(
                    $anchor_siteset{$anchor_tf_id},
                    $filtered_siteset,
                    $max_site_dist
                );

                if ($prox_sites) {
                    $tf_seq_sites{$anchor_tf_id}->{$tf_id}->{$seq_id}
                        = $prox_sites;
                    $tf_seq_sitepairs{$anchor_tf_id}->{$tf_id}->{$seq_id}
                        = $sitepairs;
                }
            }
        }
    }

    my $retval1 = %tf_seq_sites ? \%tf_seq_sites : undef;
    my $retval2 = %tf_seq_sitepairs ? \%tf_seq_sitepairs : undef;

    return ($retval1, $retval2);
}

sub compute_site_counts
{
    my ($tf_set, $seq_ids, $tf_seq_sites) = @_;

    my $tf_ids  = $tf_set->ids();

    my $counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids       => $seq_ids,
        -tf_ids         => $tf_ids
    );

    foreach my $seq_id (@$seq_ids) {
        foreach my $tf_id (@$tf_ids) {
            my $sites = $tf_seq_sites->{$tf_id}->{$seq_id};

            if ($sites) {
                # Note set size could be 0.
                $counts->gene_tfbs_count($seq_id, $tf_id, scalar @$sites);
            } else {
                $counts->gene_tfbs_count($seq_id, $tf_id, 0);
            }
        }
    }

    return $counts;
}

#
# As above but works on multiple anchor TF IDs. Note for this the $tf_seq_sites
# is actually a hash ref with the structure:
#     $tf_seq_sites->{$anchor_tf_id}->{$tf_id}->{$seq_id}->sites
#
# This routine is not really necessary. Just call above routine in a loop.
#
sub compute_site_counts_multi_anchor
{
    my ($tf_set, $anchor_tf_ids, $seq_ids, $tf_seq_sites) = @_;

    my $tf_ids  = $tf_set->ids();

    my %anchor_tf_counts;
    foreach my $anchor_tf_id (@$anchor_tf_ids) {
        my $counts = OPOSSUM::Analysis::Counts->new(
            -gene_ids       => $seq_ids,
            -tf_ids         => $tf_ids
        );

        foreach my $seq_id (@$seq_ids) {
            foreach my $tf_id (@$tf_ids) {
                my $sites = $tf_seq_sites->{$anchor_tf_id}->{$tf_id}->{$seq_id};

                if ($sites) {
                    # Note set size could be 0.
                    $counts->gene_tfbs_count($seq_id, $tf_id, scalar @$sites);
                } else {
                    $counts->gene_tfbs_count($seq_id, $tf_id, 0);
                }
            }
        }

        $anchor_tf_counts{$anchor_tf_id} = $counts;
    }

    return %anchor_tf_counts ? \%anchor_tf_counts : undef;
}


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

    #
    # Single line (tab-delimited) header format
    # NOTE: rearranged columns
    #
    my $text = "TF Name\tTF ID\tClass\tFamily\tTax Group\tIC\tGC Content\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

    foreach my $result (@$results) {
        my $tf = $tf_set->get_tf($result->id());

        my $total_ic;
        if ($tf->isa("TFBS::Matrix::PFM")) {
            $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
        } else {
            $total_ic = 'NA';
        }

        my $gc_content = 'NA';
        if ($tf->tag('gc_content')) {
            $gc_content = sprintf("%.3f", $tf->tag('gc_content'));
        }

        $text .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",
            $tf->name(),
            $tf->ID() || 'NA',
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

	print FH $text . "\n";

    close(FH);

    return $filename;
}



#
# For each TF/seq, write the details of the putative TFBSs out to text and
# html files.
#
sub write_tfbs_details
{
    my ($tf_seq_sitepairs, $seq_id_display_ids, $tf_set, $anchor_matrix,
	   	$abs_results_dir, $rel_results_dir, $web, %job_args) = @_;

    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $seq_sitepairs = $tf_seq_sitepairs->{$tf_id};

        if ($seq_sitepairs) {
            my $tf = $tf_set->get_tf($tf_id);

            my @seq_ids = keys %$seq_sitepairs;

            my $fname = $tf_id;
            $fname =~ s/\//_/g;

            my $text_filename = "$abs_results_dir/$fname.txt";
            my $html_filename = "$abs_results_dir/$fname.html";
        
            write_tfbs_details_text(
                $text_filename, $tf, $anchor_matrix, 
                \@seq_ids, $seq_id_display_ids, $seq_sitepairs,
				%job_args
            );

            write_tfbs_details_html(
                $html_filename, $rel_results_dir, $tf, $anchor_matrix,
                \@seq_ids, $seq_id_display_ids, $seq_sitepairs,
				%job_args
            ) if $web;
        }
    }
}

#
# Write the details of the putative TFBSs for each TF/seq. Create a
# text file for each TF.
#
sub write_tfbs_details_text
{
    my ($filename, $tf, $anchor_matrix, $seq_ids, $seq_id_display_ids,
        $seq_sitepairs, %job_args) = @_;

	my $logger = $job_args{-logger};

    my $tf_db = $job_args{-tf_db};
    my $anchor_tf_db = $job_args{-anchor_tf_db};

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details text file $filename", %job_args
    );

    my $tf_id       = $tf->ID();
    my $tf_name     = $tf->name();
    my $anchor_id   = $anchor_matrix->ID();
    my $anchor_name = $anchor_matrix->name();

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

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    print  FH "$anchor_name\n\n";

    printf FH "TF ID:\t$anchor_id\n";
    printf FH "Class:\t%s\n", $anchor_matrix->class() || 'NA',
    printf FH "Family:\t%s\n", $anchor_matrix->tag('family') || 'NA',
    printf FH "Tax group:\t%s\n", $anchor_matrix->tag('tax_group') || 'NA',
    printf FH "IC:\t%s\n", $anchor_total_ic;

    print FH "\n\n";

    print  FH "$tf_name\n\n";

    print  FH "TF ID:\t$tf_id\n";
    printf FH "Class:\t%s\n", $tf->class() || 'NA',
    printf FH "Family:\t%s\n", $tf->tag('family') || 'NA',
    printf FH "Tax group:\t%s\n", $tf->tag('tax_group') || 'NA',
    printf FH "IC:\t%s\n", $tf_total_ic;

    print FH "\n\n$anchor_name - $tf_name Binding Site Combinations\n\n";

    #printf FH "\n\n%-31s\t%s\t%s\t%s\t%s\t%s\%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
    printf FH "\n\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        'Sequence ID', 'Anchoring TF', 'Start', 'End', 'Strand', 'Score', '%Score', 'TFBS Sequence', 'Anchored TF', 'Start', 'End', 'Strand', 'Score', '%Score', 'TFBS Sequence', 'Distance';

    foreach my $seq_id (@$seq_ids) {
        my $sitepairs = $seq_sitepairs->{$seq_id};
        next unless $sitepairs && $sitepairs->[0];

        my $display_id = $seq_id_display_ids->{$seq_id};

        my $first = 1;
        #printf FH "%-31s", $display_id;
        printf FH "%s", $display_id;

        foreach my $pair (@$sitepairs) {
            my $anchor_site = $pair->{anchor_site};
            my $tf_site     = $pair->{tf_site};
            my $distance    = $pair->{distance};

            # Do not reprint sequence ID for subsequent sites
            #unless ($first) {
            #    #printf FH "%-31s", '';
            #    printf FH "\t";
            #}
            $first = 0;

            #printf FH "\t%s\t%7d\t%7d\t%7d\t%7.3f\t%6.1f%%\t%s\t%s\t%7d\t%7d\t%7d\t%7.3f\t%6.1f%%\t%s\t%7d\n",
            printf FH "\t%s\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\t%s\t%d\t%d\t%s\t%.3f\t%.1f%%\t%s\t%d\n",
                $anchor_name,
                $anchor_site->{start},
                $anchor_site->{end},
                $anchor_site->{strand} == 1 ? '+' : '-',
                $anchor_site->{score},
                $anchor_site->{rel_score} * 100,
                $anchor_site->{seq},
                $tf_name,
                $tf_site->{start},
                $tf_site->{end},
                $tf_site->{strand} == 1 ? '+' : '-',
                $tf_site->{score},
                $tf_site->{rel_score} * 100,
                $tf_site->{seq},
                $distance;
        }
    }

    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/seq. Create an
# HTML file for each TF.
#
sub write_tfbs_details_html
{
    my ($filename, $rel_results_dir, $tf, $anchor_matrix,
	   	$seq_ids, $seq_id_display_ids, $seq_sitepairs,
		%job_args) = @_;

	my $heading = $job_args{-heading};
	my $title = $job_args{-title};
	my $logger = $job_args{-logger};

    my $tf_id   = $tf->ID();
    my $tf_name = $tf->name();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details HTML file $filename", %job_args
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    my $section = sprintf "%s - %s Binding Site Combinations",
        $anchor_matrix->name(), $tf_name;

    my $vars = {
        abs_htdocs_path         => ABS_HTDOCS_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_htdocs_path         => REL_HTDOCS_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        version                 => VERSION,
        devel_version           => DEVEL_VERSION,
        jaspar_url              => JASPAR_URL,
        heading                 => $heading,
        title                   => $title,
        section                 => $section,
        bg_color_class          => BG_COLOR_CLASS,
        low_matrix_ic           => LOW_MATRIX_IC,
        high_matrix_ic          => HIGH_MATRIX_IC,
        low_matrix_gc           => LOW_MATRIX_GC,
        high_matrix_gc          => HIGH_MATRIX_GC,
        low_seq_gc              => LOW_SEQ_GC,
        high_seq_gc             => HIGH_SEQ_GC,


        formatf                 => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'NA'
                                   },

        anchor_tf_db            => $job_args{-anchor_tf_db},
        tf_db                   => $job_args{-tf_db},
        tf                      => $tf,
        anchor_matrix           => $anchor_matrix,
        seq_ids                 => $seq_ids,
        seq_id_display_ids      => $seq_id_display_ids,
        #seq_sites               => \%seq_sites,
        seq_sitepairs           => $seq_sitepairs,
        rel_results_dir         => $rel_results_dir,
        tfbs_details_file       => "$tf_id.txt",
        var_template            => "tfbs_details_seq_acsa.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}




1;

