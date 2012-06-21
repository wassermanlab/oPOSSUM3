=head1 NAME

 OPOSSUM::Include::SeqSSAInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the sequence-based SSA
  scripts and web modules.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use OPOSSUM::Include::SeqInclude;

use OPOSSUM::Analysis::Values;

use constant BG_COLOR_CLASS => 'bgc_seq_ssa';

use strict;


#
# Search seqs with all the TFs in the TF set. Return two hashrefs.
#
# The first hashref refers to a hash of hashes of sites indexed on seq ID
# and TF ID. In the case where there were no TFBSs for a seq/TF combo,
# the value of the hash is undef.
#
# The second hasref refers to a hash of listrefs of sequence IDs, indexed
# on TF ID. Only seq IDs for which there were TFBSs for this TF are stored
# in the list.
#
# i.e.:
# $hashref1->{$seq_id}->{$tf_id} = $siteset
# $hashref2->{$tf_id} = \@seq_ids
#
sub tf_set_search_seqs
{
    my ($tf_set, $seq_id_seqs, $threshold, %job_args) = @_;

    my $logger = $job_args{-logger};

    my $tf_ids = $tf_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

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

    my %tf_seq_siteset;

	foreach my $tf_id (@$tf_ids) {
	    my $matrix = $tf_set->get_matrix($tf_id);

        my $pwm;
        if ($matrix->isa("TFBS::Matrix::PFM")) {
            $pwm = $matrix->to_PWM();
        } else {
            $pwm = $matrix;
        }

        foreach my $seq_id (@seq_ids) {
            my $seq = $seq_id_seqs->{$seq_id};

            #$logger->info("processing sequence $seq_id");

            my $siteset = $pwm->search_seq(
                -seqobj     => $seq,
                -threshold  => $threshold
            );

            my $filtered_siteset;
            if ($siteset && $siteset->size > 0) {
                $filtered_siteset = filter_overlapping_sites($siteset);

                if ($filtered_siteset && $filtered_siteset->size > 0) {
                    # Only set if seqs had sites for this TF
                    #push @{$tf_seqs{$tf_id}}, $seq_id;
                    $tf_seq_siteset{$tf_id}->{$seq_id} = $filtered_siteset;
                }
            }
        }
    }

    my $retval1 = %tf_seq_siteset ? \%tf_seq_siteset : undef;
    #my $retval2 = %tf_seqs ? \%tf_seqs : undef;

    #return ($retval1, $retval2);
    return $retval1;
}

sub compute_site_counts
{
    my ($tf_set, $seq_ids, $tf_seq_siteset) = @_;

    my $tf_ids  = $tf_set->ids();

    my $counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids       => $seq_ids,
        -tf_ids         => $tf_ids
    );

    foreach my $seq_id (@$seq_ids) {
        foreach my $tf_id (@$tf_ids) {
            my $siteset = $tf_seq_siteset->{$tf_id}->{$seq_id};

            if ($siteset) {
                # Note set size could be 0.
                $counts->gene_tfbs_count($seq_id, $tf_id, $siteset->size());
            } else {
                $counts->gene_tfbs_count($seq_id, $tf_id, 0);
            }
        }
    }

    return $counts;
}

=head2 compute_tfbs_peak_distances

 Title   : compute_tfbs_peak_distances

 Function: Calculate the distance of each TFBS from the peak max height
           location and return an array of distances. For the purpose of
           distance computation, use the center of the TFBS. If
           seq_peak_max_pos is provided, use this as the location of the peak
           maxes otherwise assume the center point of the peak is the max peak
           height location.

 Args    : tf_set           - an OPOSSUM::TFSet object
           seq_id_seqs      - hashref where key is the sequence ID and the
                              value is the sequence object
           tf_seq_site      - hashref where the key is the TF ID and the
                              value is a hashref where key is the sequence
                              ID and value is an arrayref of hashes
                              describing TFBSs - NOTE these are NOT
                              TFBS::Site objects
           seq_peak_max_pos - Optional hashref where key is sequence ID and
                              value is the position of the sequence peak
                              height max

 Returns : A hashref of TFBS distances from the peak height maximum keyed on
           the TF ID.

 Note: for KS test, should normalize to (0,1).
 
=cut

sub compute_tfbs_peak_distances
{
    my ($tf_set, $seq_id_seqs, $seq_id_display_ids, $tf_seq_siteset, $seq_peak_max_pos) = @_;

    my $tf_ids = $tf_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

    #my %distances;
    my $distances = OPOSSUM::Analysis::Values->new(
        -seq_ids => \@seq_ids,
        -tf_ids => $tf_ids,
    );
    
    foreach my $tf_id (@$tf_ids) {
        #$distances{$tf_id} = [];

        foreach my $seq_id (@seq_ids) {
            my $siteset = $tf_seq_siteset->{$tf_id}->{$seq_id};

            if ($siteset) {
                my $seq = $seq_id_seqs->{$seq_id};
				
                my $display_id = $$seq_id_display_ids{$seq_id};

				my $peak_max_loc = 0;
                
                # if display_id does not contain coordinates, cannot proceed
                if ($display_id =~ /^(\w+):(\d+)-(\d+)/ and $seq_peak_max_pos) {
				    
					# convert peak_max_loc to sequence frame
					my ($chr, $pos) = split ":", $display_id;
					my ($seq_start, $seq_end) = split "-", $pos;
					
                    $peak_max_loc = $seq_peak_max_pos->{$seq_id} - $seq_start + 1;
                } else {
                    # Assume centre of peak sequence is the peak max location
                    $peak_max_loc = $seq->length / 2;
                }
                
                my $iter = $siteset->Iterator(-sort_by => 'start');
                while (my $site = $iter->next) {
#                foreach my $site (@$sites) {0
                    #print $site->pattern->ID . "\t" . $site->start . "\t" . $site->seq->length . "\n";
                    #my $site_start  = $site->{start};
                    #my $site_length = $site->{seq}->length;

                    my $site_loc = $site->start + $site->seq->length / 2;

                    my $dist = abs($site_loc - $peak_max_loc);
                    my $dist_n = $dist / $seq->length;
                    
                    $distances->add_seq_tfbs_value($seq_id, $tf_id, $dist_n);
                    #push @{$distances{$tf_id}}, $dist;
                }
            }
        }
    }

    #return %distances ? \%distances : undef;
    return $distances;
}

#
# Ouput combined Z-score/Fisher/KS-test results to a text file
#
sub write_results_text
{
    my ($filename, $results, $tf_set, $tf_file, $job_args) = @_;

    return unless $results && $results->[0];

    #my $filename = "$results_dir/" . RESULTS_TEXT_FILENAME;
    my $logger = $$job_args{-logger};
    $logger->info("Writing analysis results to $filename\n");

    my $text;
    if (!$tf_file)
    {
        #
        # Single line (tab-delimited) header format
        # NOTE: rearranged columns
        #
        $text = "TF\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tGC Content\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\tKS score\n";

        foreach my $result (@$results) {
            my $tf = $tf_set->get_tf($result->id());

            my $total_ic;
            if ($tf->isa("TFBS::Matrix::PFM")) {
                $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
            } else {
                $total_ic = 'NA';
            }
            
            my $gc_content = sprintf("%.3f", $tf->tag('gc_content'));

            $text .= sprintf 
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
                $tf->name(),
                $tf->ID(),
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
                    ? sprintf("%.3f", $result->t_tfbs_rate()) : 'NA',
                $result->bg_tfbs_hits() || 0,
                defined $result->bg_tfbs_rate()
                    ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'NA',
                defined $result->zscore()
                    ? sprintf("%.3f", $result->zscore()) : 'NA',
                defined $result->fisher_p_value()
                    ? sprintf("%.3f", $result->fisher_p_value()) : 'NA',
                defined $result->ks_p_value()
                    ? sprintf("%.3f", $result->ks_p_value()) : 'NA';
        }
    } else {
        $text = "TF Name\tClass\tFamily\tTax Group\tIC\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\tPeakDist p-value\n";

        foreach my $result (@$results) {
            my $tf = $tf_set->get_tf($result->id());

            my $total_ic;
            if ($tf->isa("TFBS::Matrix::PFM")) {
                $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
            } else {
                $total_ic = 'NA';
            }

            $text .= sprintf 
                "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
                $tf->name(),
                $tf->class() || 'NA',
                $tf->tag('family') || 'NA',
                $tf->tag('tax_group') || 'NA',
                $total_ic,
                $result->t_gene_hits() || 0,
                $result->t_gene_no_hits() || 0,
                $result->bg_gene_hits() || 0,
                $result->bg_gene_no_hits() || 0,
                $result->t_tfbs_hits() || 0,
                defined $result->t_tfbs_rate()
                    ? sprintf("%.3f", $result->t_tfbs_rate()) : 'NA',
                $result->bg_tfbs_hits() || 0,
                defined $result->bg_tfbs_rate()
                    ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'NA',
                defined $result->zscore()
                    ? sprintf("%.3f", $result->zscore()) : 'NA',
                defined $result->fisher_p_value()
                    ? sprintf("%.3f", $result->fisher_p_value()) : 'NA',
                defined $result->ks_p_value()
                    ? sprintf("%.3f", $result->ks_p_value()) : 'NA';
        }
    }
    
    unless (open(FH, ">$filename")) {
        fatal("Unable to create results text file $filename", %$job_args);
        return;
    }
    
    print FH $text;
    
    close(FH);

    return $filename;
}



#
# For each TF/seq, write the details of the putative TFBSs out to text and
# html files.
#
sub write_tfbs_details
{
    my ($tf_set, $tf_db, $tf_seq_siteset,
        $seq_id_display_ids, 
        $abs_results_dir, $rel_results_dir, $web,
		%job_args) = @_;

    my $logger = $job_args{-logger};
	#$logger->info("Writing TFBS details");
	
    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $seq_siteset = $tf_seq_siteset->{$tf_id};

        if ($seq_siteset) {
            my @seq_ids = keys %$seq_siteset;

            my $tf = $tf_set->get_tf($tf_id);

            my $text_filename = sprintf "$abs_results_dir/$tf_id.txt";
            my $html_filename = sprintf "$abs_results_dir/$tf_id.html";
        
            write_tfbs_details_text(
                $text_filename,
                $tf, $tf_db,
                \@seq_ids,
                $seq_id_display_ids,
                $seq_siteset,
                %job_args
            );
            
            write_tfbs_details_html(
                $html_filename, $rel_results_dir,
                $tf, $tf_db,
                \@seq_ids,
                $seq_id_display_ids,
                $seq_siteset,
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
    my ($filename, $tf, $tf_db,
        $seq_ids, $seq_id_display_ids, $seq_siteset,
        %job_args) = @_;
    
    my $logger = $job_args{-logger};
    
    my $tf_id   = $tf->ID();
    my $tf_name = $tf->name();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details text file $filename", %job_args
    );

    #$logger->info("Writing '$tf_name' TFBS details to $filename");

    my $total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $total_ic = 'NA';
    }

    print  FH "$tf_name\n\n";

    if ($tf_db) {
        print  FH "JASPAR ID:\t$tf_id\n";
    }
    printf FH "Class:\t%s\n", $tf->class() || 'NA',
    printf FH "Family:\t%s\n", $tf->tag('family') || 'NA',
    printf FH "Sysgroup:\t%s\n", $tf->tag('tax_group') || 'NA',
    printf FH "IC:\t%s\n", $total_ic;

    print FH "\n\n$tf_name Binding Sites\n\n";

    printf FH "\n\n%-31s\t%7s\t%7s\t%7s\t%7s\t%7s\t%s\n",
        'Sequence ID', 'Start', 'End', 'Strand', 'Score', '%Score', 'TFBS Sequence';

    foreach my $seq_id (@$seq_ids) {
        my $siteset = $seq_siteset->{$seq_id};

        my $display_id = $seq_id_display_ids->{$seq_id};

        next if !$siteset || $siteset->size == 0;

        my $first = 1;
        printf FH "%-31s", $display_id;

        my $iter = $siteset->Iterator(-sort_by => 'start');
        while (my $site = $iter->next()) {
            # Do not reprint sequence ID for subsequent sites
            unless ($first) {
                printf FH "%-31s", '';
            }
            $first = 0;

            #printf FH "%s\t%d\t%d\t%d\t%.3f\t%.1f%%\t%s\n",
            printf FH "\t%7d\t%7d\t%7d\t%7.3f\t%6.1f%%\t%s\n",
                $site->start(),
                $site->end(),
                $site->strand(),
                $site->score(),
                $site->rel_score() * 100,
                $site->seq->seq();
        }

        print FH "\n";
    }

    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/seq. Create an
# HTML file for each TF.
#
sub write_tfbs_details_html
{
    my ($filename, $rel_results_dir,
        $tf, $tf_db,
        $seq_ids, $seq_id_display_ids, $seq_siteset,
        %job_args) = @_;
    
    my $job_id = $job_args{-job_id};
    my $heading = $job_args{-heading};
    my $email = $job_args{-email};
    my $logger = $job_args{-logger};

    my $tf_id   = $tf->ID();
    my $tf_name = $tf->name();

    my %seq_sites;
    foreach my $seq_id (@$seq_ids) {
        my $siteset = $seq_siteset->{$seq_id};

        next if !$siteset || $siteset->size == 0;

        my @site_list;
        my $iter = $siteset->Iterator(-sort_by => 'start');
        while (my $site = $iter->next()) {
            push @site_list, $site;
        }

        $seq_sites{$seq_id} = \@site_list;
    }

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details HTML file $filename", %job_args
    );

    #$logger->info("Writing '$tf_name' TFBS details to $filename");
    
    my $title = "oPOSSUM $heading";
    my $section = sprintf("%s Binding Site Details", $tf_name);

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

        tf_db                   => $tf_db,
        tf                      => $tf,
        seq_ids                 => $seq_ids,
        seq_id_display_ids      => $seq_id_display_ids,
        seq_sites               => \%seq_sites,
        rel_results_dir         => $rel_results_dir,
        tfbs_details_file       => "$tf_id.txt",
        var_template            => "tfbs_details_seq_ssa.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}

1;
