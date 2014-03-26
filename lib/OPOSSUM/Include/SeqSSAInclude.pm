=head1 NAME

 OPOSSUM::Include::SeqSSAInclude.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

  Contains all options and routines that are common to all the
  sequence-based SSA scripts and web modules.

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
# New routine to search sequences for binding sites, compute the TFBS counts
# and the site to peak max distances. The binding site destails are written
# to temporary data files for later conversion to text and html files to
# avoid having to store all the binding sites in memory.
# DJA 2012/11/07
#
# Search seqs with all the TFs in the TF set. Return two values, an
# OPOSSUM::Analysis::Counts object of the seq-TF binding site counts and an
# OPOSSUM::Analysis::Values object which stores the seq-TF binding site to
# peak max. distances.
#
sub search_seqs_compute_tfbs_counts_and_distances
{
    my (
        $tf_set, $seq_id_seqs, $seq_id_display_ids, $seq_peak_max_pos,
        $threshold, $write_details, %job_args
    ) = @_;

    my $any_tfbs_found = 0;

    my $logger = $job_args{-logger};
    my $results_dir = $job_args{-results_dir};

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

    my $counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids   => \@seq_ids,
        -tf_ids     => $tf_ids
    );

    my $distances = OPOSSUM::Analysis::Values->new(
        -seq_ids    => \@seq_ids,
        -tf_ids     => $tf_ids,
    );

	foreach my $tf_id (@$tf_ids) {
        my $fh;

        if ($write_details && $results_dir) {
            my $fname = $tf_id;
            $fname =~ s/\//_/g;

            my $data_file = "$results_dir/$fname.data";

            if (open($fh, ">$data_file")) {
                #
                # Specific header format recognized by Datafile plugin of the
                # Template Toolkit.
                #
                printf $fh
                    "seq_id : start : end : strand : score : rel_score : seq\n";
            } else {
                carp "Error opening output TFBS details data file"
                    . " $data_file - $!\n";
            }
        }

	    my $matrix = $tf_set->get_matrix($tf_id);

        my $pwm;
        if ($matrix->isa("TFBS::Matrix::PFM")) {
            $pwm = $matrix->to_PWM();
        } else {
            $pwm = $matrix;
        }

        foreach my $seq_id (sort @seq_ids) {
            my $tfbs_count = 0;

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
                    $any_tfbs_found = 1;

                    if ($fh) {
                        write_tfbs_details_data(
                            $fh, $seq_id, $filtered_siteset, %job_args
                        );
                    }

                    $tfbs_count = $filtered_siteset->size;

                    compute_seq_tf_site_peak_distances(
                        $distances, $seq_id, $tf_id,
                        $seq_id_display_ids->{$seq_id},
                        $seq, $seq_peak_max_pos->{$seq_id},
                        $filtered_siteset
                    );
                }
            }

            $counts->gene_tfbs_count($seq_id, $tf_id, $tfbs_count);
        }

        close($fh) if $fh;
    }

    return ($any_tfbs_found ? $counts : undef, $distances);
}

#
# New routine called by search_seqs_compute_tfbs_counts_and_distances to
# compute the binding site to peak max distances for a single TF/sequence
# combination.
# DJA 2012/11/07
#
sub compute_seq_tf_site_peak_distances
{
    my (
        $distances, $seq_id, $tf_id, $display_id, $seq, $peak_max_pos,
        $siteset)
    = @_;

    my $peak_max_loc = 0;
    
    # if display_id does not contain coordinates, cannot proceed
    if ($display_id =~ /^(\w+):(\d+)-(\d+)/ and $peak_max_pos) {
        
        # convert peak_max_loc to sequence frame
        my ($chr, $pos) = split ":", $display_id;
        my ($seq_start, $seq_end) = split "-", $pos;
        
        $peak_max_loc = $peak_max_pos->{$seq_id} - $seq_start + 1;
    } else {
        # Assume centre of peak sequence is the peak max location
        $peak_max_loc = $seq->length / 2;
    }
    
    my $iter = $siteset->Iterator(-sort_by => 'start');
    while (my $site = $iter->next) {
        my $site_loc = $site->start + $site->seq->length / 2;

        my $dist = abs($site_loc - $peak_max_loc);
        my $dist_n = $dist / $seq->length;
        
        $distances->add_seq_tfbs_value($seq_id, $tf_id, $dist_n);
    }
}

#
# Deprecated routine. The site counts is now performed within the
# search_seqs_compute_tfbs_counts_and_distances routine.
#
# DJA 2012/11/07
#
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


#
# Deprecated routine. The TFBS to peak max distances are now performed within
# the search_seqs_compute_tfbs_counts_and_distances routine.
#
# DJA 2012/11/07
#

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

#
#    This special case existed because if the profiles are custom profiles
#    from a user specified file ($tf_file is defined) the TF objects didn't
#    (necessarily) have an ID defined. Therefore the output should not include
#    an TF ID header/data column. This changed though so all TF objects
#    have the ID set (even if it just duplicates the name) and we should
#    probably maintain consistency in the output columns.
#
#    if (!$tf_file) {
        #
        # Single line (tab-delimited) header format
        # NOTE: rearranged columns
        #
        $text = "TF\tID\tClass\tFamily\tTax Group\tIC\tGC Content\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\tKS score\n";

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
                    ? sprintf("%.3f", $result->fisher_p_value()) : 'NA',
                defined $result->ks_p_value()
                    ? sprintf("%.3f", $result->ks_p_value()) : 'NA';
        }
#    } else {
#        $text = "TF Name\tClass\tFamily\tTax Group\tIC\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS hits\tBackground TFBS nucleotide rate\tZ-score\tFisher score\tKS score\n";
#
#        foreach my $result (@$results) {
#            my $tf = $tf_set->get_tf($result->id());
#
#            my $total_ic;
#            if ($tf->isa("TFBS::Matrix::PFM")) {
#                $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
#            } else {
#                $total_ic = 'NA';
#            }
#
#            $text .= sprintf 
#                "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
#                $tf->name(),
#                $tf->class() || 'NA',
#                $tf->tag('family') || 'NA',
#                $tf->tag('tax_group') || 'NA',
#                $total_ic,
#                $result->t_gene_hits() || 0,
#                $result->t_gene_no_hits() || 0,
#                $result->bg_gene_hits() || 0,
#                $result->bg_gene_no_hits() || 0,
#                $result->t_tfbs_hits() || 0,
#                defined $result->t_tfbs_rate()
#                    ? sprintf("%.3f", $result->t_tfbs_rate()) : 'NA',
#                $result->bg_tfbs_hits() || 0,
#                defined $result->bg_tfbs_rate()
#                    ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'NA',
#                defined $result->zscore()
#                    ? sprintf("%.3f", $result->zscore()) : 'NA',
#                defined $result->fisher_p_value()
#                    ? sprintf("%.3f", $result->fisher_p_value()) : 'NA',
#                defined $result->ks_p_value()
#                    ? sprintf("%.3f", $result->ks_p_value()) : 'NA';
#        }
#    }
    
    unless (open(FH, ">$filename")) {
        fatal("Unable to create results text file $filename", %$job_args);
        return;
    }
    
    print FH $text;
    
    close(FH);

    return $filename;
}

#
# This routine has been modified to call the write_tfbs_details_text_from_data
# and write_tfbs_details_html_from_data routines 
# DJA 2012/11/07
#
# For each TF/seq, write the details of the putative TFBSs out to text and
# html files.
#
sub write_tfbs_details
{
    my ($tf_set, $seq_id_display_ids, %job_args) = @_;

    my $logger      = $job_args{-logger};
    my $web         = $job_args{-web};
    my $results_dir = $job_args{-results_dir};

	#$logger->info("Writing TFBS details");
	
    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $tf = $tf_set->get_tf($tf_id);

        my $fname = $tf_id;
        $fname =~ s/\//_/g;

        my $data_file = "$results_dir/$fname.data";

        my $text_file = "$results_dir/$fname.txt";

        write_tfbs_details_text_from_data(
            $tf, $seq_id_display_ids, $data_file, $text_file, %job_args
        );

        if ($web) {
            my $html_file = "$results_dir/$fname.html";

            write_tfbs_details_html_from_data(
                $tf, $seq_id_display_ids, $data_file, $html_file, %job_args
            );
        }

        #
        # Remove data file
        #
        unlink $data_file;
    }
}

#
# Write the details of the putative TFBSs to the given file in a format
# recognized by the Datafile plugin of the Template Toolkit, for later
# conversion into text and html files.
#
sub write_tfbs_details_data
{
    my ($fh, $seq_id, $siteset) = @_;
    
    my $iter = $siteset->Iterator(-sort_by => 'start');

    while (my $site = $iter->next()) {
        printf $fh "%s : %d : %d : %s : %.3f : %.1f%% : %s\n",
            $seq_id,
            $site->start,
            $site->end,
            $site->strand == -1 ? '-' : '+',
            $site->score,
            $site->rel_score * 100,
            $site->seq->seq();
    }
}

sub write_tfbs_details_text_from_data
{
    my ($tf, $seq_id_display_ids, $data_file, $text_file, %job_args) = @_;

    my $logger = $job_args{-logger};
    my $tf_db  = $job_args{-tf_db};

    my $tf_id   = $tf->ID;
    my $tf_name = $tf->name;

    unless (open(DFH, $data_file)) {
        $logger->error("Error opening TFBS detail data file $data_file - $!");
        return;
    }

    unless (open(OFH, ">$text_file")) {
        $logger->error("Error creating TFBS detail text file $text_file - $!");
        return;
    }

    my $total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $total_ic = 'NA';
    }

    print  OFH "$tf_name\n\n";

    if ($tf_db) {
        print  FH "JASPAR ID:\t$tf_id\n";
    }
    printf OFH "Class:\t%s\n", $tf->class() || 'NA',
    printf OFH "Family:\t%s\n", $tf->tag('family') || 'NA',
    printf OFH "Sysgroup:\t%s\n", $tf->tag('tax_group') || 'NA',
    printf OFH "IC:\t%s\n", $total_ic;

    print OFH "\n\n$tf_name Binding Sites\n\n";

    printf OFH "\n\n%-31s\t%7s\t%7s\t%7s\t%7s\t%7s\t%s\n",
        'Sequence ID', 'Start', 'End', 'Strand', 'Score', '%Score',
        'TFBS Sequence';

    my $last_seq_id = '';
    while (my $line = <DFH>) {
        chomp $line;

        my @cols = split /\s+:\s+/, $line;

        next if $cols[0] eq 'seq_id';  # skip header

        my $seq_id = $cols[0];

        if ($seq_id eq $last_seq_id) {
            printf OFH "%-31s", '';
        } else {
            my $display_id = $seq_id_display_ids->{$seq_id};

            printf OFH "%-31s", $display_id;

            $last_seq_id = $seq_id;
        }

        printf(OFH "\t%7s\t%7s\t%7s\t%7s\t%7s\t%s\n",
            $cols[1],
            $cols[2],
            $cols[3],
            $cols[4],
            $cols[5],
            $cols[6]
        );
    }
}

sub write_tfbs_details_html_from_data
{
    my ($tf, $seq_id_display_ids, $data_file, $html_file, %job_args) = @_;

    my $tf_name = $tf->name();
    my $tf_id   = $tf->ID();

    my $job_id          = $job_args{-job_id};
    my $heading         = $job_args{-heading};
    my $email           = $job_args{-email};
    my $tf_db           = $job_args{-tf_db};
    my $logger          = $job_args{-logger};
    my $rel_results_dir = $job_args{-rel_results_dir};

    open(FH, ">$html_file") || fatal(
        "Could not create TFBS details html file $html_file", %job_args
    );

    $logger->info("Writing '$tf_id - $tf_name' TFBS details to $html_file");
    
    my $title = "oPOSSUM $heading";
    my $section = sprintf("%s Binding Site Details", $tf_name);

    my $fname = $tf_id;
    $fname =~ s/\//_/g;

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

        tf_db                   => $tf_db,
        tf                      => $tf,
        seq_id_display_ids      => $seq_id_display_ids,
        rel_results_dir         => $rel_results_dir,
        data_file               => $data_file,
        tfbs_details_file       => "$fname.txt",

        formatf                 => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'NA'
                                   },

        var_template            => "tfbs_details_seq_ssa.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}

#
# Deprecated routine. This routine is replaced by the
# write_tfbs_details_text_from_data routine.
# DJA 2012/11/07
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
# Deprecated routine. This routine is replaced by the
# write_tfbs_details_html_from_data routine.
# DJA 2012/11/07
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
