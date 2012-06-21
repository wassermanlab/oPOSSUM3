=head1 NAME

 OPOSSUM::Include::SeqSSAInclude.pm

=head1 SYNOPSIS

=head1 DESCRIPTION

  Contains all options and routines that are common to all the sequence-based
  SSA scripts and web modules.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use OPOSSUM::Include::SeqInclude;
use OPOSSUM::Include::TCAInclude;

use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::Values;

use constant BG_COLOR_CLASS => 'bgc_seq_tca';

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
sub tf_cluster_set_search_seqs
{
    my ($tf_cluster_set, $tf_set, $seq_id_seqs, $threshold) = @_;

    my $cluster_ids = $tf_cluster_set->ids();

    my @seq_ids = keys %$seq_id_seqs;
                
    #
    # If threshold is specified as a decimal, convert it to a
    # percentage, otherwise the TFBS::Matrix::PWM::search_seq
    # method treats the number as an absolute matrix score which
    # is not what we intended. DJA 2012/06/07
    #
    unless ($threshold =~ /(.+)%$/ || $threshold > 1) {
        $threshold *= 100;
        $threshold .= '%';
    }

    my %cluster_seq_sites;

    my $count = 1;
    foreach my $seq_id (@seq_ids)
    {
        my $seq = $seq_id_seqs->{$seq_id};
        
        foreach my $cl_id (@$cluster_ids)
        {
            my $cluster = $tf_cluster_set->get_tf_cluster($cl_id);
            my $tf_ids = $cluster->tf_ids();
            my $cluster_siteset = TFBS::SiteSet->new();
            
            foreach my $tf_id (@$tf_ids)
            {
                my $matrix = $tf_set->get_matrix($tf_id);    
                next if !$matrix;
        
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
                
                # No real need to filter - to be merged later anyways
                #my $filtered_siteset;
                #if ($siteset && $siteset->size > 0) {
                #    $filtered_siteset = filter_overlapping_sites($siteset);    
                #}
                
                $cluster_siteset->add_siteset($siteset);
            }
            
            # now, merge the sites for this cluster
            # listref of OPOSSUM::ConservedTFBS objects returned
            my $merged_sites = merge_cluster_tfbs_siteset($cluster_siteset, $cl_id, $seq_id);
            
            # Only set if seqs had sites for this TF
            #push @{$tf_seqs{$tf_id}}, $seq_id;
            if ($merged_sites and scalar(@$merged_sites) > 0) {
                $cluster_seq_sites{$cl_id}->{$seq_id} = $merged_sites;
            }
        }
        $count++;
    }

    my $retval1 = %cluster_seq_sites ? \%cluster_seq_sites : undef;

    return $retval1;
}

sub compute_cluster_counts
{
    my ($cluster_set, $seq_ids, $cluster_seq_sites) = @_;

    my $cl_ids  = $cluster_set->ids();

    my $counts = OPOSSUM::Analysis::Cluster::Counts->new(
        -gene_ids       => $seq_ids,
        -cluster_ids    => $cl_ids
    );

    foreach my $seq_id (@$seq_ids) {
        foreach my $cl_id (@$cl_ids) {
            my $sites = $cluster_seq_sites->{$cl_id}->{$seq_id};

            if ($sites) {
                # Note set size could be 0.
                $counts->gene_cluster_count($seq_id, $cl_id, scalar @$sites);
                my $length = 0;
                #my $iter = $siteset->Iterator(-sort_by => 'start');
                #while (my $site = $it->next) {
                foreach my $site (@$sites) {
                    $length += length($site->seq); # site => OPOSSUM::ConservedTFBS
                }
                $counts->gene_cluster_length($seq_id, $cl_id, $length);
            } else {
                $counts->gene_cluster_count($seq_id, $cl_id, 0);
                $counts->gene_cluster_length($seq_id, $cl_id, 0);
            }
        }
    }

    return $counts;
}

=head2 compute_cluster_peak_distances

 Title   : compute_cluster_peak_distances

 Function: Calculate the distance of each TFBS cluster from the peak max height
           location and return an array of distances. For the purpose of
           distance computation, use the center of the TFBS clusters. If
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

sub compute_cluster_peak_distances
{
    my ($cluster_set, $seq_id_seqs, $seq_id_display_ids, $cluster_seq_sites, $seq_peak_max_pos) = @_;

    my $cluster_ids = $cluster_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

    #my %distances;
    my $distances = OPOSSUM::Analysis::Cluster::Values->new(
        -seq_ids => \@seq_ids,
        -cluster_ids => $cluster_ids,
    );
    
    foreach my $cl_id (@$cluster_ids) {

        foreach my $seq_id (@seq_ids) {
            my $sites = $cluster_seq_sites->{$cl_id}->{$seq_id};

            if ($sites) {
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
                
                #my $iter = $siteset->Iterator(-sort_by => 'start');
                #while (my $site = $iter->next) {
                foreach my $site (@$sites) {
                    #print $site->pattern->ID . "\t" . $site->start . "\t" . $site->seq->length . "\n";
                    #my $site_start  = $site->{start};
                    #my $site_length = $site->{seq}->length;

                    my $site_loc = $site->start + length($site->seq) / 2;

                    my $dist = abs($site_loc - $peak_max_loc);
                    my $dist_n = $dist / $seq->length;
                    
                    $distances->add_seq_cluster_value($seq_id, $cl_id, $dist_n);
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
    my ($filename, $results, $tf_cluster_set, %job_args) = @_;

    return unless $results && $results->[0];

    my $text = "TFBS Cluster ID\tClass\tFamily";
    $text .= "\tTarget seq hits\tTarget seq non-hits";
    $text .= "\tBackground seq hits\tBackground seq non-hits";
    $text .= "\tTarget cluster hits\tTarget cluster nucleotide rate";
    $text .= "\tBackground cluster hits\tBackground cluster nucleotide rate";
    $text .= "\tZ-score\tFisher score\tKS score\n";

    foreach my $result (@$results) {
        my $cl = $tf_cluster_set->get_tf_cluster($result->id());

        $text .= sprintf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
            'c' . $cl->id(),
            $cl->class() || 'NA',
            $cl->family() || 'NA',
            $result->t_gene_hits() || 0,
            $result->t_gene_no_hits() || 0,
            $result->bg_gene_hits() || 0,
            $result->bg_gene_no_hits() || 0,
            $result->t_cluster_hits() || 0,
            defined $result->t_cluster_rate()
                ? sprintf("%.3f", $result->t_cluster_rate()) : 'NA',
            $result->bg_cluster_hits() || 0,
            defined $result->bg_cluster_rate()
                ? sprintf("%.3f", $result->bg_cluster_rate()) : 'NA',
            defined $result->zscore()
                ? sprintf("%.3f", $result->zscore()) : 'NA',
            defined $result->fisher_p_value()
                ? sprintf("%.3f", $result->fisher_p_value()) : 'NA',
            defined $result->ks_p_value()
                ? sprintf("%.3f", $result->ks_p_value()) : 'NA';
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
# Write the details of the putative TFBS clusters for each TFCluster/seq.
# Create a text file for each TFCluster.
#
sub write_tfbs_cluster_details_text
{
    my ($filename, $tfcl, $cl_db,
        $seq_ids, $seq_id_display_ids, $seq_sites,
        %job_args) = @_;
    
    
    unless (open(FH, ">$filename")) {
        fatal("Unable to create TFBS cluster details file $filename", %job_args);
        return;
    }

    printf FH "C%s\n\n", $tfcl->id();
    printf FH "Class:\t%s\n", $tfcl->class || 'NA';
    printf FH "Family:\t%s\n", $tfcl->family || 'NA';
    printf FH "\n\nConserved C%s Binding Sites:\n\n", $tfcl->id();
    printf FH "\n\n%-31s\t%7s\t%7s\t%7s\t%7s\t%7s\t%s\n",
        'Sequence ID', 'Start', 'End', 'Strand', 'Score', '%Score', 'TFBS Sequence';

    foreach my $seq_id (@$seq_ids)
    {
        my $sites = $seq_sites->{$seq_id};

        my $display_id = $seq_id_display_ids->{$seq_id};

        next if !$sites || scalar @$sites == 0;

        my $first = 1;
        printf FH "%-31s", $display_id;

        #my $iter = $siteset->Iterator(-sort_by => 'start');
        #while (my $site = $iter->next()) {
        foreach my $site (@$sites)
        {
            # Do not reprint sequence ID for subsequent sites
            unless ($first) {
                printf FH "%-31s", '';
            }
            $first = 0;

            printf FH "\t%7d\t%7d\t%7d\t%7.3f\t%6.1f%%\t%s\n",
                $site->start(),
                $site->end(),
                $site->strand(),
                $site->score(),
                $site->rel_score() * 100,
                #$site->seq->seq();
                $site->seq;
        }

        printf FH "\n";
    }

    close(FH);
}

#
# Write the details of the putative TFBSs for each TF/seq. Create an
# HTML file for each TF.
#
sub write_tfbs_cluster_details_html
{
    my ($filename, $rel_results_dir,
        $tfcl, $cl_db,
        $seq_ids, $seq_id_display_ids, $seq_sites,
        %job_args) = @_;
    
    my $job_id = $job_args{-job_id};
    my $heading = $job_args{-heading};
    my $email = $job_args{-email};
    my $logger = $job_args{-logger};

    my $tfcl_id   = $tfcl->id();
    
    open(FH, ">$filename") || fatal(
        "Could not create TFBS cluster details HTML file $filename",
        %job_args
    );

    $logger->info("Writing '$tfcl_id' TFBS details to $filename");
    
    my $title = "oPOSSUM $heading";
    my $section = sprintf("C%s Binding Site Details", $tfcl_id);

    my %seq_sites;
    foreach my $seq_id (@$seq_ids) {
        my $sites = $seq_sites->{$seq_id};

        next if !$sites || scalar @$sites == 0;
        $seq_sites{$seq_id} = $sites;
    }

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

        #tf_db                   => $tf_db,
        tf_cluster              => $tfcl,
        seq_ids                 => $seq_ids,
        seq_id_display_ids      => $seq_id_display_ids,
        seq_sites               => \%seq_sites,
        rel_results_dir         => $rel_results_dir,
        tfbs_cluster_details_file       => "c$tfcl_id.txt",
        var_template            => "tfbs_details_seq_tca.html"
    };

    my $output = process_template('master.html', $vars, %job_args);

    print FH $output;

    close(FH);
}

1;
