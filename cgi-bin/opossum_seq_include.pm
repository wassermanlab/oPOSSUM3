use Carp;

sub matrix_set_compute_gc_content
{
    my $matrix_set = shift;

    my $iter = $matrix_set->Iterator;
    while (my $matrix = $iter->next) {
        my $gc_content = matrix_compute_gc_content($matrix);

        $matrix->tag('gc_content', $gc_content);
    }
}

sub matrix_compute_gc_content
{
    my $pfm = shift;

    unless ($pfm->isa("TFBS::Matrix::PFM")) {
        carp "Cannot compute GC content for non-PFM matrix\n"; 
    }

    my $matrix = $pfm->matrix();

    my $gc_count = 0;
    my $total_count = 0;
    my $row_num = 0;
    foreach my $row (@$matrix) {
        $row_num++;
        foreach my $val (@$row) {
            if ($row_num == 2 || $row_num == 3) {
                $gc_count += $val;
            }

            $total_count += $val;
        }
    }

    my $gc_content = $gc_count / $total_count;

    return $gc_content;
}

sub compute_tf_peak_distances
{
    my ($tf_set, $seq_id_seqs, $tf_seq_sites, $seq_peak_max_pos) = @_;

    my $tf_ids = $tf_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

    my %tfbs_distribution;
    foreach my $tf_id (@$tf_ids) {
        my $num_sites = 0;
        my $num_seqs  = 0;
        my $min_dist  = 9999999;
        my $max_dist  = 0;
        my %dist_count;

        foreach my $seq_id (@seq_ids) {
            my $sites = $tf_seq_sites->{$tf_id}->{$seq_id};

            if ($sites) {
                $num_seqs++;

                my $seq = $seq_id_seqs->{$seq_id};

                my $seq_length = $seq->length;

                my $peak_max_loc = 0;
                if ($seq_peak_max_pos) {
                    $peak_max_loc = $seq_peak_max_pos->{$seq_id};
                } else {
                    # Assume centre of peak sequence is the peak max location
                    $peak_max_loc = int ($seq->length / 2);
                }

                foreach my $site (@$sites) {
                    $num_sites++;

                    my $site_start = $site->start;
                    my $site_end   = $site->end;

                    my $dist = 0;
                    if ($site_end < $peak_max_loc) {
                        $dist = $peak_max_loc - $site_end;
                    } elsif ($site_start > $peak_max_loc) {
                        $dist = $site_start - $peak_max_loc;
                    }

                    if ($dist < $min_dist) {
                        $min_dist = $dist;
                    }

                    if ($dist > $max_dist) {
                        $max_dist = $dist;
                    }

                    $dist_count{$dist}++;
                }
            }
        }

        $tfbs_distribution{$tf_id} = {
            -num_sites  => $num_sites,
            -num_seqs   => $num_seqs,
            -min_dist   => $min_dist,
            -max_dist   => $max_dist,
            -counts     => $dist_count
        }
    }

    return %tfbs_distribution ? \%tfbs_distribution : undef;
}

1;
