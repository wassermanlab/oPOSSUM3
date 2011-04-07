#!/usr/local/bin/perl -w

=head1 NAME

sge_run_opossum_seq.pl

=head1 SYNOPSIS

  sge_run_opossum_seq.pl -s seq_file
                  -b back_seq_file
                  -d results_dir
                  [-db tf_database]
                  [-co collection]
                  [-tax tax_groups]
                  [-ic min_ic]
                  [-ids tf_ids]
                  [-tf tf_file]
                  [-th threshold]
                  [-n num_results | -zcutoff cutoff -fcutoff cutoff]
                  [-sr sort_by]
                  [-m email]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -s seq_file      = Input test sequences file
   -b seq_file      = Input background sequences file
   -d results_dir   = Name of directory used for input sequence and TF
                      files and output results files
   -db tf_database  = Specify which TF database to use
                      (default = JASPAR_2010)
   -co collection   = Specify a JASPAR collection
   -tax tax_groups  = Specify a comma separated string of tax groups
   -ic min_ic       = Specify the minimum IC
   -ids tf_ids      = Specify a comma separated string of TF IDs
   -tf tf_file      = Input TFBS matrix file. If specified overrides above
                      TFBS parameters.
   -th threshold    = Minimum relative TFBS position weight matrix
                      (PWM) score to report in the analysis. The
                      thresold may be spefifies as a percentage, i.e.
                      '80%' or a decimal, i.e. 0.8.
                      Default = 80%
                      Min. = 75%
   -n num_results   = Number of results to display
   -zcutoff cutoff  = Z-score cutoff of results to display
   -fcutoff cutoff  = Fisher p-value cutoff of results to display
   -sr sort_by      = Sort results by this value ('zscore', 'fisher')
   -m email         = E-mail address of user

=head1 DESCRIPTION

Run the opossum_seq.pl script on the SGE cluster. Copy neccesary input files
to cluster, run opossum_seq.pl command and copy back results files.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use Getopt::Long;
use Pod::Usage;
use File::Temp qw/ tempfile tempdir /;
use Log::Log4perl qw(get_logger :levels);

use oPossumWebOpt.pm;

use OPOSSUM_LIB_PATH;

use constant DEBUG          => 0;

my $t_seq_file;
my $bg_seq_file;
my $results_dir;
my $tf_db;
my $tf_collection;
my $tax_groups_str;
my $min_ic;
my $tf_ids_str;
my $tf_file;
my $threshold;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $email;
GetOptions(
    's=s'       => \$t_seq_file,
    'b=s'       => \$bg_seq_file,
    'd=s'       => \$results_dir,
    'db=s'      => \$tf_db,
    'co=s'      => \$tf_collection,
    'tax=s'     => \$tax_groups_str,
    'ic=s'      => \$min_ic,
    'ids=s'     => \$tf_ids_str,
    'tf=s'      => \$tf_file,
    'th=s'      => \$threshold,
    'n=s'       => \$num_results,   # integer or string 'All'
    'zcutoff=f' => \$zscore_cutoff,
    'fcutoff=f' => \$fisher_cutoff,
    'sr=s'      => \$sort_by,
    'm=s'       => \$email
);

die "No results directory specified\n" if !$results_dir;

# Create relative results dir name from abs results dir
my $results_subdir = $results_dir;

# Remove absolute path
$results_subdir =~ s/.*\///; 

# Add relative path
my $rel_results_dir = REL_HTDOCS_RESULTS_PATH . "/$results_subdir";

#
# Initialize logging
#
my $log_file = "$results_dir/sge_run_opossum_seq.log";

my $logger = get_logger();
if (DEBUG) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}

#my $layout = Log::Log4perl::Layout::PatternLayout->new("%M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("[%d] %p\t%m%n");

my $appender = Log::Log4perl::Appender->new(
    "Log::Dispatch::File",
    filename    => $log_file,
    mode        => "append"
);

$appender->layout($layout);
$logger->add_appender($appender);

$logger->info("SGE Run Initialized");

$logger->logdie("No sequences file specified") unless $t_seq_file;
$logger->logdie("Sequences file $t_seq_file does not exist")
    unless -e $t_seq_file;

$logger->logdie("No background sequences file specified") unless $bg_seq_file;
$logger->logdie("Background sequences file $bg_seq_file does not exist")
    unless -e $bg_seq_file;
 
system("scp -p $t_seq_file 'watson.cmmt.ubc.ca:~/$results_subdir/$'");

#my $tf_file = "$results_dir/matrices.txt";

my $matrix_set;
if ($tf_file) {
    unless (-e $tf_file) {
        $logger->logdie(
            "Specified input TF profile file $tf_file does not exist"
        );
    }

    $matrix_set = read_PFMs($tf_file);

    unless ($matrix_set && $matrix_set->size > 0) {
        $logger->logdie("Error reading TF profiles from $tf_file");
    }
} else {
    #$tf_file = undef;

    unless ($tf_db) {
        $logger->logdie("No JASPAR database specified");
    }

    unless ($tf_collection) {
        $logger->logdie("No JASPAR collection specified");
    }

    #
    # Connect to JASPAR database
    #
    my $jdb = TFBS::DB::JASPAR5->connect(
        "dbi:mysql:" . $tf_db . ":" . JASPAR_DB_HOST,
        JASPAR_DB_USER,
        JASPAR_DB_PASS
    );

    $logger->logdie("Could not connect to JASPAR database $tf_db") if !$jdb;

    my %get_matrix_args = (
        -collection => $tf_collection,
        -matrixtype => 'PFM'
    );

    if (@tf_ids) {
        # This takes precendence over tax groups and min IC
        $get_matrix_args{-ID} = \@tf_ids;
    } else {
        if (@tax_groups) {
            $get_matrix_args{-tax_group} = \@tax_groups;
        }

        if (defined $min_ic) {
            $get_matrix_args{-min_ic} = $min_ic;
        }
    }

    $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

    unless ($matrix_set && $matrix_set->size > 0) {
        $logger->logdie("Error fetching TF profiles from JASPAR DB");
    }
}

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

my ($t_seq_tf_siteset, $t_tf_seqs)
    = tf_set_search_seqs($tf_set, $seqs, $threshold);

$logger->logdie("No TFBSs found in test sequences") if !$t_seq_tf_siteset;

my ($c_seq_tf_siteset, $c_tf_seqs)
    = tf_set_search_seqs($tf_set, $back_seqs, $threshold);

$logger->logdie("No TFBSs found in control sequences") if !$c_seq_tf_siteset;

my $t_counts = compute_site_counts($tf_set, $t_seq_tf_siteset);
my $c_counts = compute_site_counts($tf_set, $c_seq_tf_siteset);

#$t_counts->conserved_region_length_set($t_cr_len_set);
#$c_counts->conserved_region_length_set($c_cr_len_set);

#my $fresults = fisher_analysis($c_counts, $t_counts, \%tfs,
#                               scalar @$seqs, scalar @$back_seqs);
my $fisher = OPOSSUM::Analysis::Fisher->new();
my $fresults = $fisher->calculate_Fisher_probability($c_counts, $t_counts);

#my $zresults = zscore_analysis($c_counts, $t_counts, \%tfs,
#                               $t_seq_len, $c_seq_len);
my $zscore = OPOSSUM::Analysis::Zscore->new();
my $zresults = $zscore->calculate_Zscore(
    $c_counts, $t_counts, $c_seq_len, $t_seq_len, $tf_set
);

# Sort results by z-score
#$zresults->sort_by('z_score', 1);
#
#my $formatted_results = format_results($fresults, $zresults, $tf_set);
#

#
# Use new OPOSSUM::Analysis::CombinedResultSet to combine Fisher and
# Z-score result sets.
#
my $cresults = OPOSSUM::Analysis::CombinedResultSet->new(
    -fisher_result_set  => $fresults,
    -zscore_result_set  => $zresults
);

#
# Get results as a list
#
my %result_params;
$result_params{-num_results} = $num_results if defined $num_results;
$result_params{-zscore_cutoff} = $zscore_cutoff if defined $zscore_cutoff;
$result_params{-fisher_cutoff} = $fisher_cutoff if defined $fisher_cutoff;
if (defined $sort_by) {
    if ($sort_by =~ /^fisher/) {
        $sort_by = 'fisher_p_value';
    } elsif ($sort_by =~ /^z_score/ || $sort_by =~ /^z-score/) {
        $sort_by = 'zscore';
    }

    $result_params{-sort_by} = $sort_by;

    if ($sort_by eq 'zscore') {
        # Sort z-score from highest to lowest
        $result_params{-reverse} = 1;
    }
}

my $cresult_list = $cresults->get_list(%result_params);

write_results_text($tf_set, $cresult_list);

write_results_html($tf_set, $cresult_list);

write_tfbs_details_text($tf_set, $t_seq_tf_siteset, $t_tf_seqs);

write_tfbs_details_html($tf_set, $c_seq_tf_siteset, $c_tf_seqs);

$logger->info("Finished analysis");

send_email($email) if $email;

exit;

sub read_seqs
{
    my $file = shift;

    my $seqIO = Bio::SeqIO->new(-file => $file, -format => 'fasta');

    my @seqs;
    while (my $seq = $seqIO->next_seq()) {
        push @seqs, $seq;
    }
    $seqIO->close();

    return @seqs ? \@seqs : undef;
}

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
    my ($tf_set, $seqs, $threshold) = @_;

    my $tf_ids = $tf_set->ids();

    my %seq_tf_siteset;
    my %tf_seqs;

	foreach my $tf_id (@$tf_ids) {
	    my $pfm = $tf_set->get_matrix($tf_id);
	    my $pwm = $pfm->to_PWM();

        foreach my $seq (@$seqs) {
            my $seq_id;

            my $display_id = $seq->display_id();
            if ($display_id =~ /^(chr\w+):(\d+)-(\d+)/) {
                $seq_id = "$1:$2-$3";
            } else {
                $seq_id = $display_id;
            }

            my $siteset = $pwm->search_seq(
                -seqobj     => $seq,
                -threshold  => $threshold
            );

            my $filtered_siteset;
            if ($siteset && $siteset->size > 0) {
                $filtered_siteset = filter_overlapping_sites($siteset);

                if ($filtered_siteset && $filtered_siteset->size > 0) {
                    # Only set if seqs had sites for this TF
                    push @{$tf_seqs{$tf_id}}, $seq_id;
                }
            }

            # Note: this will be undef if there were no sites for this seq/TF
            $seq_tf_siteset{$seq_id}->{$tf_id} = $filtered_siteset;
        }
    }

    my $retval1 = %seq_tf_siteset ? \%seq_tf_siteset : undef;
    my $retval2 = %tf_seqs ? \%tf_seqs : undef;

    return ($retval1, $retval2);
}

sub compute_site_counts
{
    my ($tf_set, $seq_tf_siteset) = @_;

    my $tf_ids  = $tf_set->ids();
    my @seq_ids = keys %$seq_tf_siteset;

    my $counts = OPOSSUM::Analysis::Counts->new(
        -gene_ids       => \@seq_ids,
        -tf_ids         => $tf_ids
    );

    foreach my $seq_id (@seq_ids) {
        foreach my $tf_id (@$tf_ids) {
            my $siteset = $seq_tf_siteset->{$seq_id}->{$tf_id};

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
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
    my ($tf_set, $results) = @_;

    my $filename = "$results_dir/" . RESULTS_TEXT_FILENAME;

    open(FH, ">$filename")
        || $logger->logdie("Could not create analysis results file $filename");

    $logger->info("Writing analysis results to $filename\n");

    #
    # Multi-line (tab-delimited) header format
    #
    #printf FH "TF\tTF\tTF\tIC\tBackgrd\tBackgrd\tTarget\tTarget\tBackgrd\tBackgrd\tTarget\tTarget\tZ-score\tFisher\n";
    #printf FH "\tClass\tS.group\t\tgene\tgene\tgene\tgene\tTFBS\tTFBS\tTFBS\tTFBS\t\tscore\n";
    #printf FH "\t\t\t\thits\tno-hits\thits\tno-hits\thits\trate\thits\trate\t\t\n";
    #
    #
    # Single line (tab-delimited) header format
    # NOTE: rearranged columns
    #
    if ($tf_db) {
        printf FH "TF Name\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tBackground TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

        foreach my $result (@$results) {
            my $tf = $tf_set->get_tf($result->id());

            printf FH 
                "%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                $tf->name(),
                $tf->ID(),
                $tf->class() || 'N/A',
                $tf->tag('family') || 'N/A',
                $tf->tag('tax_group') || 'N/A',
                $tf->to_ICM->total_ic(),
                $result->t_gene_hits() || 0,
                $result->t_gene_no_hits() || 0,
                $result->bg_gene_hits() || 0,
                $result->bg_gene_no_hits() || 0,
                $result->t_tfbs_hits() || 0,
                $result->bg_tfbs_hits() || 0,
                defined $result->t_tfbs_rate()
                    ? sprintf("%.3f", $result->t_tfbs_rate()) : 'N/A',
                defined $result->bg_tfbs_rate()
                    ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'N/A',
                defined $result->zscore()
                    ? sprintf("%.3f", $result->zscore()) : 'N/A',
                defined $result->fisher_p_value()
                    ? sprintf("%.3g", $result->fisher_p_value()) : 'N/A';
        }
    } else {
        printf FH "TF Name\tClass\tFamily\tTax Group\tIC\tTarget gene hits\tTarget gene non-hits\tBackground gene hits\tBackground gene non-hits\tTarget TFBS hits\tBackground TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

        foreach my $result (@$results) {
            my $tf = $tf_set->get_tf($result->id());

            printf FH 
                "%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                $tf->name(),
                $tf->class() || 'N/A',
                $tf->tag('family') || 'N/A',
                $tf->tag('tax_group') || 'N/A',
                $tf->to_ICM->total_ic(),
                $result->t_gene_hits() || 0,
                $result->t_gene_no_hits() || 0,
                $result->bg_gene_hits() || 0,
                $result->bg_gene_no_hits() || 0,
                $result->t_tfbs_hits() || 0,
                $result->bg_tfbs_hits() || 0,
                defined $result->t_tfbs_rate()
                    ? sprintf("%.3f", $result->t_tfbs_rate()) : 'N/A',
                defined $result->bg_tfbs_rate()
                    ? sprintf("%.3f", $result->bg_tfbs_rate()) : 'N/A',
                defined $result->zscore()
                    ? sprintf("%.3f", $result->zscore()) : 'N/A',
                defined $result->fisher_p_value()
                    ? sprintf("%.3g", $result->fisher_p_value()) : 'N/A';
        }
    }
    close(FH);

    return $filename;
}


#
# Ouput combined z-score/Fisher results as HTML
#
sub write_results_html
{
    my ($tf_set, $results) = @_;

    my $tf_ids = $tf_set->ids();

    my $warn_zero_bg_gene_hits;
    foreach my $result (@$results) {
        if (!$result->bg_gene_hits()) {
            $warn_zero_bg_gene_hits = 1;
            last;
        }
    }

    my $result_type;
    if (defined $zscore_cutoff || defined $fisher_cutoff) {
        $result_type = 'significant_hits';
    } else {
        $result_type = 'top_x_results';
    }

    my $vars = {
        abs_htdocs_path         => ABS_HTDOCS_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_htdocs_path         => REL_HTDOCS_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        version                 => VERSION,
        devel_version           => DEVEL_VERSION,
        jaspar_url              => JASPAR_URL,
        heading                 => HEADING,
        result_retain_days      => REMOVE_RESULTFILES_OLDER_THAN,
        title                   => "oPOSSUM Results",

        #
        # XXX These are the local file names on the server and not particularly
        # informative to the user.
        #
        #seq_file            => $t_seq_file,
        #back_file           => $bg_seq_file,

        formatf             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*f", $dec, $f)
                                        : 'N/A'
                               },

        formatg             => sub {
                                    my $dec = shift;
                                    my $f = shift;
                                    return ($f || $f eq '0')
                                        ? sprintf("%.*g", $dec, $f)
                                        : 'N/A'
                               },

        num_t_seqs          => scalar @$seqs,
        num_c_seqs          => scalar @$back_seqs,
        tf_file             => $tf_file,
        tf_db               => $tf_db,
        tf_collection       => $tf_collection,
        tf_ids              => $tf_ids,
        tax_groups          => \@tax_groups,
        threshold           => $threshold,
        tf_set              => $tf_set,
        results             => $results,
        rel_results_dir     => $rel_results_dir,
        result_type         => $result_type,
        num_results         => $num_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        sort_by             => $sort_by,
        warn_zero_bg_gene_hits  => $warn_zero_bg_gene_hits,
        results_file        => RESULTS_TEXT_FILENAME,
        var_template        => "results_seq.html"
    };

    my $output = process_template('master.html', $vars);

    my $html_filename = "$results_dir/" . RESULTS_HTDOCS_FILENAME;

    open(OUT, ">$html_filename")
        || $logger->logdie("Could not create HTML results file $html_filename");

    print OUT $output;

    close(OUT);

    $logger->info("Wrote HTML formatted results to $html_filename");

    return $html_filename;
}

#
# Write the details of the putative TFBSs for each TF/gene. Create a
# text file for each TF.
#
sub write_tfbs_details_text
{
    my ($tf_set, $seq_tf_siteset, $tf_seqs) = @_;

    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $tf = $tf_set->get_tf($tf_id);

        my $tf_name = $tf->name();

        my $filename = "$results_dir/$tf_id.txt";

        open(FH, ">$filename") || $logger->logdie(
            "Could not create TFBS details text file $filename"
        );

        $logger->info("Writing '$tf_name' TFBS details to $filename");

        print  FH "Name     :\t$tf_name\n";
        if ($tf_db) {
            print  FH "JASPAR ID:\t$tf_id\n";
        }
        printf FH "Class    :\t%s\n", $tf->class() || 'N/A',
        printf FH "Family   :\t%s\n", $tf->tag('family') || 'N/A',
        printf FH "Sysgroup :\t%s\n", $tf->tag('tax_group') || 'N/A',
        printf FH "IC       :\t%.3f\n", $tf->to_ICM->total_ic();

        print FH "\n\n";

        if ($tf_seqs->{$tf_id}) {
            #printf FH "\n\n%s\t%s\t%s\t%s\t%s\t%s\n",
            printf FH "\n\n%-31s\t%7s\t%7s\t%7s\t%7s\t%7s\t%s\n",
                'Sequence ID', 'Start', 'End', 'Strand', 'Score', '%Score', 'Bound sequence';

            my @seq_ids = sort @{$tf_seqs->{$tf_id}};
            foreach my $seq_id (@seq_ids) {
                my $siteset = $seq_tf_siteset->{$seq_id}->{$tf_id};

                next if !$siteset || $siteset->size == 0;

                my $first = 1;
                printf FH "%-31s", $seq_id;

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
        } else {
            print FH "No sites found for this TF\n";
        }

        close(FH);
    }
}

#
# Write the details of the putative TFBSs for each TF/gene. Create an
# HTML file for each TF.
#
sub write_tfbs_details_html
{
    my ($tf_set, $seq_tf_siteset, $tf_seqs) = @_;

    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $tf = $tf_set->get_tf($tf_id);

        my $tf_name = $tf->name();

        my %seq_sites;
        my @seq_ids;
        if ($tf_seqs->{$tf_id}) {
            @seq_ids = sort @{$tf_seqs->{$tf_id}};
            foreach my $seq_id (@seq_ids) {
                my $siteset = $seq_tf_siteset->{$seq_id}->{$tf_id};

                next if !$siteset || $siteset->size == 0;

                my @site_list;
                my $iter = $siteset->Iterator(-sort_by => 'start');
                while (my $site = $iter->next()) {
                    push @site_list, $site;
                }

                $seq_sites{$seq_id} = \@site_list;
            }
        }

        my $html_filename = "$results_dir/$tf_id.html";

        open(FH, ">$html_filename") || $logger->logdie(
            "Could not create TFBS details HTML file $html_filename"
        );

        $logger->info("Writing '$tf_name' TFBS details to $html_filename");

        my $vars = {
            abs_htdocs_path         => ABS_HTDOCS_PATH,
            abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
            rel_htdocs_path         => REL_HTDOCS_PATH,
            rel_cgi_bin_path        => REL_CGI_BIN_PATH,
            version                 => VERSION,
            devel_version           => DEVEL_VERSION,
            jaspar_url              => JASPAR_URL,
            heading                 => HEADING,
            title                   => "oPOSSUM TFBS Detail",

            formatf                 => sub {
                                        my $dec = shift;
                                        my $f = shift;
                                        return ($f || $f eq '0')
                                            ? sprintf("%.*f", $dec, $f)
                                            : 'N/A'
                                       },

            tf_db                   => $tf_db,
            tf                      => $tf,
            seq_ids                 => \@seq_ids,
            seq_sites               => \%seq_sites,
            rel_results_dir         => $rel_results_dir,
            tfbs_details_file       => "$tf_id.txt",
            var_template            => "tfbs_details_seq.html"
        };

        my $output = process_template('master.html', $vars);

        print FH $output;

        close(FH);
    }
}

sub read_PFMs
{
    my ($file) = @_;

    open(FH, $file) || $logger->logdie("Could not open PFM file $file - $!");

    my $matrix_set = TFBS::MatrixSet->new();

    my $name            = '';
    my $matrix_string   = '';
    my $line_count      = 0;
    my $matrix_count    = 0;
    while (my $line = <FH>) {
        chomp $line;

        next if !$line;

        if ($line =~ /^>\s*(\S+)/) {
            $name = $1;
        } else {
            if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            } elsif ($line =~ /^\s*\d+/) {
                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            } else {
                next;
            }
            $line_count++;

            if ($line_count == 4) {
                my $id = sprintf "matrix%d", $matrix_count + 1;

                unless ($name) {
                    $name = $id;
                }

                my $pfm = TFBS::Matrix::PFM->new(
                    -ID           => $id,
                    -name         => $name,
                    -matrixstring => $matrix_string
                );

                $matrix_set->add_matrix($pfm);

                $matrix_string = '';
                $name = '';
                $line_count = 0;
                $matrix_count++;
            }
        }
    }
    close(FH);

    return $matrix_set;
}

#
# This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SitePairSet where each site pair in the set corresponds to the
# same transcription factor and filter overlapping site pairs such that only
# the highest scoring site pair of any mutually overlapping site pairs is kept.
# In the event that site pairs score equally, the first site pair is kept, i.e.
# bias is towards the site pair with the lowest starting position.
#
sub filter_overlapping_sites
{
    my ($siteset) = @_;

    return if !defined $siteset || $siteset->size == 0;

    my $filtered_set = TFBS::SiteSet->new();

    my $iter = $siteset->Iterator(-sort_by => 'start');
    my $prev_site = $iter->next;
    if ($prev_site) {
        while (my $site = $iter->next) {
            if ($site->overlaps($prev_site)) {
                #
                # Bias is toward the site pair with the lower start
                # site (i.e. if the scores are equal).
                # 
                if ($site->score > $prev_site->score) {
                    $prev_site = $site;
                }
            } else {
                $filtered_set->add_site($prev_site);
                $prev_site = $site;
            }
        }
        $filtered_set->add_site($prev_site);
    }

    return $filtered_set;
}

sub process_template
{
    my ($template_name, $vars) = @_;

    my $config = {
        ABSOLUTE        => 1,
        INCLUDE_PATH    => ABS_HTDOCS_TEMPLATE_PATH . "/",  # or list ref
        INTERPOLATE     => 1,   # expand "$var" in plain text
        POST_CHOMP      => 1,   # cleanup whitespace
        #PRE_PROCESS     => 'header',   # prefix each template
        EVAL_PERL => 1    # evaluate Perl code blocks
    };

    my $string   = '';
    my $template = Template->new($config);
    my $input    = ABS_HTDOCS_TEMPLATE_PATH . "/$template_name";
    $template->process($input, $vars, \$string)
        || $logger->logdie($template->error());

    return $string;
}

sub send_email
{
    my ($email) = @_;

    return if !$email;

    my $results_url = sprintf "%s%s/%s",
        WEB_SERVER_URL,
        "$rel_results_dir",
        RESULTS_HTDOCS_FILENAME;

    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg = "Your oPOSSUM sequence analysis results are now available.\n\n";
    $msg .= "They can be viewed at $results_url\n\n";
    $msg .= "Thank-you,\n";
    $msg .= "the oPOSSUM development team\n";
    $msg .= ADMIN_EMAIL . "\n";

    if (!open(SM, "|" . $cmd)) {
        $logger->error("Could not open sendmail - $!");
        return;
    }

    printf SM "To: %s\n", $email;
    print SM "Subject: oPOSSUM sequence-based analysis results\n\n";
    print SM "$msg" ;

    close(SM);
}

sub error
{
    my ($error) = @_;

    $error = 'Unknown' if !$error;

    my $admin_msg = "oPOSSUM sequence-based analysis failed\n";
    $admin_msg .= "\nError: $error\n";
    $admin_msg .= "\nResults directory: $results_dir\n" if $results_dir;
    $admin_msg .= "\nUser e-mail: $email\n" if $email;

    if ($logger) {
        $logger->logdie("$error");
    }
}
