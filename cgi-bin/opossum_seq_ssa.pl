#!/usr/local/bin/perl -w

=head1 NAME

opossum_seq_ssa.pl

=head1 SYNOPSIS

  opossum_seq_ssa.pl
                  -j job_id
                  -s t_seq_file
                  -b bg_seq_file
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
                 [-usf user_seq_file]
                 [-ubsf user_bg_seq_file]
                 [-utff user_tf_file]
                 [-bss bg_seq_set_name]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -j job_id        = The oPOSSUM job ID. This is also the same as the
                      temporary subdirectory created in the results
                      directory.
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
   -usf user_seq_file       = User supplied sequence upload file name
   -ubsf user_bg_seq_file   = User supplied background sequence upload file
                              name
   -utff user_tf_file       = User supplied TFBS matrix upload file name
   -bss bg_seq_set_name     = Background sequence set name

=head1 DESCRIPTION

Take a set of peak sequences as provided in the input peaks file and 
a set of control peak sequences. Optionally take a subset of transcription
factors (TFs) either specified in an input file, or limited by external
(JASPAR) database name and information content or taxonomic supergroup or
all TFs in the oPOSSUM database. Also optionally specify PWM score threshold.

Count the number of TFBSs for each TF which was found at the given
PWM score threshold for both the test and control set of peaks. Perform
Fisher exact test and z-score analysis and output these results to the
output file. Optionally write details of TFBSs found in test set to detailed
TFBS hits file.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use oPossumWebOpt;
use oPossumSeqWebOpt;
use oPossumSSAWebOpt;

use lib OPOSSUM_LIB_PATH;

use Getopt::Long;
use Pod::Usage;
use File::Temp;
use Carp;
#use CGI::Carp qw(carpout);
use Template;
use File::Temp qw/ tempfile tempdir /;
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

use Bio::SeqIO;

use TFBS::DB::JASPAR5;

use OPOSSUM::TFSet;
use OPOSSUM::ConservedRegionLength;
use OPOSSUM::ConservedRegionLengthSet;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Counts;
use OPOSSUM::Analysis::CombinedResultSet;

use Statistics::Distributions;

use constant DEBUG          => 0;

use constant HEADING        => 'Sequence-based Single Site Analysis';
use constant TITLE          => 'oPOSSUM Sequence-based Single Site Analysis';
use constant BG_COLOR_CLASS => 'bgc_seq_ssa';

my $job_id;
my $t_seq_file;
my $bg_seq_file;
my $results_dir;
my $tf_db;
my $collection;
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
my $user_seq_file;
my $user_bg_seq_file;
my $user_tf_file;
my $bg_seq_set_name;
GetOptions(
    'j=s'       => \$job_id,
    's=s'       => \$t_seq_file,
    'b=s'       => \$bg_seq_file,
    'd=s'       => \$results_dir,
    'db=s'      => \$tf_db,
    'co=s'      => \$collection,
    'tax=s'     => \$tax_groups_str,
    'ic=s'      => \$min_ic,
    'ids=s'     => \$tf_ids_str,
    'tf=s'      => \$tf_file,
    'th=s'      => \$threshold,
    'n=s'       => \$num_results,   # integer or string 'All'
    'zcutoff=f' => \$zscore_cutoff,
    'fcutoff=f' => \$fisher_cutoff,
    'sr=s'      => \$sort_by,
    'm=s'       => \$email,
    'usf=s'     => \$user_seq_file,
    'ubsf=s'    => \$user_bg_seq_file,
    'utff=s'    => \$user_tf_file,
    'bss=s'     => \$bg_seq_set_name
);

die "No results directory specified\n" if !$results_dir;

# Create relative results dir name from abs results dir
my $rel_results_dir = $results_dir;

# Remove absolute path
$rel_results_dir =~ s/.*\///; 

# Add relative path
$rel_results_dir = REL_HTDOCS_RESULTS_PATH . "/$rel_results_dir";

#
# Initialize logging
#
my $log_file = "$results_dir/opossum_seq_ssa.log";

my $logger = get_logger();
if (DEBUG) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}

#my $layout = Log::Log4perl::Layout::PatternLayout->new("%M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("[%d] %L %p\t%m%n");

my $appender = Log::Log4perl::Appender->new(
    "Log::Dispatch::File",
    filename    => $log_file,
    mode        => "append"
);

$appender->layout($layout);
$logger->add_appender($appender);

#open(ERR, ">>$log_file") || die "Could not open log file $log_file\n";

#carpout(\*ERR);

$logger->info("Starting analysis");

fatal("No sequences file specified") unless $t_seq_file;
fatal("Sequences file $t_seq_file does not exist")
    unless -e $t_seq_file;

fatal("No background sequences file specified") unless $bg_seq_file;
fatal("Background sequences file $bg_seq_file does not exist")
    unless -e $bg_seq_file;
 
#$tf_db          = JASPAR_DB_NAME if !$tf_db;
#$collection  = 'CORE' if !$collection;

$threshold      = DFLT_THRESHOLD . "%" if !$threshold;

my @tf_ids      = split /\s*,\s*/, $tf_ids_str if $tf_ids_str;
my @tax_groups  = split /\s*,\s*/, $tax_groups_str if $tax_groups_str;

$logger->info("Reading target sequences from $t_seq_file");
my $t_seqs = read_seqs($t_seq_file)
    || fatal("Could not read sequences");

my @t_seq_ids;
my %t_seq_id_seqs;
my %t_seq_id_display_ids;
my $seq_num = 1;
foreach my $seq (@$t_seqs) {
    my $seq_id = "seq$seq_num";

    push @t_seq_ids, $seq_id;

    my $display_id = $seq->display_id();
    if ($display_id) {
        if ($display_id =~ /^(chr\w+):(\d+)-(\d+)/) {
            $display_id = "$1:$2-$3";
        }
    } else {
        $display_id = $seq_id;
    }

    $t_seq_id_display_ids{$seq_id} = $display_id;
    $t_seq_id_seqs{$seq_id} = $seq;

    $seq_num++;
}

my ($t_seq_gc_content, $t_seq_len) = calc_seqs_total_gc_content($t_seqs);

$logger->debug("Target seq GC content: $t_seq_gc_content");
$logger->debug("Target seq total length: $t_seq_len");

$logger->info("Reading background sequences from $bg_seq_file");
my $bg_seqs = read_seqs($bg_seq_file)
    || fatal("Could not read background sequences");

my @bg_seq_ids;
my %bg_seq_id_seqs;
my %bg_seq_id_display_ids;
$seq_num = 1;
foreach my $seq (@$bg_seqs) {
    my $seq_id = "seq$seq_num";

    push @bg_seq_ids, $seq_id;

    my $display_id = $seq->display_id();
    if ($display_id) {
        if ($display_id =~ /^(chr\w+):(\d+)-(\d+)/) {
            $display_id = "$1:$2-$3";
        }
    } else {
        $display_id = $seq_id;
    }

    $bg_seq_id_display_ids{$seq_id} = $display_id;
    $bg_seq_id_seqs{$seq_id} = $seq;

    $seq_num++;
}

my ($bg_seq_gc_content, $bg_seq_len) = calc_seqs_total_gc_content($bg_seqs);

$logger->debug("Background seq GC content: $bg_seq_gc_content");
$logger->debug("Background seq total length: $bg_seq_len");

#my $tf_file = "$results_dir/matrices.txt";

my $matrix_set;
my $tf_select_criteria;
if ($tf_file) {
    $logger->info("Reading matrices from $tf_file");
    unless (-e $tf_file) {
        fatal(
            "Specified input TF profile file $tf_file does not exist"
        );
    }

    $matrix_set = read_matrices($tf_file);

    unless ($matrix_set && $matrix_set->size > 0) {
        fatal("Error reading TF profiles from $tf_file");
    }
} else {
    #$tf_file = undef;
    $logger->info("Fetching matrices from JASPAR");

    unless ($tf_db) {
        fatal("No JASPAR database specified");
    }

    unless ($collection) {
        fatal("No JASPAR collection specified");
    }

    #
    # Connect to JASPAR database
    #
    my $jdb = TFBS::DB::JASPAR5->connect(
        "dbi:mysql:" . $tf_db . ":" . JASPAR_DB_HOST,
        JASPAR_DB_USER,
        JASPAR_DB_PASS
    );

    fatal("Could not connect to JASPAR database $tf_db") if !$jdb;

    my %get_matrix_args = (
        -collection => $collection,
        -matrixtype => 'PFM'
    );

    if (@tf_ids) {
        # This takes precendence over tax groups and min IC
        $tf_select_criteria = 'specific';
        $get_matrix_args{-ID} = \@tf_ids;
    } else {
        $tf_select_criteria = 'min_ic';

        $get_matrix_args{-min_ic} = $min_ic || DFLT_CORE_MIN_IC;

        if (@tax_groups) {
            $get_matrix_args{-tax_group} = \@tax_groups;
        }
    }

    $logger->debug("Fetch matrix arguments:\n"
        . Data::Dumper::Dumper(%get_matrix_args));

    $matrix_set = $jdb->get_MatrixSet(%get_matrix_args);

    unless ($matrix_set && $matrix_set->size > 0) {
        fatal("Error fetching TF profiles from JASPAR DB");
    }
}

$logger->debug("Matrix set:\n" . Data::Dumper::Dumper($matrix_set));

my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

$logger->info("Searching target sequences for TFBSs");
#my ($t_tf_seq_siteset, $t_tf_seqs)
my $t_tf_seq_siteset
    = tf_set_search_seqs($tf_set, \%t_seq_id_seqs, $threshold);

fatal("No TFBSs found in test sequences") if !$t_tf_seq_siteset;

$logger->info("Searching background sequences for TFBSs");
#my ($bg_tf_seq_siteset, $bg_tf_seqs)
my $bg_tf_seq_siteset
    = tf_set_search_seqs($tf_set, \%bg_seq_id_seqs, $threshold);

fatal("No TFBSs found in control sequences") if !$bg_tf_seq_siteset;

$logger->info("Computing target sequence TFBS counts");
my $t_counts = compute_site_counts($tf_set, \@t_seq_ids, $t_tf_seq_siteset);

$logger->debug("Target TFBS counts:\n" . Data::Dumper::Dumper($t_counts));

$logger->info("Computing background sequence TFBS counts");
my $bg_counts = compute_site_counts($tf_set, \@bg_seq_ids, $bg_tf_seq_siteset);

#$t_counts->conserved_region_length_set($t_cr_len_set);
#$bg_counts->conserved_region_length_set($bg_cr_len_set);

#my $fresults = fisher_analysis($bg_counts, $t_counts, \%tfs,
#                               scalar @$t_seqs, scalar @$bg_seqs);
my $fisher = OPOSSUM::Analysis::Fisher->new();
$logger->info("Computing Fisher scores");
my $fresults = $fisher->calculate_Fisher_probability($bg_counts, $t_counts);

$logger->debug("Fisher results:\n" . Data::Dumper::Dumper($fresults));

#my $zresults = zscore_analysis($bg_counts, $t_counts, \%tfs,
#                               $t_seq_len, $bg_seq_len);
my $zscore = OPOSSUM::Analysis::Zscore->new();
$logger->info("Computing z-scores");
my $zresults = $zscore->calculate_Zscore(
    $bg_counts, $t_counts, $bg_seq_len, $t_seq_len, $tf_set
);

$logger->debug("Z-score results:\n" . Data::Dumper::Dumper($zresults));

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

    #
    # Fisher scores are now -ln() values so all results are reverse sorted
    # (largest value first).
    # DJA 14/03/2011
    #
    #if ($sort_by eq 'zscore') {
    #    # Sort z-score from highest to lowest
    #    $result_params{-reverse} = 1;
    #}
    $result_params{-reverse} = 1;
}

my $cresult_list = $cresults->get_list(%result_params);

$logger->info("Writing html results");
write_results_html($tf_set, $cresult_list);

if ($cresult_list && $cresult_list->[0]) {
    $logger->info("Writing text results");
    write_results_text($tf_set, $cresult_list);

    $logger->info("Writing TFBS details");
    write_tfbs_details($t_tf_seq_siteset, \%t_seq_id_display_ids);
}

$logger->info("Sending notification email to $email");
send_email($email) if $email;

$logger->info("Finished analysis");

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
    my ($tf_set, $seq_id_seqs, $threshold) = @_;

    my $tf_ids = $tf_set->ids();

    my @seq_ids = keys %$seq_id_seqs;

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

sub calc_seqs_total_gc_content
{
    my ($seqs) = @_;

    my %count = (
        'A' => 0,
        'C' => 0,
        'G' => 0,
        'T' => 0,
        'N' => 0
    );

    my $total_length = 0;
    foreach my $seq (@$seqs) {
        my @nts = split //, $seq->seq();

        foreach my $nt (@nts) {
            $count{uc $nt}++;
        }
    }

    my $gc_count    = $count{'G'} + $count{'C'};
    my $total_count = $gc_count + $count{'A'} + $count{'T'};
    my $gc_content  = $gc_count / $total_count;

    return ($gc_content, $total_count);
}

#
# Ouput combined Z-score/Fisher results to a text file
#
sub write_results_text
{
    my ($tf_set, $results) = @_;

    return unless $results && $results->[0];

    my $filename = "$results_dir/" . RESULTS_TEXT_FILENAME;

    open(FH, ">$filename")
        || fatal("Could not create analysis results file $filename");

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
        printf FH "TF Name\tJASPAR ID\tClass\tFamily\tTax Group\tIC\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tBackground TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

        foreach my $result (@$results) {
            my $tf = $tf_set->get_tf($result->id());

            my $total_ic;
            if ($tf->isa("TFBS::Matrix::PFM")) {
                $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
            } else {
                $total_ic = 'N/A';
            }

            printf FH 
                "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                $tf->name(),
                $tf->ID(),
                $tf->class() || 'N/A',
                $tf->tag('family') || 'N/A',
                $tf->tag('tax_group') || 'N/A',
                $total_ic,
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
        printf FH "TF Name\tClass\tFamily\tTax Group\tIC\tTarget seq hits\tTarget seq non-hits\tBackground seq hits\tBackground seq non-hits\tTarget TFBS hits\tBackground TFBS hits\tTarget TFBS nucleotide rate\tBackground TFBS nucleotide rate\tZ-score\tFisher score\n";

        foreach my $result (@$results) {
            my $tf = $tf_set->get_tf($result->id());

            my $total_ic;
            if ($tf->isa("TFBS::Matrix::PFM")) {
                $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
            } else {
                $total_ic = 'N/A';
            }

            printf FH 
                "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
                $tf->name(),
                $tf->class() || 'N/A',
                $tf->tag('family') || 'N/A',
                $tf->tag('tax_group') || 'N/A',
                $total_ic,
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
        bg_color_class          => BG_COLOR_CLASS,
        heading                 => HEADING,
        title                   => TITLE,
        section                 => 'Analysis Results',
        result_retain_days      => REMOVE_RESULTFILES_OLDER_THAN,
        low_matrix_ic           => LOW_MATRIX_IC,
        high_matrix_ic          => HIGH_MATRIX_IC,
        low_matrix_gc           => LOW_MATRIX_GC,
        high_matrix_gc          => HIGH_MATRIX_GC,
        low_seq_gc              => LOW_SEQ_GC,
        high_seq_gc             => HIGH_SEQ_GC,
        job_id                  => $job_id,
        num_t_seqs              => scalar @$t_seqs,
        num_bg_seqs             => scalar @$bg_seqs,
        t_seq_len               => $t_seq_len,
        bg_seq_len              => $bg_seq_len,
        t_seq_gc_content        => sprintf("%.3f", $t_seq_gc_content),
        bg_seq_gc_content       => sprintf("%.3f", $bg_seq_gc_content),
        tf_db                   => $tf_db,
        tf_select_criteria      => $tf_select_criteria,
        tf_file                 => $tf_file,
        tf_set                  => $tf_set,
        collection              => $collection,
        min_ic                  => $min_ic,
        tax_groups              => \@tax_groups,
        threshold               => $threshold,
        email                   => $email,
        results                 => $results,
        rel_results_dir         => $rel_results_dir,
        result_type             => $result_type,
        num_display_results     => $num_results,
        zscore_cutoff           => $zscore_cutoff,
        fisher_cutoff           => $fisher_cutoff,
        result_sort_by          => $sort_by,
        warn_zero_bg_gene_hits  => $warn_zero_bg_gene_hits,
        results_file            => RESULTS_TEXT_FILENAME,
        user_seq_file           => $user_seq_file,
        user_bg_seq_file        => $user_bg_seq_file,
        user_tf_file            => $user_tf_file,
        bg_seq_set_name         => $bg_seq_set_name,

        formatf                 => sub {
                                        my $dec = shift;
                                        my $f = shift;
                                        return ($f || $f eq '0')
                                            ? sprintf("%.*f", $dec, $f)
                                            : 'N/A'
                                   },

        formatg                 => sub {
                                        my $dec = shift;
                                        my $f = shift;
                                        return ($f || $f eq '0')
                                            ? sprintf("%.*g", $dec, $f)
                                            : 'N/A'
                                   },

        var_template            => "results_seq_ssa.html"
    };

    my $output = process_template('master.html', $vars);

    my $html_filename = "$results_dir/" . RESULTS_HTDOCS_FILENAME;

    open(OUT, ">$html_filename")
        || fatal("Could not create HTML results file $html_filename");

    print OUT $output;

    close(OUT);

    $logger->info("Wrote HTML formatted results to $html_filename");

    return $html_filename;
}

#
# For each TF/seq, write the details of the putative TFBSs out to text and
# html files.
#
sub write_tfbs_details
{
    my ($tf_seq_siteset, $seq_id_display_ids) = @_;

    my $tf_ids = $tf_set->ids();

    foreach my $tf_id (@$tf_ids) {
        my $seq_siteset = $tf_seq_siteset->{$tf_id};

        if ($seq_siteset) {
            my @seq_ids = keys %$seq_siteset;

            my $tf = $tf_set->get_tf($tf_id);

            my $text_filename = sprintf "$results_dir/$tf_id.txt";
            my $html_filename = sprintf "$results_dir/$tf_id.html";
        
            write_tfbs_details_text(
                $text_filename,
                $tf,
                \@seq_ids,
                $seq_id_display_ids,
                $seq_siteset
            );

            write_tfbs_details_html(
                $html_filename,
                $tf,
                \@seq_ids,
                $seq_id_display_ids,
                $seq_siteset
            );
        }
    }
}

#
# Write the details of the putative TFBSs for each TF/seq. Create a
# text file for each TF.
#
sub write_tfbs_details_text
{
    my ($filename, $tf, $seq_ids, $seq_id_display_ids, $seq_siteset) = @_;

    my $tf_id   = $tf->ID();
    my $tf_name = $tf->name();

    open(FH, ">$filename") || fatal(
        "Could not create TFBS details text file $filename"
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    my $total_ic;
    if ($tf->isa("TFBS::Matrix::PFM")) {
        $total_ic = sprintf("%.3f", $tf->to_ICM->total_ic());
    } else {
        $total_ic = 'N/A';
    }

    print  FH "$tf_name\n\n";

    if ($tf_db) {
        print  FH "JASPAR ID:\t$tf_id\n";
    }
    printf FH "Class:\t%s\n", $tf->class() || 'N/A',
    printf FH "Family:\t%s\n", $tf->tag('family') || 'N/A',
    printf FH "Sysgroup:\t%s\n", $tf->tag('tax_group') || 'N/A',
    printf FH "IC:\t%s\n", $total_ic;

    print FH "\n\nConserved $tf_name Binding Sites\n\n";

    #printf FH "\n\n%-31s\t%7s\t%7s\t%7s\t%7s\t%7s\t%s\n",
    printf FH "\n\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        'Sequence ID', 'Start', 'End', 'Strand', 'Score', '%Score', 'TFBS Sequence';

    foreach my $seq_id (@$seq_ids) {
        my $siteset = $seq_siteset->{$seq_id};

        my $display_id = $seq_id_display_ids->{$seq_id};

        next if !$siteset || $siteset->size == 0;

        my $first = 1;
        printf FH "%s", $display_id;

        my $iter = $siteset->Iterator(-sort_by => 'start');
        while (my $site = $iter->next()) {
            # Do not reprint sequence ID for subsequent sites
            unless ($first) {
                printf FH "%s", '';
            }
            $first = 0;

            #printf FH "\t%7d\t%7d\t%7d\t%7.3f\t%6.1f%%\t%s\n",
            printf FH "\t%d\t%d\t%d\t%.3f\t%.1f%%\t%s\n",
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
    my ($filename, $tf, $seq_ids, $seq_id_display_ids, $seq_siteset) = @_;

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
        "Could not create TFBS details HTML file $filename"
    );

    $logger->info("Writing '$tf_name' TFBS details to $filename");

    my $vars = {
        abs_htdocs_path         => ABS_HTDOCS_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_htdocs_path         => REL_HTDOCS_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        version                 => VERSION,
        devel_version           => DEVEL_VERSION,
        jaspar_url              => JASPAR_URL,
        heading                 => HEADING,
        title                   => TITLE,
        section                 => "$tf_name Binding Site Details",
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
                                        : 'N/A'
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

    my $output = process_template('master.html', $vars);

    print FH $output;

    close(FH);
}

sub read_matrices
{
    my ($file) = @_;

    open(FH, $file) || fatal("Could not open matrix file $file - $!");

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

                #
                # Simplistic determination of whether matrix looks more like
                # a PWM than a PFM.
                #
                my $matrix_type = 'PFM';
                if ($matrix_string =~ /\d*\.\d+/) {
                    $matrix_type = 'PWM';
                }

                my $matrix;
                if ($matrix_type eq 'PWM') {
                    $matrix = TFBS::Matrix::PWM->new(
                        -ID           => $id,
                        -name         => $name,
                        -matrixstring => $matrix_string
                    );
                } else {
                    $matrix = TFBS::Matrix::PFM->new(
                        -ID           => $id,
                        -name         => $name,
                        -matrixstring => $matrix_string
                    );
                }

                $matrix_set->add_matrix($matrix);

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
        EVAL_PERL       => 1,   # evaluate Perl code blocks
        DEBUG           => DEBUG
    };

    my $string   = '';
    my $template = Template->new($config);
    my $input    = ABS_HTDOCS_TEMPLATE_PATH . "/$template_name";

    $template->process($input, $vars, \$string)
        || fatal(
            "Error processing template $input\n" . $template->error() . "\n\n"
        );

    return $string;
}

sub send_email
{
    my ($email) = @_;

    return if !$email;

    my $num_seqs    = scalar @$t_seqs;
    my $num_bg_seqs = scalar @$bg_seqs;

    my $results_url = sprintf "%s%s/%s",
        WEB_SERVER_URL,
        "$rel_results_dir",
        RESULTS_HTDOCS_FILENAME;

    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg = "Your oPOSSUM sequence SSA results are now available at\n\n";
    $msg .= "$results_url\n";

    $msg .= "\nAnalysis Summary\n\n";

    $msg .= "Job ID:                         $job_id\n";
    $msg .= "Target sequence file            $user_seq_file\n"
        if $user_seq_file;
    $msg .= "Background sequence file        $user_bg_seq_file\n"
        if $user_bg_seq_file;
    $msg .= "Number of target sequences:     $num_seqs\n";
    $msg .= "Number of background sequences: $num_bg_seqs\n";

    if ($tf_db) {
        $msg .= "TFBS profile source:            JASPAR\n";
        $msg .= "JASPAR collection:              $collection\n" if $collection;
        $msg .= "Min. IC                         $min_ic\n" if $min_ic;
        $msg .= "Taxonomic supergroups           $tax_groups_str\n"
            if $tax_groups_str;
    } else {
        $msg .= "TFBS profile source:            User supplied matrices\n";
    }

    $msg .= "TFBS matrix score threshold:      $threshold\n";

    $msg .= "Results returned:               ";
    if (defined $zscore_cutoff || defined $fisher_cutoff) {
        $msg .= "All results with a z-score >= $zscore_cutoff\n";
        if (defined $fisher_cutoff) {
            $msg .= " and a Fisher score <= $fisher_cutoff\n";
        }
    } else {
        if ($num_results eq 'All') {
            $msg .= "All results";
        } else {
            $msg .= "Top $num_results results";
        }

        $msg .= " sorted by";
        if ($sort_by =~ /zscore/) {
            $msg .= " z-score\n";
        } elsif ($sort_by =~ /fisher/) {
            $msg .= " Fisher score\n";
        }
    }

    $msg .= "\nYour analysis results will be kept on our server for "
            . REMOVE_RESULTFILES_OLDER_THAN . " days.\n";
    $msg .= "\nThank-you,\n";
    $msg .= "The oPOSSUM development team\n";
    $msg .= ADMIN_EMAIL . "\n";

    if (!open(SM, "|" . $cmd)) {
        $logger->error("Could not open sendmail - $!");
        return;
    }

    printf SM "To: %s\n", $email;
    printf SM "From: %s\n", ADMIN_EMAIL;
    print SM "Subject: oPOSSUM Sequence SSA results $job_id\n\n";
    print SM "$msg" ;

    close(SM);
}

sub fatal
{
    my ($error) = @_;

    $error = 'Unknown error' if !$error;
    
    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg = "oPOSSUM sequence-based SSA analysis failed\n";
    $msg .= "\nError: $error\n";
    $msg .= "\nJob ID: $job_id\n";
    $msg .= "\nUser e-mail: $email\n" if $email;

    if (open(SM, "|" . $cmd)) {
        printf SM "To: %s\n", ADMIN_EMAIL;
        print SM "Subject: oPOSSUM sequence-based SSA fatal error\n\n";
        print SM "$msg" ;

        close(SM);
    } else {
        $logger->error("Could not open sendmail - $!") if $logger;
    }

    $logger->logdie("$error") if $logger;
}
