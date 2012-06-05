=head1 NAME

 OPOSSUM::Include::BaseInclude.pm

=head1 SYNOPSIS

=head1 DESCRIPTION

  Contains all options and routines that are common to all the analyses.

=head1 AUTHOR

  Andrew Kwon & David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca, dave@cmmt.ubc.ca

=cut

use strict;

use OPOSSUM::Web::Opt::BaseOpt;

#use Data::Dumper;    # for debugging only
use Carp;

use Template;
use CGI;

use TFBS::DB::JASPAR5;

use OPOSSUM::DBSQL::DBAdaptor;
use OPOSSUM::TFSet;
use OPOSSUM::ConservedTFBS;
use OPOSSUM::Analysis::Counts;


#
# Get log file
#
sub get_log_filename
{
	my ($analysis_type, $results_dir) = @_;

	my $USER = $ENV{'USER'};

	#my $log_dir;
	#if ($USER && $USER ne 'nobody' && $USER ne 'apache') {
	#	$log_dir = "/tmp";
	#} else {
	#	$log_dir = $results_dir;
	#}
    my $log_dir = $results_dir;

	my $log_file = "$log_dir/$analysis_type";
	#$log_file .= "_devel" if DEVEL_VERSION;
	#$log_file .= "_$USER" if $USER;
	$log_file .= ".log";

	return $log_file;
}


#
# Connect to JASPAR database
#
sub jaspar_db_connect
{
    my ($tf_db) = @_;
    
    my $jdb = TFBS::DB::JASPAR5->connect(
        "dbi:mysql:" . $tf_db . ":" . JASPAR_DB_HOST,
        JASPAR_DB_USER,
        JASPAR_DB_PASS
    );

    return $jdb;
}

sub read_tf_ids_from_file
{
    my ($file, %job_args) = @_;

    return read_ids_from_file($file, %job_args);
}

sub read_tf_families_from_file
{
    my ($file, %job_args) = @_;

    return read_ids_from_file($file, %job_args);
}

#
# Generic read IDs (gene or TF IDs) from file.
#
sub read_ids_from_file
{
    my ($file, %job_args) = @_;

    my $text = read_file($file, %job_args);
    my $ids = parse_id_text($text);

    return $ids;
}

sub read_file
{
    my ($file, %job_args) = @_;

    open(FH, $file) || fatal("Could not open file $file", %job_args);

    my $text = "";
    while (my $line = <FH>) {
        $text .= $line;
    }

    close(FH);

    return $text;
}

sub parse_id_text
{
    my $text = shift;

    #
    # Strip anything out that is NOT a part of the ID or a valid
    # separator
    #
    $text =~ s/[^\w\.\/_\-,;:\s\n]+//g;

    #
    # Strip out leading and trailing separators
    #
    $text =~ s/^[,;:\s\n]+//g;
    $text =~ s/[,;:\s\n]+$//g;

    my @raw_list = split /[,;:\n\s]+/, $text;

    my %included;
    my @unique_list;
    if (@raw_list) {
        foreach my $id (@raw_list) {
            unless ($included{$id}) {
                push @unique_list, $id;
                $included{$id} = 1;
            }
        }
    }

    return @unique_list ? \@unique_list : undef;
}

sub revcom
{
	my ($seq) = @_;

	my $rc_seq = reverse $seq;

	$rc_seq =~ tr/[acgtACGT]/[tgcaTGCA]/;

	return $rc_seq;
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

=head2 matrix_set_compute_gc_content

 Title   : matrix_set_compute_gc_content

 Function: Compute the GC content of each of the matrices in a set of TFBS
           matrices and set the matrix tag value 'gc_content' to the value
           computed.

 Args    : matrix_set   - a TFBS::MatrixSet object

 Returns : 1 on success, otherwise undef.

=cut

sub matrix_set_compute_gc_content
{
    my $matrix_set = shift;

    unless ($matrix_set && $matrix_set->size) {
        return undef;
    }

    my $iter = $matrix_set->Iterator;
    while (my $matrix = $iter->next) {
        my $gc_content = matrix_compute_gc_content($matrix);

        if (defined $gc_content) {
            $matrix->tag('gc_content', $gc_content);
        }
    }

    return 1;
}

=head2 matrix_compute_gc_content

 Title   : matrix_compute_gc_content

 Function: Compute the GC content of a TFBS matrix.

 Args    : pfm  - a TFBS::Matrix::PFM object

 Returns : On success, the GC content of the matrix in the range 0 - 1,
           otherwise undef.

=cut

sub matrix_compute_gc_content
{
    my $pfm = shift;

    unless ($pfm->isa("TFBS::Matrix::PFM")) {
        carp "Cannot compute GC content for non-PFM matrix\n"; 
        return undef;
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

sub tfbss_to_conserved_tfbss
{
    my ($sites, $cluster_id, $seq_id) = @_;
    
    return if !defined $sites || scalar($sites) == 0;
    
    my @ctfbss;
    foreach my $site (@$sites)
    {
        my $ctfbs = OPOSSUM::ConservedTFBS->new(
            -tf_id      => $cluster_id,
            -gene_id    => $seq_id,
            -start      => $site->start,
            -end        => $site->end,
            -strand     => $site->strand,
            -score      => $site->score,
            -rel_score  => $site->rel_score,
            -seq        => $site->seq->seq
        );
        
        push @ctfbss, $ctfbs;
    }
    
    return @ctfbss ? \@ctfbss : undef;
}


sub process_template
{
    my ($template_name, $vars, %job_args) = @_;

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
            "Error processing template $input\n" . $template->error() . "\n\n",
            %job_args
        );

    return $string;
}

sub send_email
{
    my (%args) = @_;

    my $job_id = $args{-job_id};
    my $heading = $args{-heading};
    my $email = $args{-email};
    my $web = $args{-web};
    my $abs_results_dir = $args{-abs_results_dir};
    my $rel_results_dir = $args{-rel_results_dir};
    my $user_t_file = $args{-user_t_file};
    my $user_bg_file = $args{-user_bg_peak_file};
    my $user_t_peak_file = $args{-user_t_peak_file};
    my $user_bg_peak_file = $args{-user_bg_file};
    my $t_num = $args{-t_num};
    my $bg_num = $args{-bg_num};
    my $tf_db = $args{-tf_db};
    my $collections = $args{-collections};
    my $tax_groups = $args{-tax_groups};
    my $fam_file = $args{-families};
    my $tf_ids = $args{-tf_ids};
    my $min_ic = $args{-min_ic};
    my $threshold = $args{-threshold};
    my $z_cutoff = $args{-z_cutoff};
    my $f_cutoff = $args{-f_cutoff};
    my $ks_cutoff = $args{-ks_cutoff};
    my $num_results = $args{-num_results};
    my $sort_by = $args{-sort_by};
    my $logger = $args{-logger};

    return if !$email;


    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg .= "\n";
    $msg = "Your oPOSSUM $heading results are now available at\n\n";
    if ($web) {
        my $results_url = sprintf "%s%s/%s",
            WEB_SERVER_URL,
            "$rel_results_dir",
            RESULTS_HTDOCS_FILENAME;
        $msg .= "$results_url\n\n";
    } else {
        $msg .= "$abs_results_dir\n\n";
    }
    
    $msg .= "\nAnalysis Summary\n\n";

    $msg .= "Job ID:                                $job_id\n";
    $msg .= "Target gene/sequence file:             $user_t_file\n"
        if $user_t_file;
    $msg .= "Background gene/sequence file:         $user_bg_file\n"
        if $user_bg_file;
    $msg .= "Target sequence peak position file:    $user_t_peak_file\n"
        if $user_t_peak_file;
    $msg .= "Background sequence peak positionfile: $user_bg_peak_file\n"
        if $user_bg_peak_file;
    $msg .= "Number of target genes/sequences:      $t_num\n";
    $msg .= "Number of background genes/sequences:  $bg_num\n";

    if ($tf_db) {
        $msg .= "TFBS profile source:                   JASPAR\n";
        $msg .= "JASPAR collection:                     $collections\n"
            if $collections;
        $msg .= "Taxonomic supergroups:                 $tax_groups\n"
            if $tax_groups;
        $msg .= "Structural family file:                $fam_file\n"
            if $fam_file;
        $msg .= "TF ID's:                               $tf_ids\n"
            if $tf_ids;
        $msg .= "Min. IC                                $min_ic\n"
            if $min_ic;
    } else {
        $msg .= "TFBS profile source:                   User supplied matrices\n";
    }

    $msg .= "TFBS matrix score threshold:    $threshold\n" if $threshold;

    $msg .= "Results returned:                      ";

    if (defined $z_cutoff || defined $f_cutoff || defined $ks_cutoff) {
        $msg .= "All results with a z-score >= $z_cutoff\n";
        if (defined $f_cutoff) {
            $msg .= " and a Fisher score >= $f_cutoff\n";
        }
        if (defined $ks_cutoff) {
            $msg .= " and a KS p-value <= $ks_cutoff\n";
        }
    } else {
        if (!$num_results || $num_results =~ /^all/i) {
            $msg .= "All results";
        } else {
            $msg .= "Top $num_results results";
        }
        
        if (defined $sort_by) {
            $msg .= " sorted by";
            if ($sort_by =~ /zscore/) {
                $msg .= " z-score\n";
            } elsif ($sort_by =~ /fisher/) {
                $msg .= " Fisher score\n";
            } elsif ($sort_by =~ /ks/) {
                $msg .= " KS p-value\n";
            }
        }
    }

    $msg .= "\n";
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
    print SM "Subject: oPOSSUM $heading results\n\n";
    print SM "$msg" ;

    close(SM);
}

sub fatal
{
    my ($error, %args) = @_;
#    my ($error, $job_id, $heading, $email, $logger) = @_;

    my $job_id = $args{-job_id};
    my $heading = $args{-heading};
    my $email = $args{-email};
    my $logger = $args{-logger};
    my $web = $args{-web};
    my $abs_results_dir = $args{-abs_results_dir};
    
    $error = 'Unknown error' unless $error;

    my $cmd = "/usr/sbin/sendmail -i -t";

    my $msg = "oPOSSUM $heading analysis failed\n";
    $msg .= "\nJob ID: $job_id\n";
    $msg .= "\nError: $error\n";

    if (open(SM, "|" . $cmd)) {
        printf SM "To: %s\n", ADMIN_EMAIL;
        print SM "Subject: oPOSSUM $heading fatal error\n\n";
        print SM "$msg" ;
        print SM "\nUser e-mail: $email\n" if $email;

        close(SM);
    } else {
        $logger->error("Could not open sendmail - $!") if $logger;
    }

    if ($email) {
        if (open(SM, "|" . $cmd)) {
            printf SM "To: %s\n", $email;
            printf SM "From: %s\n", ADMIN_EMAIL;
            print SM "Subject: oPOSSUM $heading fatal error\n\n";
            print SM "$msg" ;

            close(SM);
        }
    }

    #
    # If in web context and we know the results dir, then create an error
    # html page in the results directory in place of the results html.
    # Currently this is only useful for gene based SSA which displays the
    # results html upon completion.
    # DJA 2012/04/25
    #
    if ($web && $abs_results_dir) {
        my $vars = {
            abs_htdocs_path  => ABS_HTDOCS_PATH,
            rel_htdocs_path  => REL_HTDOCS_PATH,
            abs_cgi_bin_path => ABS_CGI_BIN_PATH,
            rel_cgi_bin_path => REL_CGI_BIN_PATH,
            #bg_color_class   => BG_COLOR_CLASS,
            version          => VERSION,
            devel_version    => DEVEL_VERSION,
            title            => 'Error',
            error            => $error,
            var_template     => "error.html"
        };

        #
        # DO NOT call proccess_template routine here as that routine may call
        # this fatal routine, resulting in an infinite loop!
        #
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
        my $input    = ABS_HTDOCS_TEMPLATE_PATH . "/master.html";

        if ($template->process($input, $vars, \$string)) {
            my $error_html_file = "$abs_results_dir/" . RESULTS_HTDOCS_FILENAME;

            if (open(EHFH, ">$error_html_file")) {
                print EHFH $string;

                close(EHFH);
            }
        }
    }

    $logger->logdie("$error") if $logger;
}

1;

