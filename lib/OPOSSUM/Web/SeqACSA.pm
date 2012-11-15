package SeqACSA;

use base 'CGI::Application';

use OPOSSUM::Opt::ACSAOpt;

use OPOSSUM::Web::Include::BaseInclude;
use OPOSSUM::Web::Include::SeqInclude;

use Data::Dumper;    # for debugging only
use Template;
use CGI::Carp qw(carpout);    # fatalsToBrowser;
use File::Temp qw/ tempdir /;

use OPOSSUM::Web::State;

use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_seq_acsa';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_seq_acsa";
$log_file .= "_devel" if DEVEL_VERSION;
$log_file .= "_$USER" if $USER;
$log_file .= ".log";

open(LOG, ">>$log_file") || die "Error opening log file $log_file - $!\n";

carpout(\*LOG);

sub setup
{
    my $self = shift;

    #print STDERR "setup\n";

    $self->start_mode('input');
    $self->mode_param('rm');
    $self->run_modes(
        'input'             => 'input',
        'process'           => 'process'
    );

    my $q = $self->query();

    my $sid = $q->param('sid');

    my $state;
    if ($sid) {
        #
        # Existing session. Load state from file.
        #
        my $filename = _session_tmp_file($sid);
        $state = OPOSSUM::Web::State->new(__Fn => $filename);
    } else {
        #
        # New session. Create new session ID and state object.
        #
        $sid = $$ . time;
        my $filename = _session_tmp_file($sid);

        $state = OPOSSUM::Web::State->new(
            -sid => $sid,
            __Fn => $filename
        );

        $self->initialize_state($state);
    }

    $self->state($state);

    $self->errors(undef);
    $self->warnings(undef);

    #printf STDERR sprintf("\n\noPOSSUM State:\n%s\n\n",
    #    Data::Dumper::Dumper($self->state())
    #);

	#unless ($self->opossum_db_connect()) {
	#    return $self->error("Could not connect to oPOSSUM DB");
	#}

    #printf STDERR "\n\nrun mode = %s\n\n", $q->param('rm');
}

sub teardown
{
    my $self = shift;

    #print STDERR "teardown\n";

    my $state = $self->state();
    if ($state) {
        $state->dumper->Purity(1);
        $state->dumper->Deepcopy(1);
        $state->commit();
    }

    $self->_clean_tempfiles;
    $self->_clean_resultfiles;
}

sub input
{
    my $self = shift;

    #print STDERR "input\n";

    my $q = $self->query;
    #print STDERR "input query:\n"
    #		    . Data::Dumper::Dumper($q);

    my $state = $self->state();

	my $min_ic		= DFLT_CORE_MIN_IC;
	my $tax_groups	= JASPAR_CORE_TAX_GROUPS;

	my $num_tax_groups = scalar @$tax_groups;

    #
    # Connect to JASPAR DB and retrieve TF info
    #
    $self->jaspar_db_connect();

    my $tf_set = $self->fetch_tf_set(
        -collection => 'CORE',
        -tax_group  => $tax_groups,
        -min_ic     => $min_ic
    );

    my %tax_tf_sets;
    foreach my $tax_group (@$tax_groups) {
        $tax_tf_sets{$tax_group} = $tf_set->subset(-tax_groups => $tax_group);
    }

    my $sample_seqs = qq{>seq1
AAGAAAATAAAGGTAACTCAATTTAGCTCCCTGGGATCCTGACTTACATACAACTGCATA
TTTGTGATGACTGTTTTCCATCCAATAATTATATACTTTTTTTCTAGGAAACAGAAGGCA
ACTTATTGCGCTTTTCTTTTTTCCTTGTAAGTGTTTTCAGGCAGAGATTTTTTTTTTGGA
AGGTATTTCCCTTCCCCTGTG
>seq2
CATATAGTCAATACAATGTGTTTAATAAAATCAATTTATTTCCTTTTAGTCAACTAATAA
TCATGCCCTCACACAAATTATTTGCAACATTTGATTGAAAATATTCATATTCCAAGCATC
TCTGTTCTCTGACATATGTTTAAAGTGGGTTCTAAGCTATAAAACGAAACTGCTAACATA
TTAGCTTATTACAATCTGACA
>seq3
GTAGCCAAGTTGAGTGCTGGGCCAGCACTGAGTCATATGCTGATTGCATCAAAAGTAACT
TCACAAATGCTTTTCTTTCCCTGTGTGTGCTCAGTCATGCTTACATTACTGCTTTCTGAT
AATCCTGGCTAAAAAGATATAAGAGACACTCATGGAAATTAAAAGGTGATTTACAGTGCA
CCAAAGATGATGAATAGAATA
>seq4
GAAATCATATCTCTGTGCCTATTTCCTCATCTCTAGAACAAATGAATGCGTTGCCAGGCA
GCACTTAATGAGAGTACAACTGCCTTGTAAGGAAAGCATGACCAGGCAGTTTCACTGTTG
GGTGAACCTCACAGAGGCATTGAGGCATTGGCAGATACAGCCTCCTAATACTGCCGAAAG
ATACAATAAGAATATATATAT
>seq5
CATTTGGGACAAAGACTGCAGAGGAAGAAGATAATGGAGATTGGAGCCTGACTCAGCACA
GACTCTTTGGAGCATCCCTGTGAGACATTGATGCTGACTCATTGCTCGGGGCTAGTTTCC
AAGCATATAAACAGAAACAGTAAGGTCGGCTGAAAATCCACAGAGCTTTAGTGTACAAAG
CACTAAAGAAAAACAAGGGCA};
	
    #
    # Format a tax group list for display on the web page. The $tax_group
    # variable is obtained from the db and may be a single value or a comma
    # separated list (no space after comma).
    #
    my @s_tax_groups;
    foreach my $tg (@$tax_groups) {
        # remove trailing 's' (un-pluralize)
        my $stg = $tg;
        $stg =~ s/s$//;
        push @s_tax_groups, $stg;
    }

    my $tax_group_list = join ' / ', @s_tax_groups;

    my $vars = {
        abs_htdocs_path         => ABS_HTDOCS_PATH,
        rel_htdocs_path         => REL_HTDOCS_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        bg_color_class          => $state->bg_color_class(),
        title                   => $state->title(),
        heading                 => $state->heading(),
        section                 => 'Select Analysis Parameters',
        version                 => VERSION,
        devel_version           => DEVEL_VERSION,
        dflt_inter_binding_dist => DFLT_INTER_BINDING_DIST,
        max_inter_binding_dist  => MAX_INTER_BINDING_DIST,
        nresults                => NUM_RESULTS,
        dflt_nresults           => DFLT_NUM_RESULTS,
        zcutoffs                => ZSCORE_CUTOFFS,
        fcutoffs                => FISHER_CUTOFFS,
        dflt_zcutoff            => DFLT_ZSCORE_CUTOFF,
        dflt_fcutoff            => DFLT_FISHER_CUTOFF,
		#thl_dflt_level          => DFLT_THRESHOLD_LEVEL,
		dflt_threshold			=> DFLT_THRESHOLD,
		min_threshold			=> MIN_THRESHOLD,
		#collections				=> JASPAR_COLLECTIONS,
		bg_seq_set_keys			=> BG_SEQ_SET_KEYS,
		bg_seq_set_names		=> BG_SEQ_SET_NAMES,
        sid                     => $state->sid(),
        tax_groups              => $tax_groups,
        num_tax_groups          => $num_tax_groups,
        tax_group_list          => $tax_group_list,
		min_ic					=> $min_ic,
        tf_set                  => $tf_set,
        tax_tf_sets             => \%tax_tf_sets,
		sample_sequences		=> $sample_seqs,
        var_template            => "input_seq_acsa.html"
    };

    my $output = $self->process_template('master.html', $vars);
    #print STDERR "input results:\n"
    #		    . Data::Dumper::Dumper($output);

    return $output;
}

sub process
{
    my $self = shift;

    my $state = $self->state();

    my $q = $self->query;

    my $dflt_min_ic = DFLT_CORE_MIN_IC;
	my $dflt_tax_groups	= JASPAR_CORE_TAX_GROUPS;
    
	my $num_tax_groups = scalar @$dflt_tax_groups;

	my $tf_db;
	my $collection;
    my @tf_ids;
	my $tf_ids_str;
    my $min_ic;
    my @tax_groups;
    my $tax_groups_str;
	my $matrix_set;
	my $tf_set;

    my $tf_select_method = $q->param("tf_select_method");
    if (!$tf_select_method) {
        return $self->error("TFBS profiles not specified");
    }

	if ($tf_select_method =~ /(\w+?)_/) {
		$collection = uc $1;
		$tf_db = JASPAR_DB_NAME;
	}

	my $tf_select_criteria;
    if ($tf_select_method =~ /min_ic/) {
	    $tf_select_criteria = 'min_ic';
    } elsif ($tf_select_method =~ /specific/) {
	    $tf_select_criteria = 'specific';
    } elsif ($tf_select_method =~ /paste/) {
	    $tf_select_criteria = 'paste';
    } elsif ($tf_select_method =~ /upload/) {
 	    $tf_select_criteria = 'upload';
    } else {
	    $tf_select_criteria = 'all';
    }

	if ($tf_select_criteria eq 'min_ic') {
		if ($collection eq 'CORE') {
			$min_ic = $q->param('core_min_ic');

			if ($num_tax_groups > 1) {
				@tax_groups = $q->param('core_tax_groups');
				if (scalar @tax_groups == 0) {
					return $self->error("No JASPAR CORE tax group selected");
				}
			} else {
				@tax_groups = @$dflt_tax_groups;
			}
		}

		if (!$min_ic) {
			return $self->error("No JASPAR $collection minimum IC specified");
		}

		$tax_groups_str = join(',', @tax_groups);
	} elsif ($tf_select_criteria eq 'specific') {
		if ($collection eq 'CORE') {
			push @tf_ids, $q->param('core_tfs');
			if (scalar(@tf_ids) == 0) {
				return $self->error("No specific $collection TFs selected");
			}
		}

		$tf_ids_str = join ",", @tf_ids;

		# This is just to get the TF names for display on the results
		# summary page
		$self->jaspar_db_connect();
		$tf_set = $self->fetch_tf_set(
			-ID => \@tf_ids
		);
	} elsif ($tf_select_criteria eq 'paste') {
		my $matrix_text = $q->param('matrix_paste_text');

		$matrix_set = $self->parse_matrix_text($matrix_text);

		if (!$matrix_set || $matrix_set->size == 0) {
			return $self->error(
				"Error parsing TFBS profile matrices from pasted text");
		}
        my $matrix_upload_file = $q->param('matrix_upload_file');
        my $matrix_upload_fh   = $q->upload('matrix_upload_file');

        my $matrix_text = $self->upload_matrix_file($matrix_upload_fh);
        if (!$matrix_text) {
            $self->_error("Uploading TFBS profile matrices");
            return;
        }

        $matrix_set = $self->parse_matrix_text($matrix_text);

        if (!$matrix_set || $matrix_set->size == 0) {
            return $self->error(
                "Error parsing TFBS profile matrices from uploaded file"
            );
        }
    }

    #
    # XXX
    # Right now if the selected TFs are JASPAR TFs then the anchoring TF
    # is also considered to be from JASPAR. This may change in the future.
    # XXX
    #
    my $anchor_tf_select_method = $q->param("anchor_tf_select_method");
    if (!$anchor_tf_select_method) {
        return $self->error("No anchoring TFBS profile selected");
    }

    my $anchor_collection;
    my $anchor_tf_db;
    if ($anchor_tf_select_method =~ /(\w+?)_/) {
        $anchor_collection = uc $1;
        $anchor_tf_db = JASPAR_DB_NAME;
    }

    my $anchor_tf_select_criteria;
    if ($anchor_tf_select_method =~ /specific/) {
        $anchor_tf_select_criteria = 'specific';
    } elsif ($anchor_tf_select_method =~ /paste/) {
        $anchor_tf_select_criteria = 'paste';
    } elsif ($anchor_tf_select_method =~ /upload/) {
        $anchor_tf_select_criteria = 'upload';
    }

    my $anchor_tf;
    my $anchor_tf_id;
    my $anchor_matrix;
    if ($anchor_tf_select_criteria eq 'specific') {
        $anchor_tf_id = $q->param('anchor_tf_id');

        unless ($anchor_tf_id) {
            return $self->error("Anchoring TF not specified");
        }

        $self->jaspar_db_connect();

        $anchor_matrix = $self->jdbh->get_Matrix_by_ID(
            $anchor_tf_id, 'PFM'
        );
    } elsif ($anchor_tf_select_criteria eq 'paste') {
        my $anchor_matrix_text = $q->param('anchor_matrix_paste_text');

        #printf STDERR "anchor_matrix_text:\n$anchor_matrix_text\n\n";
    
        my $anchor_matrix_set = $self->parse_matrix_text($anchor_matrix_text);

        if (!$anchor_matrix_set || $anchor_matrix_set->size == 0) {
            return $self->error(
                "Error parsing anchor TFBS profile matrix from pasted text");
        }

        if ($anchor_matrix_set->size > 1) {
            $self->_warning(
                "Pasted more than one custom anchor TFBS profile"
                . " matrix - using the first one read"
            );
        }

        my $iter = $anchor_matrix_set->Iterator();
        $anchor_matrix = $iter->next();
    } elsif ($anchor_tf_select_criteria eq 'upload') {
        my $anchor_matrix_upload_file = $q->param('anchor_matrix_upload_file');
        my $anchor_matrix_upload_fh   = $q->upload('anchor_matrix_upload_file');

        my $anchor_matrix_text = $self->upload_matrix_file(
            $anchor_matrix_upload_fh
        );
        if (!$anchor_matrix_text) {
            $self->_error("Uploading anchor TFBS profile matrix");
            return;
        }

        #printf STDERR "anchor_matrix_text:\n$matrix_text\n\n";

        my $anchor_matrix_set = $self->parse_matrix_text($anchor_matrix_text);

        if (!$anchor_matrix_set || $anchor_matrix_set->size == 0) {
            return $self->error(
                "Error parsing anchor TFBS profile matrix from uploaded file"
            );
        }

        if ($anchor_matrix_set->size > 1) {
            $self->_warning(
                "Uploaded more than one custom anchor TFBS profile"
                . " matrix - using the first one read"
            );
        }

        my $iter = $anchor_matrix_set->Iterator();
        $anchor_matrix = $iter->next();
    }

    my $threshold       = $q->param('threshold');

    my $inter_binding_dist = $q->param('inter_binding_dist');
    if ($inter_binding_dist > MAX_INTER_BINDING_DIST) {
        return $self->error(
            sprintf(
                "Specified inter-binding distance exceeds maximum of %d allowed",
                MAX_INTER_BINDING_DIST
            )
        );
    }

    my $num_display_results;
    my $zscore_cutoff;
    my $fisher_cutoff;
    my $result_type = $q->param('result_type');
    $state->result_type($result_type);
    if ($result_type eq 'top_x_results') {
        $num_display_results = $q->param('num_display_results');
    } elsif ($result_type eq 'significant_hits') {
        $zscore_cutoff = $q->param('zscore_cutoff');
        $fisher_cutoff = $q->param('fisher_cutoff');
    }

    my $result_sort_by = $q->param('result_sort_by');
    $state->result_sort_by($result_sort_by);

    my $email = $q->param('email');
    $state->email($email);

    my $seq_input_method = $q->param('seq_input_method');
    if (!$seq_input_method) {
        return $self->error("No target sequences specified");
    }

    my $bg_seq_input_method = $q->param('bg_seq_input_method');
    if (!$bg_seq_input_method) {
        return $self->error("No background sequence specified");
    }

    my $user_seq_file;
    if ($seq_input_method eq 'upload') {
        $user_seq_file = $q->param('seq_file');

        $user_seq_file =~ s/.*\///;
    }

    my $user_bg_seq_file;
    my $bg_seq_set_name;
    if ($bg_seq_input_method eq 'upload') {
        $user_bg_seq_file = $q->param('bg_seq_file');

        $user_bg_seq_file =~ s/.*\///;
    } elsif ($bg_seq_input_method eq 'default') {
        my $bg_seq_set_names = BG_SEQ_SET_NAMES;
        my $bg_seq_set_key = $q->param('bg_seq_set_key');

        $bg_seq_set_name = $bg_seq_set_names->{$bg_seq_set_key};
    }

    my $user_matrix_file;
    if ($tf_select_criteria eq "upload") {
        $user_matrix_file = $q->param('matrix_upload_file');

        $user_matrix_file =~ s/.*\///;
    }

    my $user_anchor_matrix_file;
    if ($anchor_tf_select_criteria eq "upload") {
        $user_anchor_matrix_file = $q->param('anchor_matrix_upload_file');

        $user_anchor_matrix_file =~ s/.*\///;
    }

    #
    # Create a temporary working directory for all the temp. input files and
    # output results file as a sub-directory under the defined temp. dir.
    #
    my $tempdir = tempdir(DIR => ABS_HTDOCS_RESULTS_PATH);

    my $job_id = $tempdir;
    $job_id =~ s/.*\///;

    my $seq_filename  = $self->get_seq_file(
        $seq_input_method,
        $tempdir
    );

    unless ($seq_filename) {
        if ($seq_input_method eq 'upload') {
            return $self->error(
                "There was a problem with the provided target sequence file"
            );
        } elsif ($seq_input_method eq 'paste') {
            return $self->error(
                "There was a problem with the pasted target sequence text"
            );
        }
    }

    my $bg_seq_filename = $self->get_back_seq_file(
        $bg_seq_input_method,
        $tempdir
    );

    unless ($bg_seq_filename) {
        if ($bg_seq_input_method eq 'upload') {
            return $self->error(
                "There was a problem with the provided background sequence file"
            );
        } elsif ($bg_seq_input_method eq 'paste') {
            return $self->error(
                "There was a problem with the pasted background sequence text"
            );
        } else {
            return $self->error();
        }
    }

    my $matrix_filename;
    if ($tf_select_criteria eq 'paste' || $tf_select_criteria eq 'upload') {
        $matrix_filename = $self->write_matrix_file(
            $matrix_set,
            $tempdir
        );

        return $self->error('Could not write TFBS profile matrices file')
            if !$matrix_filename;
    }

    my $anchor_matrix_filename;
    if ($anchor_tf_select_criteria eq 'paste'
        || $anchor_tf_select_criteria eq 'upload'
    ) {
        $anchor_matrix_filename = $self->write_matrix_file(
            $anchor_matrix,
            $tempdir,
            'anchor_matrices.txt'
        );

        return $self->error('Could not write anchor TFBS profile matrix file')
            if !$anchor_matrix_filename;
    }

	#
	# Call the analysis script
	# 
    my $command = OPOSSUM_SCRIPTS_PATH . "/opossum_seq_acsa.pl"
   		. " -j $job_id"
        . " -tsf $seq_filename"
        . " -bsf $bg_seq_filename"
		. " -d $tempdir"
        . " -dist $inter_binding_dist"
		. " -web"
		. " -m $email";

    if ($matrix_filename) {
        $command .= " -mf $matrix_filename";
    }

	if ($tf_db) {
		$command .= " -tfdb $tf_db";
	}

    if (defined $threshold) {
        $command .= " -th $threshold%";
    }

    if ($collection) {
        $command .= " -co $collection";
    }

    if ($tax_groups_str) {
        $command .= " -tax '$tax_groups_str'";
    }

    if ($min_ic) {
        $command .= " -ic $min_ic";
    }

    if ($tf_ids_str) {
        $command .= " -tfids '$tf_ids_str'";
    }

    if ($anchor_tf_id) {
        $command .= " -atfid '$anchor_tf_id'";
    } elsif ($anchor_matrix_filename) {
        $command .= " -amf $anchor_matrix_filename";
    } else {
        return $self->error("No anchor TF specified");
    }

    if ($result_type eq 'top_x_results') {
        $command .= " -n $num_display_results";
    } elsif ($result_type eq 'significant_hits') {
        $command .= " -zcutoff $zscore_cutoff";
        $command .= " -fcutoff $fisher_cutoff";
    }

    if ($result_sort_by) {
        $command .= " -sr $result_sort_by";
    }

    if ($user_anchor_matrix_file) {
        $command .= " -uamf '$user_anchor_matrix_file'";
    }

    if ($user_matrix_file) {
        $command .= " -umf '$user_matrix_file'";
    }

	if ($user_seq_file) {
		$command .= " -utsf '$user_seq_file'";
	}

	if ($user_bg_seq_file) {
		$command .= " -ubsf '$user_bg_seq_file'";
	}

    my $submitted_time = scalar localtime(time);
    printf LOG "\nStarting aCSA analysis at $submitted_time:\n$command\n\n";

	my $rv = system("$command >/dev/null 2>&1 &");
	if ($rv) {
		printf LOG "Analysis script returned $rv - $!\n";
	}

    my $vars = {
        abs_htdocs_path     => ABS_HTDOCS_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_htdocs_path     => REL_HTDOCS_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        rel_htdocs_tmp_path => REL_HTDOCS_TMP_PATH,
        bg_color_class      => BG_COLOR_CLASS,
        rel_htdocs_results_path => REL_HTDOCS_RESULTS_PATH,
        version             => VERSION,
        devel_version       => DEVEL_VERSION,
        heading             => $state->heading(),
        title               => $state->title(),
        section             => 'Analysis Submitted',
        jaspar_url          => JASPAR_URL,
        job_id              => $job_id,
        submitted_time      => $submitted_time,
        user_seq_file       => $user_seq_file,
        user_bg_seq_file    => $user_bg_seq_file,
        user_matrix_file    => $user_matrix_file,
        bg_seq_set_name     => $bg_seq_set_name,
        tf_select_criteria  => $tf_select_criteria,
        anchor_tf_select_criteria 
                            => $anchor_tf_select_criteria,
        tf_db               => $tf_db,
        matrix_file         => $matrix_filename,
        tf_set              => $tf_set,
        anchor_matrix       => $anchor_matrix,
        user_anchor_matrix_file => $user_anchor_matrix_file,
        collection          => $collection,
        tax_groups          => \@tax_groups,
        min_ic              => $min_ic,
        threshold           => $threshold,
        result_type         => $result_type,
        num_display_results => $num_display_results,
        zscore_cutoff       => $zscore_cutoff,
        fisher_cutoff       => $fisher_cutoff,
        result_sort_by      => $result_sort_by,
        email               => $email,
        var_template        => "analysis_summary_seq_acsa.html"
    };

    #print STDERR "results vars:\n" . Data::Dumper::Dumper($vars);

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

sub initialize_state
{
    my ($self, $state) = @_;

    my $heading = sprintf "Sequence-based Anchored Combination Site Analysis";
    $state->heading($heading);

    $state->debug(DEBUG);
    $state->title("oPOSSUM $heading");
    $state->bg_color_class(BG_COLOR_CLASS);

    #$state->errors(undef);
    #$state->warnings(undef);
}

1;
