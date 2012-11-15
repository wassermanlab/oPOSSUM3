package GeneSSA;

use base 'CGI::Application';

use OPOSSUM::Web::Include::BaseInclude;
use OPOSSUM::Web::Include::GeneInclude;

use Data::Dumper;    # for debugging only
use Template;
use CGI::Carp qw(carpout);    # fatalsToBrowser;
use File::Temp qw/ tempdir /;

use OPOSSUM::Web::State;

use OPOSSUM::Analysis::Fisher;
use OPOSSUM::Analysis::Zscore;
use OPOSSUM::Analysis::CombinedResultSet;


use constant DEBUG          => 0;
use constant BG_COLOR_CLASS => 'bgc_gene_ssa';

use strict;

my $USER = $ENV{'USER'};

my $log_dir;
if (DEVEL_VERSION || ($USER && $USER ne 'nobody' && $USER ne 'apache')) {
    $log_dir = "/tmp";
} else {
    $log_dir = OPOSSUM_LOG_PATH;
}

my $log_file = "$log_dir/oPOSSUM_gene_ssa";
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
        'results'           => 'results',
#        'tfbs_details'      => 'tfbs_details',
#        'text_results'      => 'text_results',
#        'text_tfbs_details' => 'text_tfbs_details'
    );

    my $q = $self->query();
    
    my $sid = $q->param('sid');

    my $state;
    if ($sid) {
        #
        # Existing session. Load state from file.
        #
        my $filename = _session_tmp_file($sid);

        #printf STDERR "%s: loading state\n", scalar localtime;

        $state = OPOSSUM::Web::State->new(__Fn => $filename);
    } else {
        #
        # New session. Create new session ID and state object.
        #

        #printf STDERR "%s: creating state\n", scalar localtime;

        $sid = $$ . time;
        my $filename = _session_tmp_file($sid);

        $state = OPOSSUM::Web::State->new(
            -sid => $sid,
            __Fn => $filename
        );

        #printf STDERR "%s: initializing state\n", scalar localtime;

        $self->initialize_state($state);
    }

    $self->state($state);

    #printf STDERR sprintf("\n\noPOSSUM State:\n%s\n\n",
    #    Data::Dumper::Dumper($self->state())
    #);

    #printf STDERR "%s: connecting to oPOSSUM DB\n", scalar localtime;

    $self->errors(undef);
    $self->warnings(undef);

    unless ($self->opossum_db_connect()) {
        return $self->error("Could not connect to oPOSSUM DB");
    }
    
    #printf STDERR "\n\nrun mode = %s\n\n", $q->param('rm');
}


sub teardown
{
    my $self = shift;

    #print STDERR "teardown\n";

    if ($self->opdba()) {
        if ($self->opdba()->dbc()) {
            $self->opdba()->dbc()->disconnect();
        }
        $self->{-opdba} = undef;
    }

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
    #    . Data::Dumper::Dumper($q);

    my $state = $self->state();

    my $species = $state->species() || $self->param('species');
    
    #
    # Retrieve any previously entered values from state.
    #
    #my $in_t_gene_id_type = $state->t_gene_id_type();
    my @in_t_gene_ids = @{$state->t_gene_ids()} if $state->t_gene_ids();

    #
    # If gene IDs passed in via the query string from an external
    # application, these take precedent. These are always assumed
    # to be Ensembl IDs (gene ID type = 0).
    #
    if ($q->param("id")) {
        foreach my $id ($q->param("id")) {
            push @in_t_gene_ids, $id;
        }
        #$in_t_gene_id_type = DFLT_GENE_ID_TYPE;
    }

    # this should be later modified to give the user a choice
    my $biotype = DFLT_BIOTYPE;
	#my $biotype = 'protein_coding' if $species eq 'worm';
	#$biotype = _parse_biotype($biotype);
    $state->biotype($biotype);

    my $total_genes;
    if ($biotype) {
        $total_genes = $self->fetch_gene_count("biotype = '$biotype'");
    } else {
        $total_genes = $self->fetch_gene_count();
    }
    my $db_info = $self->fetch_db_info();
    #my $xgid_types = $self->fetch_external_gene_id_types();
    my $cl_hash = $self->fetch_conservation_levels();
    my $thl_hash = $self->fetch_threshold_levels();
    my $srl_hash = $self->fetch_search_region_levels();

    # This doesn't work properly, db_info methods seem to be lost when 
    # retrieving from state later on DJA 19/10/2010
    #$state->db_info($db_info);

    my @cl_levels  = sort keys %$cl_hash;
    my @thl_levels = sort keys %$thl_hash;
    my @srl_levels = sort keys %$srl_hash;

    my $tax_group       = $db_info->tax_group();
    my $min_ic          = $db_info->min_ic();
    my $max_upstream_bp = $db_info->max_upstream_bp();

    # for worms, tax_group include nematodes,vertebrates,insects
    # you could give this option for insects as well
    # in the input menu, you need to be able to choose from 3 tax groups
    # easier to deal with if you separate them
    my @tax_groups = split /\s*,\s*/, $tax_group;
    my $num_tax_groups = scalar @tax_groups;

    $state->has_operon($db_info->has_operon());
    $state->tax_groups(\@tax_groups);
    $state->num_tax_groups($num_tax_groups);

    #
    # Connect to JASPAR DB and retrieve TF info
    #
    $self->jaspar_db_connect();

    my %core_tf_sets;
    my %pending_tf_sets;
    my %pbm_tf_sets;
    foreach my $tax_group (@tax_groups)
    {
        my $core_tf_set = $self->fetch_tf_set(
            -collection => 'CORE',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $core_tf_sets{$tax_group} = $core_tf_set;
        
        my $pending_tf_set = $self->fetch_tf_set(
            -collection => 'PENDING',
            -tax_group  => $tax_group,
            -min_ic     => $min_ic
        );
        $pending_tf_sets{$tax_group} = $pending_tf_set;
        
        unless ($species eq 'yeast') {
            my $pbm_tf_set = $self->fetch_tf_set(
                -collection => 'PBM',
                -tax_group  => $tax_group,
                -min_ic     => $min_ic
            );
            $pbm_tf_sets{$tax_group} = $pbm_tf_set;
        }
    }

    #my $fam_tf_set;
    #unless ($species eq 'yeast' or $species eq 'worm') {
    #    $fam_tf_set = $self->fetch_tf_set(-collection => 'FAM');
    #}
    
    # so that the second db access won't be necessary 
    # why do I get tainted data error?
    #$state->core_tf_sets(\%core_tf_sets);
    #$state->pbm_tf_sets(\%pbm_tf_sets);
    #$state->pending_tf_sets(\%pending_tf_sets);
    #$state->fam_tf_set($fam_tf_set);

    #printf STDERR "\ncore_tf_set:\n"
    #    . Data::Dumper::Dumper($core_tf_set) . "\n\n";

    #printf STDERR "input_gene_ids:\n" . Data::Dumper::Dumper(\@input_gene_ids);

    #
    # Format a tax group list for display on the web page. The $tax_group
    # variable is obtained from the db and may be a single value or a comma
    # separated list (no space after comma).
    #
    my @s_tax_groups;
    foreach my $tg (@tax_groups) {
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
        nresults                => NUM_RESULTS,
        dflt_nresults           => DFLT_NUM_RESULTS,
        zcutoffs                => ZSCORE_CUTOFFS,
        fcutoffs                => FISHER_CUTOFFS,
        dflt_zcutoff            => DFLT_ZSCORE_CUTOFF,
        dflt_fcutoff            => DFLT_FISHER_CUTOFF,
        cl_dflt_level           => DFLT_CONSERVATION_LEVEL,
        thl_dflt_level          => DFLT_THRESHOLD_LEVEL,
        srl_dflt_level          => DFLT_SEARCH_REGION_LEVEL,
        max_target_genes        => MAX_TARGET_GENES,
        dflt_bg_num_rand_genes  => DFLT_BG_NUM_RAND_GENES,
        sid                     => $state->sid(),
        species                 => $species,
        db_info                 => $db_info,
        total_genes             => $total_genes,
        #xgid_types              => $xgid_types,
        cl_hash                 => $cl_hash,
        thl_hash                => $thl_hash,
        srl_hash                => $srl_hash,
        cl_levels               => \@cl_levels,
        thl_levels              => \@thl_levels,
        srl_levels              => \@srl_levels,
        tax_groups              => \@tax_groups,
        num_tax_groups          => $num_tax_groups,
        tax_group_list          => $tax_group_list,
        core_tf_sets            => \%core_tf_sets,
        pbm_tf_sets             => \%pbm_tf_sets,
        pending_tf_sets         => \%pending_tf_sets,
        #fam_tf_set              => $fam_tf_set,
        #in_t_gene_id_type       => $in_t_gene_id_type,
        in_t_gene_ids           => \@in_t_gene_ids,
        var_template            => "input_gene_ssa.html"
    };

    my $output = $self->process_template('master.html', $vars);
    #print STDERR "input results:\n"
    #    . Data::Dumper::Dumper($output);

    return $output;
}


sub results
{
    my $self = shift;

    my $q = $self->query;

    my $state = $self->state();

    my $analysis_type = 'default';

    my $species = $state->species() || $self->param('species');

    if (!defined $species) {
        return $self->error("Species not specified");
    }

    my $opdba = $self->opdba();

    #
    # This doesn't work properly, db_info methods seem to get lost
    # DJA 19/10/2010
    #
    #my $db_info = $state->db_info();

    my $dbia = $opdba->get_DBInfoAdaptor();
    my $db_info = $dbia->fetch_db_info();

    #my $tax_group   = $db_info->tax_group();
    my $dflt_min_ic      = $db_info->min_ic();

	# need to check
	#my $tf_db = $q->param('tf_db');
	#my $collections_str = $q->param('collections');

    #
    # Fetch all the input form entries and store them in the state
    #

    my $t_id_input_method = $q->param('t_id_input_method');
    if (!$t_id_input_method) {
        return $self->error("Target dene ID input method not specified");
    }
    #$state->t_id_input_method($t_id_input_method);

    my $t_gene_id_type = $q->param("t_gene_id_type");
    #if (!defined $t_gene_id_type) {
    #    return $self->error("Target gene ID type not specified");
    #}
    #$state->t_gene_id_type($t_gene_id_type);

    my $bg_id_input_method = $q->param('bg_id_input_method');
    if (!$bg_id_input_method) {
        return $self->error("Background gene ID input method not specified");
    }
    #$state->bg_id_input_method($bg_id_input_method);

	my $bg_gene_id_type;
	my $bg_num_rand_genes;
	if ($bg_id_input_method eq 'paste' || $bg_id_input_method eq 'upload') {
	    $bg_gene_id_type = $q->param("bg_gene_id_type");
	    #unless (defined $bg_gene_id_type) {
	    #    return $self->error("Background gene ID type not specified");
	    #}
	} elsif ($bg_id_input_method eq 'random') {
        $bg_gene_id_type = DFLT_GENE_ID_TYPE;
        $bg_num_rand_genes = $q->param('bg_num_rand_genes');
	    unless ($bg_num_rand_genes) {
	        return $self->error(
                "Number of random background genes not specified"
            );
	    }
    }
	#$state->bg_gene_id_type($bg_gene_id_type);

	my $tf_select_method = $q->param('tf_select_method');
	if (!$tf_select_method) {
		return $self->error("TFBS profiles not specified");
	}

	#
	# Divine the collection name from the tf_select_method parameter.
	# The tf_select_method parameter should be something like, 
	# e.g. core_min_ic, fam_specific, pbm etc.
	#
	my $collection;
	if ($tf_select_method =~ /(\w+?)_/) {
		$collection = uc $1;
	} else {
		$collection = uc $tf_select_method;
	}
	#$state->collection($collection);

	my $tf_select_criteria;
	if ($tf_select_method =~ /min_ic/) {
		$tf_select_criteria = 'min_ic';
	} elsif ($tf_select_method =~ /specific/) {
		$tf_select_criteria = 'specific';
	} else {
		$tf_select_criteria = 'all';
	}
	#$state->tf_select_method($tf_select_method);

	my $min_ic;
	my @tax_groups;
	my $tax_group_str;
	my @tf_ids;
	my $tf_ids_str;

	if ($tf_select_criteria eq 'min_ic') {
		my $num_tax_groups = $state->num_tax_groups();
		if ($collection eq 'CORE') {
			$min_ic = $q->param('core_min_ic');
			if ($num_tax_groups > 1) {
				@tax_groups = $q->param("core_tax_groups");
				$tax_group_str = join(',', @tax_groups);
			} else {
				@tax_groups = @{$state->tax_groups()};
				$tax_group_str = $db_info->tax_group();
			}
		} elsif ($collection eq 'PBM') {
			$min_ic = $q->param('pbm_min_ic');
			if ($num_tax_groups > 1) {
				@tax_groups = $q->param("pbm_tax_groups");
				$tax_group_str = join(',', @tax_groups);
			} else {
				@tax_groups = @{$state->tax_groups()};
				$tax_group_str = $db_info->tax_group();
			}
		} elsif ($collection eq 'PENDING') {
			$min_ic = $q->param('pending_min_ic');
			if ($num_tax_groups > 1) {
				@tax_groups = $q->param("pending_tax_groups");
				$tax_group_str = join(',', @tax_groups);
			} else {
				@tax_groups = @{$state->tax_groups()};
				$tax_group_str = $db_info->tax_group();
			}
		}

		if (!$min_ic) {
			return $self->error("No TF minimum IC specified");
		} elsif ($min_ic < $dflt_min_ic) {
			return $self->error(
                "Minimum IC specified is less than $dflt_min_ic"
            );
		}
		
		if ($collection eq 'CORE' && scalar @tax_groups == 0) {
			return $self->error(
                "No JASPAR CORE collection tax groups selected"
            );
		}

	} elsif ($tf_select_criteria eq 'specific') {

		if ($collection eq 'CORE') {
			push @tf_ids, $q->param('core_tfs');
		} elsif ($collection eq 'PBM') {
			push @tf_ids, $q->param('pbm_tfs');
		} elsif ($collection eq 'PENDING') {
			push @tf_ids, $q->param('pending_tfs');
		}

		if (scalar @tf_ids == 0) {
			return $self->error("No specific TFs selected");
		}
		
		$tf_ids_str = join(',', @tf_ids);
	
	}

	#$state->tf_ids(\@tf_ids);
	#$state->min_ic($min_ic);
	#$state->tax_groups(\@tax_groups);

	my $conservation_level;
    if ($species eq	'yeast') {
		$conservation_level = 1;
	} else {
		$conservation_level = $q->param('conservation_level');
		if (!$conservation_level) {
			$conservation_level = DFLT_CONSERVATION_LEVEL;
		}
	}
	my $threshold_level     = $q->param('threshold_level');
	my $threshold           = $q->param('threshold');
	my $search_region_level = $q->param('search_region_level');
	my $upstream_bp         = $q->param('upstream_bp');
	my $downstream_bp       = $q->param('downstream_bp');

	#my $cl_hash = $self->fetch_conservation_levels();
	#my $thl_hash = $self->fetch_threshold_levels();
	#my $srl_hash = $self->fetch_search_region_levels();

	#my $min_conservation = $cl_hash->{$conservation_level}->min_conservation();

	#
	# No need to check whether this is default or custom analysis type here
	# Will be handled by the actual analysis script
	#

    #
    # If threshold, upstream/downstream bp are defined strip any blanks and
    # if result is empty string, set to undef as future code needs to
    # distinguish undefined values from empty strings. Note: values of 0 are OK.
    # Should really check that if these values are set, they are actually
    # numeric.
    #
    # DJA 2012/03/23
    #
    if (defined $upstream_bp) {
        $upstream_bp =~ s/\s+//g;

        if ($upstream_bp eq '') {
            $upstream_bp = undef;
        }
    }

    if (defined $downstream_bp) {
        $downstream_bp =~ s/\s+//g;

        if ($downstream_bp eq '') {
            $downstream_bp = undef;
        }
    }

    if (defined $threshold) {
        $threshold =~ s/\s+//g;

        if ($threshold eq '') {
            $threshold = undef;
        }
    }

    #
    # Convert to number in range 0-1. Rather than '%' as analysis script
    # will interpret correctly. This is also consistent with the way this
    # value would be retrieved from the threshold level hash.
    #
    # DJA 2012/03/07
    #
    if (defined $threshold) {
        $threshold /= 100;
    }
	
    #
    # If neither threshold or threshold level is set, then set threshold level
    # to default but DON'T set actual threshold as it messes up the logic
    # further down with regard to  which arguments are passed to analysis
    # script, e.g. if threshold is set then it is assumed we are doing a custom
    # analysis and threshold level is NOT passed to analysis script.
    #
    # NOTE: This should never happen as the threshold level will ALWAYS have
    # some value since it is pull down menu.
    #
    # DJA 2012/03/07
    #
	if (!defined $threshold && !$threshold_level) {
        $threshold_level = DFLT_THRESHOLD_LEVEL;
		#$threshold = $thl_hash->{$threshold_level}->threshold();
	}

    if ($species eq 'yeast') {
        if (!defined($upstream_bp) && !$search_region_level) {
			$search_region_level = DFLT_SEARCH_REGION_LEVEL;
        }
    } else {
        if (defined($upstream_bp) != defined($downstream_bp)) {
            return $self->error(
                "Both upstream and downstream search regions must be specified"
            );
        } elsif (!$upstream_bp && !$downstream_bp) {
            if (!$search_region_level) {
                $search_region_level = DFLT_SEARCH_REGION_LEVEL;
            }

            #
            # DON'T set upstream/downstream bp for the same reasons as outlined
            # for not setting the threshold when setting default threshold level
            # above.
            #
            #$upstream_bp = $srl_hash->{$search_region_level}->upstream_bp();
            #$downstream_bp = $srl_hash->{$search_region_level}->downstream_bp();
        }
	}

	my $num_display_results;
	my $zscore_cutoff;
	my $fisher_cutoff;
	my $result_type = $q->param('result_type');
	#$state->result_type($result_type);
	if ($result_type eq 'top_x_results') {
		$num_display_results = $q->param('num_display_results');
	} elsif ($result_type eq 'significant_hits') {
		$zscore_cutoff = $q->param('zscore_cutoff');
		$fisher_cutoff = $q->param('fisher_cutoff');
	}

	my $result_sort_by = $q->param('result_sort_by');
	#$state->result_sort_by($result_sort_by);

	my $user_t_gene_file;
	if ($t_id_input_method eq 'upload') {
		$user_t_gene_file = $q->param('t_gene_file');
		$user_t_gene_file =~ s/.*\///;
        $user_t_gene_file =~ s/\s+//g;
	}

	my $user_bg_gene_file;
	if ($bg_id_input_method eq 'upload') {
		$user_bg_gene_file = $q->param('bg_gene_file');
		$user_bg_gene_file = s/.*\///;
        $user_bg_gene_file =~ s/\s+//g;
	}

	#
	# Create a temporary working directory for all the temp. input files and
	# output results file as a sub-directory under the defined temp. dir.
	#
	my $tempdir = tempdir(DIR => ABS_HTDOCS_RESULTS_PATH);

	my $job_id = $tempdir;
	$job_id =~ s/.*\///;

	my $t_gene_filename = $self->get_t_gene_file(
		$t_id_input_method,
		$tempdir
	);

	return $self->error("Could not create target gene IDs file")
		if !$t_gene_filename;

	my $bg_gene_filename;
	if ($bg_id_input_method eq 'upload' || $bg_id_input_method eq 'paste') {
		$bg_gene_filename = $self->get_bg_gene_file(
			$bg_id_input_method,
			$tempdir
		);

		return $self->error("Could not create background gene IDs file")
			if !$bg_gene_filename;
	}

	#
	# Call the analysis script
	#
	my $command = OPOSSUM_SCRIPTS_PATH . "/opossum_gene_ssa.pl"
		. " -j $job_id"
		. " -tgf $t_gene_filename"
		. " -s $species"
		. " -d $tempdir"
		. " -web";

	if (defined $t_gene_id_type) {
		$command .= " -tgidt $t_gene_id_type";
	}

	if ($bg_id_input_method eq 'upload' || $bg_id_input_method eq 'paste') {
		$command .= " -bgf $bg_gene_filename";
	} elsif ($bg_id_input_method eq 'random') {
        $command .= " -bnr $bg_num_rand_genes";
	}

	if (defined $bg_gene_id_type) {
		$command .= " -bgidt $bg_gene_id_type";
	}

	if ($state->has_operon()) {
		$command .= " -has_operon";
	}

	#if (defined $tf_db) {
	#	$command .= " -tfdb $tf_db";
	#}

    if (defined $collection) {
        $command .= " -co $collection";
    }

	if (defined $tf_ids_str) {
		$command .= " -tfids $tf_ids_str";
	} else {
		if ($tax_group_str) {
			$command .= " -tax '$tax_group_str'";
		}
		if (defined $min_ic) {
			$command .= " -ic $min_ic";
		}
	}

	if (defined $conservation_level) {
		$command .= " -cl $conservation_level";
	}

    #
    # User defined threshold overrides threshold_level. DO NOT pass -th
    # argument to script if threshold is set as the threshold level is ALWAYS
    # set by the browser which would result in the analysis script ALWAYS
    # performing a default analysis.
    #
	if (defined $threshold) {
		$command .= " -th $threshold";
	} else {
		$command .= " -thl $threshold_level";
	}

    #
    # User defined search region overrides search_region_level. DO NOT pass
    # -srl argument to script if upstream/downstream bp are set for the same
    # reasons as not passing -th parameter above.
    #
	if (defined $upstream_bp) {
		$command .= " -up $upstream_bp";
    }

	if (defined $downstream_bp) {
		$command .= " -dn $downstream_bp";
    }

	if (!defined $upstream_bp && !defined $downstream_bp) {
		$command .= " -srl $search_region_level";
    }

	if ($result_type eq 'top_x_results') {
		$command .= " -n $num_display_results";
	} elsif ($result_type eq 'significant_hits') {
		$command .= " -zcutoff $zscore_cutoff";
		$command .= " -fcutoff $fisher_cutoff";
	}

	if (defined $result_sort_by) {
		$command .= " -sr $result_sort_by";
	}

	if (defined $user_t_gene_file) {
		$command .= " -utgf '$user_t_gene_file'";
	}

	if (defined $user_bg_gene_file) {
		$command .= " -ubgf '$user_bg_gene_file'";
	}

	if ($state->biotype()) {
		my $biotype = $state->biotype();
		if (ref $biotype eq 'ARRAY') {
			$command .= sprintf " -biotype '%s'", join ',', @$biotype;
		} else {
			$command .= sprintf " -biotype '%s'", $biotype;
		}
	}

	my $submitted_time = scalar localtime(time);
	printf LOG "\nStarting SSA analysis at $submitted_time:\n$command\n\n";

	#my $rv = system("$command >/dev/null 2>&1 &");
	my $out = `exec 2>&1; $command`;
	if ($out) {
		printf LOG "\nCommand returned $out\n";
        #return $self->error("$out");
	}

	my $outfile = REL_HTDOCS_RESULTS_PATH . "/$job_id/results.html";

# 	This does not work
#    my $output = $self->process_html(
#		ABS_HTDOCS_RESULTS_PATH . "/$job_id/results.html");

    my $redirect = qq{
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">";
        <html>
        <head>
        <title>Your Page Title</title>
        <meta http-equiv="REFRESH" content="0;url=$outfile"></HEAD>
        <BODY>
        </BODY>
        </HTML>
    };

	return $redirect;
}

sub initialize_state
{
	my ($self, $state) = @_;

	my $species = $self->param('species');

	unless ($species) {
		return $self->error("Species is undefined");
	}
	$state->species($species);

	my $heading = sprintf "%s Single Site Analysis",
		ucfirst $state->species();
	$state->heading($heading);

	$state->debug(DEBUG);
	$state->title("oPOSSUM $heading");
	$state->bg_color_class(BG_COLOR_CLASS);

	#$state->errors(undef);
	#$state->warnings(undef);
}

1;


