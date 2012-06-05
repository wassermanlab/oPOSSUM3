=head1 NAME

OPOSSUM::KS.pm - module to perform Kolmogorov-Smirnov test

=head1 SYNOPSIS

 $fisher = OPOSSUM::Analysis::KS->new();

 $results = $ks->calculate_KS_probability(
					    $background_vals,
					    $target_vals);

=head1 DESCRIPTION

Given the values in the target set, calculate the Kolmogorov-Smirnov test against
the uniform distribution and output the p-value. To be used for comparing the
distributions of distances between predicted TFBSs and ChIP-Seq peaks, but could
be used for other values. Must have R installed to run.

=head1 AUTHORS

 Andrew Kwon
 Modified from original code by David Arenillas, Shannan Ho Sui

 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: tjkwon@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::KS;

use strict;

use Carp;
use File::Temp qw/tempfile tempdir/;

use OPOSSUM::Analysis::Values;
use OPOSSUM::Analysis::ValuesIO;
use OPOSSUM::Analysis::KSResult;
use OPOSSUM::Analysis::KSResultSet;
use OPOSSUM::Analysis::KSResultsIO;

=head2 new

 Title    : new
 Usage    : $ks = OPOSSUM::Analysis::KS->new();
 Function : Create a new OPOSSUM::Analysis::KS object.
 Returns  : An OPOSSUM::Analysis::KS object.
 Args     : None

=cut

sub new
{
    my ($class, @args) = @_;

    # Create temp. working directory
    my $tempdir = tempdir("KSXXXXXX", TMPDIR => 1, CLEANUP => 1);
    if (!$tempdir) {
        carp "error creating temporary working directory - $!";
        return;
    }

    my $self = bless {}, ref $class || $class;

    $self->{_tempdir} = $tempdir;

    return $self;
}

=head2 calculate_KS_with_bg_values

 Title    : calculate_KS_with_bg_values
 Usage    : $results = $ks->calculate_KS_with_bg_values(
                $bg_vals,
                $t_vals
            );
 Function : Perform the Kolmogorov-Smirnov analysis on bg_vals and t_vals
 Returns  : An OPOSSUM::Analysis::KSResultSet object.
 Args     : bg_vals - OPOSSUM::Analysis::Values object
            t_vals  - OPOSSUM::Analysis::Values object 

=cut

sub calculate_KS_with_bg_values
{
    my ($self, $bg_vals, $t_vals) = @_;
	
    return if !$t_vals || !$bg_vals;

	if (!$bg_vals->isa("OPOSSUM::Analysis::Values")) {
        carp "bg vals is not an OPOSSUM::Analysis::Values object";
        return;
    }
	
	if (!$t_vals->isa("OPOSSUM::Analysis::Values")) {
        carp "test vals is not an OPOSSUM::Analysis::Values object";
        return;
    }

    #
    # Check that the background and test counts are in sync with each
    # other, i.e. that they have the same number of TF IDs.
    #
    my $t_tf_ids  = $t_vals->tf_ids();
	my $bg_tf_ids = $bg_vals->tf_ids();
	
    if (scalar @$bg_tf_ids != scalar @$t_tf_ids) {
        carp "background and test counts do not have same number of TFs"
            . " defined";
        return;
    }

    my $tf_idx = 0;
    while ($tf_idx < scalar @$bg_tf_ids) {
        if ($bg_tf_ids->[$tf_idx] ne $t_tf_ids->[$tf_idx]) {
            carp "background and test counts TF ID "
                . $tf_idx + 1
                . " differ: "
                . $bg_tf_ids->[$tf_idx] . " "
                . $t_tf_ids->[$tf_idx];
            return;
        }
        $tf_idx++;
    }
	
	my $output_file = File::Temp::tempnam($self->{_tempdir}, "");
	#my $output_file = "output_file";
	if (!$output_file) {
		carp "error generating temp. output file name";
		return;
	}
	
	if (!open(OUTPUT, ">$output_file")) {
		carp "error writing to temp. output";
		return;
	}
	print OUTPUT "TF\tp-value\tBG_distribution\n";
	close(OUTPUT);

	# probably not necessary unless needed for debugging
	my $log_file = File::Temp::tempnam($self->{_tempdir}, "");
	if (!$log_file) {
		carp "error generating temp. log file name";
		return;
	}
	
	#
	# go through each TF and perform the KS test
	# output file: p-value for each TF
    #
	$tf_idx = 0;
	while ($tf_idx < scalar @$t_tf_ids)
	{
        my $tf_id = $t_tf_ids->[$tf_idx];
        
		if (scalar @{$t_vals->all_tfbs_values($tf_id)} <= 0 ||
			scalar @{$bg_vals->all_tfbs_values($tf_id)} <= 0) {
			print "Skipping $tf_id\n";
			$tf_idx++;
			next;
		}
		#print "KS for $tf_id\n";
		
		my $R_file = File::Temp::tempnam($self->{_tempdir}, "");
		#my $R_file = "R_file_$tf_id";
		if (!$R_file) {
			carp "error generating temp. R file name";
			return;
		}
		
		my $bg_vals_file = $self->_write_vals($bg_vals, $tf_id, 'bg');
		if (!$bg_vals_file) {
			carp "error creating temporary bg vals file";
			$self->_clean_temp_files;
			return;
		}
		
		my $t_vals_file = $self->_write_vals($t_vals, $tf_id, 't');
		if (!$t_vals_file) {
			carp "error creating temporary target vals file";
			$self->_clean_temp_files;
			return;
		}
	
		if (!open(RFILE, ">$R_file")) {
			carp "error opening temporary R file $R_file";
			$self->_clean_temp_files;
			return;
		}

		my $rcmds = qq{
			t <- read.table("$t_vals_file", header=FALSE)\n
			bg <- read.table("$bg_vals_file", header=FALSE)\n
			result = ks.test(t[,2], bg[,2])\n
			pvals <- data.frame(TF="$tf_id", p=-log(result\$p.value), bg="data")\n
			write.table(pvals, "$output_file", quote=FALSE,
				row.names=FALSE, col.names=FALSE, append=TRUE, sep="\\t")\n
		};
	
		print RFILE "$rcmds\n";
	
		close(RFILE);

		my $R_cmd = "/usr/local/bin/R --no-save < $R_file 2>&1 1> $log_file";
		my $R_err = `$R_cmd`;
		if ($? >> 8) {
			carp "error running R command:\n\t$R_cmd\n$R_err";
			return;
		}

        $tf_idx++;
	}

	my $results = _read_results($output_file);	
	$self->_clean_temp_files;

    return $results;
}

=head2 calculate_KS_with_bg_distribution

 Title    : calculate_KS_with_bg_distribution
 Usage    : $results = $ks->calculate_KS_with_bg_distribution(
                $bg_distribution,
                $t_vals
            );
 Function : Perform the Kolmogorov-Smirnov analysis on t_vals, based on the specified
			background distributiont
 Returns  : An OPOSSUM::Analysis::KSResultSet object.
 Args     : bg_distribution - string code for distribution as recognized by
								R's ks.test function
            t_vals       - OPOSSUM::Analysis::Values object 

=cut

sub calculate_KS_with_bg_distribution
{
    my ($self, $bg_distribution, $t_vals) = @_; 
	
    return if !$bg_distribution || !$t_vals;
	
	if (!$t_vals->isa("OPOSSUM::Analysis::Values")) {
        carp "test vals is not an OPOSSUM::Analysis::Values object";
        return;
    }
	
    my $t_tf_ids  = $t_vals->tf_ids();

    my $tf_idx = 0;
	
	my $output_file = File::Temp::tempnam($self->{_tempdir}, "");
	#my $output_file = "output_file";
	if (!$output_file) {
		carp "error generating temp. output file name";
		return;
	}
	
	if (!open(OUTPUT, ">$output_file")) {
		carp "error writing to temp. output";
		return;
	}
	print OUTPUT "TF\tp-value\tBG_distribution\n";
	close(OUTPUT);

	# probably not necessary unless needed for debugging
	my $log_file = File::Temp::tempnam($self->{_tempdir}, "");
	if (!$log_file) {
		carp "error generating temp. log file name";
		return;
	}
	
	#
	# go through each TF and perform the KS test
	# output file: p-value for each TF
    #
	while ($tf_idx < scalar @$t_tf_ids)
	{
        my $tf_id = $t_tf_ids->[$tf_idx];
        
		if (scalar @{$t_vals->all_tfbs_values($tf_id)} <= 0) {
			print "Skipping $tf_id\n";
			$tf_idx++;
			next;
		}
		print "KS for $tf_id\n";
		
		my $R_file = File::Temp::tempnam($self->{_tempdir}, "");
		#my $R_file = "R_file_$tf_id";
		if (!$R_file) {
			carp "error generating temp. R file name";
			return;
		}
	
		my $t_vals_file = $self->_write_vals($t_vals, $tf_id);
		if (!$t_vals_file) {
			carp "error creating temporary target vals file";
			$self->_clean_temp_files;
			return;
		}
	
		if (!open(RFILE, ">$R_file")) {
			carp "error opening temporary R file $R_file";
			$self->_clean_temp_files;
			return;
		}

		my $rcmds = qq{
			t <- read.table("$t_vals_file", header=FALSE)\n
			result = ks.test(t[,2], '$bg_distribution')\n
			pvals <- data.frame(TF="$tf_id", p=-log(result\$p.value), bg="$bg_distribution")\n
			write.table(pvals, "$output_file", quote=FALSE,
				row.names=FALSE, col.names=FALSE, append=TRUE, sep="\\t")\n
		};
	
		print RFILE "$rcmds\n";
	
		close(RFILE);

		my $R_cmd = "/usr/local/bin/R --no-save < $R_file 2>&1 1> $log_file";
		my $R_err = `$R_cmd`;
		if ($? >> 8) {
			carp "error running R command:\n\t$R_cmd\n$R_err";
			return;
		}

        $tf_idx++;
	}

	my $results = _read_results($output_file);	
	$self->_clean_temp_files;

    return $results;
}

#
# Ouput the OPOSSUM::Analysis::Values to a temporary file for further processing
#
sub _write_vals
{
    my ($self, $vals, $tf_id, $header) = @_;

    my $vals_file = File::Temp::tempnam($self->{_tempdir}, "");
	#$header = 'a' if !$header;
	#my $vals_file = $header . "_values_$tf_id.txt";

    if (!$vals_file) {
        carp "error generating temp. vals file name";
        return;
    }

    my $valsIO = OPOSSUM::Analysis::ValuesIO->new(
        -file   => ">$vals_file",
        -format => 'ks'
    );
    if (!$valsIO) {
        carp "error creating temp. vals file $vals_file";
        return;
    }

    if (!$valsIO->write_vals($vals, $tf_id)) {
        carp "error writing vals to temp. file $vals_file";
        return;
    }

    return $vals_file;
}

#
# Read the results of the analysis from the given file and return them as
# an OPOSSUM::Analysis::KSResultSet object
#
sub _read_results
{
    my ($file) = @_;

    my $ksrIO = OPOSSUM::Analysis::KSResultsIO->new(-file => $file);
    if (!$ksrIO) {
        carp "error opening results file $file";
        return;
    }

    my $results = $ksrIO->read_results();

    $ksrIO->close;

    return $results;
}

#
# Clean up any temporary files which were created during the analysis.
#
sub _clean_temp_files
{
    my $self = shift;

   #
   # Explicitely remove temp. files. Temp. dir is set to autoclean anyway
   # upon object destruction but clean the temp. files to minimize chance
   # of file name clashes if calculate_KS_probability is called repeatedly
   # for the same OPOSSUM::Analysis::KS object.
   #
    my $tempdir = $self->{_tempdir};
    if ($tempdir && -d $tempdir) {
        unlink(<$tempdir/*>);
    }
}

1;
