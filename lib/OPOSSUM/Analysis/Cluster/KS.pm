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

package OPOSSUM::Analysis::Cluster::KS;

use strict;

use Carp;
use File::Temp qw/tempfile tempdir/;

use OPOSSUM::Analysis::Cluster::Values;
use OPOSSUM::Analysis::Cluster::ValuesIO;
use OPOSSUM::Analysis::Cluster::KSResult;
use OPOSSUM::Analysis::Cluster::KSResultSet;
use OPOSSUM::Analysis::Cluster::KSResultsIO;

=head2 new

 Title    : new
 Usage    : $ks = OPOSSUM::Analysis::Cluster::KS->new();
 Function : Create a new OPOSSUM::Analysis::Cluster::KS object.
 Returns  : An OPOSSUM::Analysis::Cluster::KS object.
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
 Args     : bg_vals - OPOSSUM::Analysis::Cluster::Values object
            t_vals  - OPOSSUM::Analysis::Cluster::Values object 

=cut

sub calculate_KS_with_bg_values
{
    my ($self, $bg_vals, $t_vals) = @_;
	
    return if !$t_vals || !$bg_vals;

	if (!$bg_vals->isa("OPOSSUM::Analysis::Cluster::Values")) {
        carp "bg vals is not an OPOSSUM::Analysis::Cluster::Values object";
        return;
    }
	
	if (!$t_vals->isa("OPOSSUM::Analysis::Cluster::Values")) {
        carp "test vals is not an OPOSSUM::Analysis::Cluster::Values object";
        return;
    }

    #
    # Check that the background and test counts are in sync with each
    # other, i.e. that they have the same number of TF IDs.
    #
    my $t_cluster_ids  = $t_vals->get_all_cluster_ids();
	my $bg_cluster_ids = $bg_vals->get_all_cluster_ids();
	
    if (scalar @$bg_cluster_ids != scalar @$t_cluster_ids) {
        carp "background and test counts do not have same number of TF clusters"
            . " defined";
        return;
    }

    my $cluster_idx = 0;
    while ($cluster_idx < scalar @$bg_cluster_ids) {
        if ($bg_cluster_ids->[$cluster_idx] ne $t_cluster_ids->[$cluster_idx]) {
            carp "background and test values for TFCluster ID "
                . $cluster_idx + 1
                . " differ: "
                . $bg_cluster_ids->[$cluster_idx] . " "
                . $t_cluster_ids->[$cluster_idx];
            return;
        }
        $cluster_idx++;
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
	print OUTPUT "TFCluster\tp-value\tBG_distribution\n";
	close(OUTPUT);

	# probably not necessary unless needed for debugging
	my $log_file = File::Temp::tempnam($self->{_tempdir}, "");
	if (!$log_file) {
		carp "error generating temp. log file name";
		return;
	}
	
	#
	# go through each TFCluster and perform the KS test
	# output file: p-value for each TFCluster
    #
	$cluster_idx = 0;
	while ($cluster_idx < scalar @$t_cluster_ids)
	{
        my $cluster_id = $t_cluster_ids->[$cluster_idx];
        
		if (scalar @{$t_vals->all_cluster_values($cluster_id)} <= 0 ||
			scalar @{$bg_vals->all_cluster_values($cluster_id)} <= 0) {
			print "Skipping $cluster_id\n";
			$cluster_idx++;
			next;
		}
		#print "KS for $cluster_id\n";
		
		my $R_file = File::Temp::tempnam($self->{_tempdir}, "");
		#my $R_file = "R_file_$cluster_id";
		if (!$R_file) {
			carp "error generating temp. R file name";
			return;
		}
		
		my $bg_vals_file = $self->_write_vals($bg_vals, $cluster_id, 'bg');
		if (!$bg_vals_file) {
			carp "error creating temporary bg vals file";
			$self->_clean_temp_files;
			return;
		}
		
		my $t_vals_file = $self->_write_vals($t_vals, $cluster_id, 't');
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
			pvals <- data.frame(TFCluster="$cluster_id", p=-log(result\$p.value), bg="data")\n
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

        $cluster_idx++;
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
 Returns  : An OPOSSUM::Analysis::Cluster::KSResultSet object.
 Args     : bg_distribution - string code for distribution as recognized by
								R's ks.test function
            t_vals       - OPOSSUM::Analysis::Cluster::Values object 

=cut

sub calculate_KS_with_bg_distribution
{
    my ($self, $bg_distribution, $t_vals) = @_; 
	
    return if !$bg_distribution || !$t_vals;
	
	if (!$t_vals->isa("OPOSSUM::Analysis::Cluster::Values")) {
        carp "test vals is not an OPOSSUM::Analysis::Cluster::Values object";
        return;
    }
	
    my $t_cluster_ids  = $t_vals->get_all_cluster_ids();

    my $cluster_idx = 0;
	
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
	print OUTPUT "TFCluster\tp-value\tBG_distribution\n";
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
	while ($cluster_idx < scalar @$t_cluster_ids)
	{
        my $cluster_id = $t_cluster_ids->[$cluster_idx];
        
		if (scalar @{$t_vals->all_cluster_values($cluster_id)} <= 0) {
			print "Skipping $cluster_id\n";
			$cluster_idx++;
			next;
		}
		print "KS for $cluster_id\n";
		
		my $R_file = File::Temp::tempnam($self->{_tempdir}, "");
		#my $R_file = "R_file_$cluster_id";
		if (!$R_file) {
			carp "error generating temp. R file name";
			return;
		}
	
		my $t_vals_file = $self->_write_vals($t_vals, $cluster_id);
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
			pvals <- data.frame(TFCluster="$cluster_id", p=-log(result\$p.value), bg="$bg_distribution")\n
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

        $cluster_idx++;
	}

	my $results = _read_results($output_file);	
	$self->_clean_temp_files;

    return $results;
}

#
# Ouput the OPOSSUM::Analysis::Cluster::Values to a temporary file for further processing
#
sub _write_vals
{
    my ($self, $vals, $cluster_id, $header) = @_;

    my $vals_file = File::Temp::tempnam($self->{_tempdir}, "");
	#$header = 'a' if !$header;
	#my $vals_file = $header . "_values_$cluster_id.txt";

    if (!$vals_file) {
        carp "error generating temp. vals file name";
        return;
    }

    my $valsIO = OPOSSUM::Analysis::Cluster::ValuesIO->new(
        -file   => ">$vals_file",
        -format => 'ks'
    );
    if (!$valsIO) {
        carp "error creating temp. vals file $vals_file";
        return;
    }

    if (!$valsIO->write_vals($vals, $cluster_id)) {
        carp "error writing vals to temp. file $vals_file";
        return;
    }

    return $vals_file;
}

#
# Read the results of the analysis from the given file and return them as
# an OPOSSUM::Analysis::Cluster::KSResultSet object
#
sub _read_results
{
    my ($file) = @_;

    my $ksrIO = OPOSSUM::Analysis::Cluster::KSResultsIO->new(-file => $file);
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
