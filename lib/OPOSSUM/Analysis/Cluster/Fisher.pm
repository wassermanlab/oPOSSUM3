=head1 NAME

OPOSSUM::Fisher.pm - module to perform Fisher's exact probability analysis

=head1 SYNOPSIS

 $fisher = OPOSSUM::Analysis::Cluster::Fisher->new();

 $results = $fisher->calculate_Fisher_probability(
					    $background_counts,
					    $target_counts);

=head1 DESCRIPTION

Given TFBS cluster counts for each gene in a target and
background set, calculate the one-tailed Fisher Exact Probability that each
binding site is overrepresented in the set of target genes.  Must have R
installed to run.

=head1 AUTHORS

 Original coding: Shannan Ho Sui

 Modified by: David Arenillas

 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::Analysis::Cluster::Fisher;

use strict;

use Carp;
use File::Temp qw/tempfile tempdir/;
use OPOSSUM::Analysis::Cluster::Counts;
use OPOSSUM::Analysis::Cluster::CountsIO;
use OPOSSUM::Analysis::Cluster::FisherResult;
use OPOSSUM::Analysis::Cluster::FisherResultSet;
use OPOSSUM::Analysis::Cluster::FisherResultsIO;

=head2 new

 Title    : new
 Usage    : $fisher = OPOSSUM::Analysis::Cluster::Fisher->new();
 Function : Create a new OPOSSUM::Analysis::Cluster::Fisher object.
 Returns  : An OPOSSUM::Analysis::Cluster::Fisher object.
 Args     : None

=cut

sub new
{
    my ($class, @args) = @_;

    # Create temp. working directory
    my $tempdir = tempdir("FisherXXXXXX", TMPDIR => 1, CLEANUP => 1);
    if (!$tempdir) {
        carp "error creating temporary working directory - $!";
        return;
    }

    my $self = bless {}, ref $class || $class;

    $self->{_tempdir} = $tempdir;

    return $self;
}

=head2 calculate_Fisher_probability

 Title    : calculate_Fisher_probability
 Usage    : $results = $fisher->calculate_Fisher_probability(
                $background_counts,
                $target_counts
            );
 Function : Perform the Fisher exact probability analysis to determine
            over-representation of TF profile clusters and return the results.
 Returns  : An OPOSSUM::Analysis::Cluster::FisherResultSet object.
 Args     : background_counts   - An OPOSSUM::Analysis::Cluster::Counts object
                                  containing the TFBS cluster counts in the
                                  background set of genes
            target_counts       - An OPOSSUM::Analysis::Cluster::Counts object
                                  containing the TFBS cluster counts in the
                                  target set of genes

=cut

sub calculate_Fisher_probability
{
    my ($self, $bg_counts, $t_counts) = @_;

    return if !$bg_counts || !$t_counts;
    if (!$bg_counts->isa("OPOSSUM::Analysis::Cluster::Counts")) {
        carp "background counts is not an OPOSSUM::Analysis::Cluster::Counts object";
        return;
    }

    if (!$t_counts->isa("OPOSSUM::Analysis::Cluster::Counts")) {
        carp "test counts is not an OPOSSUM::Analysis::Cluster::Counts object";
        return;
    }

    #
    # Check that the background and test counts are in sync with each
    # other, i.e. that they have the same number of TF IDs.
    #
    my $bg_cluster_ids = $bg_counts->get_all_cluster_ids;

    my $t_cluster_ids = $bg_counts->get_all_cluster_ids;

    if (scalar @$bg_cluster_ids != scalar @$t_cluster_ids) {
        carp "background and test counts do not have same number of TF clusters"
            . " defined";
        return;
    }

    my $cluster_idx = 0;
    while ($cluster_idx < scalar @$bg_cluster_ids) {
        if ($bg_cluster_ids->[$cluster_idx] ne $t_cluster_ids->[$cluster_idx]) {
            carp "background and test counts TFCluster ID "
                . $cluster_idx + 1
                . " differ: "
                . $bg_cluster_ids->[$cluster_idx] . " "
                . $t_cluster_ids->[$cluster_idx];
            return;
        }
        $cluster_idx++;
    }

    my $R_file = File::Temp::tempnam($self->{_tempdir}, "");
    if (!$R_file) {
        carp "error generating temp. R file name";
        return;
    }

    my $output_file = File::Temp::tempnam($self->{_tempdir}, "");
    if (!$output_file) {
        carp "error generating temp. output file name";
        return;
    }

    # probably not necessary unless needed for debugging
    my $log_file = File::Temp::tempnam($self->{_tempdir}, "");
    if (!$log_file) {
        carp "error generating temp. log file name";
        return;
    }

	# convert the count objects to cluster counts
	#my $bg_cluster_counts = convert_to_cluster_count($bg_counts);
	#my $t_cluster_counts = convert_to_cluster_count($t_counts);
	
    my $bg_counts_file = $self->_write_counts($bg_counts);
    if (!$bg_counts_file) {
        carp "error creating temporary background counts file";
        return;
    }

    my $t_counts_file = $self->_write_counts($t_counts);
    if (!$t_counts_file) {
        carp "error creating temporary target counts file";
        $self->_clean_temp_files;
        return;
    }

    if (!open(RFILE, ">$R_file")) {
        carp "error opening temporary R file $R_file";
        $self->_clean_temp_files;
        return;
    }

    my $rcmds = qq{
        bg <- scan("$bg_counts_file", list(tf="",hits=0,nohits=0),
            flush=TRUE)\n
        t <- scan("$t_counts_file", list(tf="",hits=0,nohits=0),
            flush=TRUE)\n
        results = c()\n
        for (i in 1:length(bg\$tf)) {\n
        testmat <- matrix(c(t\$nohits[i], t\$hits[i],
            bg\$nohits[i], bg\$hits[i]), nr = 2,
            dimnames = list(c("Site", "NoSite"),
            c("Background", "Target")))\n
        fisher = fisher.test(testmat, alternative="l")\n
        fscore = -log(fisher\$p.value)\n
        results = rbind(results, cbind(bg\$tf[i], t\$hits[i],
            t\$nohits[i], bg\$hits[i],
            bg\$nohits[i], fscore))\n
        }\n
        scores = c()\n
        index <- sort(as.numeric(results[,6]), index=TRUE)\$ix\n
        for(i in index) {\n
            scores <- rbind(scores, results[i,])\n
        }\n
        tablenames = c("TFCluster", "tHits", "tNoHits", "bHits",
                       "bNoHits", "p-value")\n
        scores2 <- data.frame(scores)\n
        names(scores2) = tablenames\n
        write.table(scores2, "$output_file", quote=FALSE,
            row.names=FALSE, sep="\\t")\n
        q()\n
    };
#            bg\$nohits[i], fisher\$p.value))\n

    print RFILE "$rcmds\n";

    close(RFILE);

    my $R_cmd = "/usr/local/bin/R --no-save < $R_file 2>&1 1> $log_file";
    my $R_err = `$R_cmd`;
    if ($? >> 8) {
        carp "error running R command:\n\t$R_cmd\n$R_err";
        return;
    }

    my $results = _read_results($output_file);

    #$self->_clean_temp_files;

    return $results;
}


#
# Ouput the OPOSSUM::Analysis::Cluster::Counts to a temporary file for further processing
#
sub _write_counts
{
    my ($self, $counts) = @_;

    my $counts_file = File::Temp::tempnam($self->{_tempdir}, "");
	#my $counts_file = "/tmp/Fisher_counts_test" . $counts->num_genes;

    if (!$counts_file) {
        carp "error generating temp. counts file name";
        return;
    }

    my $countsIO = OPOSSUM::Analysis::Cluster::CountsIO->new(
        -file   => ">$counts_file",
        -format => 'fisher'
    );
    if (!$countsIO) {
        carp "error creating temp. counts file $counts_file";
        return;
    }

    if (!$countsIO->write_counts($counts)) {
        carp "error writing counts to temp. file $counts_file";
        return;
    }

    return $counts_file;
}

#
# Read the results of the analysis from the given file and return them as
# an OPOSSUM::Analysis::Cluster::FisherResultSet object
#
sub _read_results
{
    my ($file) = @_;

    my $frIO = OPOSSUM::Analysis::Cluster::FisherResultsIO->new(-file => $file);
    if (!$frIO) {
        carp "error opening results file $file";
        return;
    }

    my $results = $frIO->read_results();

    $frIO->close;

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
   # of file name clashes if calculate_Fisher_probability is called repeatedly
   # for the same OPOSSUM::Analysis::Cluster::Fisher object.
   #
    my $tempdir = $self->{_tempdir};
    if ($tempdir && -d $tempdir) {
        unlink(<$tempdir/*>);
    }
}

1;
