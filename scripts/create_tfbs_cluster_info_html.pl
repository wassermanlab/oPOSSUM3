#!/usr/local/bin/perl -w

=head1 NAME

create_tfbs_cluster_info_html.pl

=head1 SYNOPSIS

  opossum_seq_tca.pl
                  -d results_dir
                 [-tdb tf_database]
                 [-cdb tf_cluster_database]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -d results_dir   = Name of directory used for output html files
   -tdb tf_database = Specify which TFBS database to use
                      (default = JASPAR_2010)
   -cdb tf_cluster_database = Specify which TFBSCluster database to use
                      (default = oPOSSUM_cluster)

=head1 DESCRIPTION

Given the TFBS cluster database, create info html pages for each cluster,
with detailed list of the member TFs and their attributes
These pages will be used by oPOSSUM TCA.

=head1 AUTHOR

  Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca

=cut

use strict;

use lib '/space/devel/oPOSSUM3/cgi-bin';
use oPossumWebOpt;
use oPossumTCAWebOpt;

use lib OPOSSUM_LIB_PATH;
use lib TFBS_CLUSTER_LIB_PATH;

use Getopt::Long;
use Pod::Usage;
use Carp;

use Template;
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

use TFBS::DB::JASPAR5;
use TFBSCluster::DBSQL::DBAdaptor;

use constant DEBUG			=> 0;
use constant HEADING        => 'TFBS Cluster Information';
use constant TITLE          => 'TFBS Cluster Information';
use constant BG_COLOR_CLASS => 'bgc_seq_tca';

my $results_dir;
my $tf_db;
my $cl_db;
my $log_file;
GetOptions(
    'd=s'       => \$results_dir,
    'tdb=s'     => \$tf_db,
    'cdb=s'     => \$cl_db,
	'l=s'		=> \$log_file
);

die "No results directory specified\n" if !$results_dir;

#
# Initialize logging
#
if (!$log_file) {
	$log_file = "create_tfbs_cluster_info_html.log";
}

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

$tf_db = JASPAR_DB_NAME if !$tf_db;
$cl_db = TFBS_CLUSTER_DB_NAME if !$cl_db;

#
# Connect to JASPAR and TFBS cluster databases
#
my $jdb = TFBS::DB::JASPAR5->connect(
    "dbi:mysql:" . $tf_db . ":" . JASPAR_DB_HOST,
    JASPAR_DB_USER,
    JASPAR_DB_PASS
);
$logger->error("Could not connect to JASPAR database $tf_db") if !$jdb;

my $cdb = TFBSCluster::DBSQL::DBAdaptor->new(
    -host     => TFBS_CLUSTER_DB_HOST,
    -dbname   => $cl_db,
    -user     => TFBS_CLUSTER_DB_USER,
    -password => TFBS_CLUSTER_DB_PASS
);
$logger->error("Could not connect to TFBS cluster database $cl_db") if !$cdb;

my %get_matrix_args = (
    -matrixtype => 'PFM'
);

#
# retrieve the TFBS cluster set
#
my $tf_clusters = fetch_tf_clusters(%get_matrix_args);

foreach my $cluster (@$tf_clusters)
{
	write_tfbs_cluster_info_html($cluster);
}

sub fetch_tf_clusters
{
    my (%args) = @_;

    my %matrix_args = %args;

    unless ($matrix_args{-matrixtype}) {
        $matrix_args{-matrixtype} = 'PFM';
    }

    #printf STDERR "fetch_tf_cluster_set: matrix_args = \n"
    #    . Data::Dumper::Dumper(%matrix_args) . "\n";

    unless ($cdb) {
        fatal("No oPOSSUM_cluster DB found");
    }
    
    my $tfca = $cdb->get_TFClusterAdaptor;

    my $clusters;
    #print STDERR "Getting cluster set\n";
    if ($matrix_args{-family}) {
        $clusters = $tfca->fetch_by_tf_families($matrix_args{-family});
    } else {
        $clusters = $tfca->fetch_all();    
    }
	
    #print STDERR "# clusters = " . scalar(@$clusters) . "\n";

    $logger->error("Could not fetch TFBS clusters\n")
        if !$clusters || scalar(@$clusters) == 0;

    return $clusters;
}

sub write_tfbs_cluster_info_html
{
    my $tfc = shift;
	
	my $cid = $tfc->id;
	
	my @tfc_tfs;
	my %collections;
	my %tax_groups;
	my %tf_ic;
	foreach my $tfid (@{$tfc->tf_ids})
	{
		my $tf = $jdb->get_Matrix_by_ID($tfid);
		push @tfc_tfs, $tf;
		$collections{$tf->ID} = $tf->tag('collection');
		$tax_groups{$tf->ID} = $tf->tag('tax_group');
		$tf_ic{$tf->ID} = sprintf "%.2f", $tf->to_ICM->total_ic;
	}
	
	my $filename = $results_dir . "/c$cid" . "_info.html";
	open (FH, ">$filename") || fatal(
		"Could not create TFBS cluster info HTML file $filename"
	);
	
    my $vars = {
        abs_htdocs_path    => ABS_HTDOCS_PATH,
        rel_htdocs_path    => REL_HTDOCS_PATH,
        abs_cgi_bin_path   => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path   => REL_CGI_BIN_PATH,
        bg_color_class     => BG_COLOR_CLASS,
        title              => TITLE,
        heading            => HEADING,
        #section            => 'TFBS Cluster Information',
        version            => VERSION,
        devel_version      => DEVEL_VERSION,
		jaspar_url         => JASPAR_URL,
		tf_cluster         => $tfc,
		cluster_tfs        => \@tfc_tfs,
		collections        => \%collections,
		tax_groups         => \%tax_groups,
		tf_ic              => \%tf_ic,
        var_template       => "tfbs_cluster_info.html"
    };

    my $output = process_template('master.html', $vars);

    print FH $output;
	
	close (FH);
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
        || $logger->error(
            "Error processing template $input\n" . $template->error() . "\n\n"
        );

    return $string;
}


1;