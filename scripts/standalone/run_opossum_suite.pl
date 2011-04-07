#!/usr/local/bin/perl

=head1 NAME

  run_opossum_suite.pl

=head1 SYNOPSIS

  opossum_gene_acsa.pl
        -s species
        -g t_gene_file
        -id tf_id
        [-cl conservation_level]
        [-dist site_distance]
        [-b bg_gene_file]
        [-gt t_gene_id_type]
        [-bt bg_gene_id_type]
        [-bnr num_rand_bg_genes]
        [-has_operon]
        [-biotype biotype]
        [-db tf_database]
        [-co collection]
        [-tf tf_file]
        [-tax tax_groups]
        [-ic min_ic]
        [-th threshold]
        [-up upstream_bp]
        [-dn downstream_bp]
        [-n num_results | -zcutoff cutoff -fcutoff cutoff]
        [-sr sort_by]
        [-d results_dir]
        [-hits]

=head1 ARGUMENTS

Argument switches may be abbreviated where unique. Arguments enclosed by
brackets [] are optional.

   -s species           = Species
   -id tf_id            = Anchor TF ID
   -cl conservation_level = Conservation level at which to perform analysis
   -dist site_distance  = Maximum inter-binding site distance
   -g t_gene_file       = File containing list of target gene IDs
   -b bg_gene_file      = File containing background gene IDs
   -gt t_gene_id_type   = Target gene ID type
   -bt bg ID type       = Background gene ID type
   -bnr num rand bg genes
                        = Number of random background genes
   -has_operon          = Boolean indicating whether this species has
                          operons
   -biotype             = String containing biotype name(s) of genes to use
                          as background
   -db tf_database      = Specify which TF database to use
                         (default = JASPAR_2010)
   -co collections      = Specify which JASPAR collections to use
                         (default = CORE)
   -tf tf_file          = File containing list of TF IDs
   -tax tax_groups      = Specify a comma separated string of tax groups
   -ic min_ic           = Specify minimum IC of TFs
   -th threshold        = Minimum relative TFBS position weight matrix
                          (PWM) score to report in the analysis. The
                          thresold may be spefifies as a percentage, i.e.
                          '80%' or a decimal, i.e. 0.8.
                          Default = 80%
                          Min. = 75%
   -up upstream         = Amount of sequence upstream of TSS to include in
                          analysis
   -dn downstream       = Amount of sequence downstream of TSS to include in
                          analysis
   -n num_results       = Number of results to display
   -zcutoff cutoff      = Z-score cutoff of results to display
   -fcutoff cutoff      = Fisher p-value cutoff of results to display
   -sr sort_by          = Sort results by this value ('zscore', 'fisher')
   -d results_dir       = Name of directory used for input gene ID and TF ID
                          files and output results files
   -h hits              = save TFBS hits
   -l log_file          = log file
   
=head1 DESCRIPTION

Take a list of gene IDs, optional background gene IDs and optional subset
of transcription factors (TFs) either specified in an input file, or limited
by external (JASPAR) database name and information content or taxonomic
supergroup or all TFs in the oPOSSUM database. Also optionally specify PWM
score threshold.

Count the number of TFBSs for each TF which was found at the given
PWM score threshold for both the test and background set of genes. Perform
Fisher exact test and z-score analysis and output these results to the
output file. Optionally write details of TFBSs found in test set to detailed
TFBS hits file.

=head1 AUTHOR

  Andrew Kwon
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: tjkwon@cmmt.ubc.ca

=cut

use strict;

#use lib '/home/tjkwon/OPOSSUM/oPOSSUM3/scripts/standalone';
use lib '/space/devel/oPOSSUM3/scripts/standalone';

use oPossumGeneInclude;

use constant SSA => 'opossum_gene_ssa.pl';
use constant TCA => 'opossum_gene_tca.pl';
use constant ACSA => 'opossum_gene_acsa.pl';
use constant ACTCA => 'opossum_gene_actca.pl';

use constant SSA_PREFIX => 'gene_ssa';
use constant TCA_PREFIX => 'gene_tca';
use constant ACSA_PREFIX => 'gene_acsa';
use constant ACTCA_PREFIX => 'gene_actca';

use constant DFLT_MIN_SR => 3;
use constant DFLT_MAX_SR => 3;
use constant DFLT_MIN_CL => 3;
use constant DFLT_MAX_CL => 3;

use constant BASE_DIR => '/space/devel/oPOSSUM3/scripts/standalone'; # burgundy
#use constant BASE_DIR => '/home/tjkwon/OPOSSUM/oPOSSUM3/scripts/standalone'; # watson

my $species;
my $t_gene_file;
my $t_gene_id_type;
my $bg_gene_file;
my $bg_gene_id_type;
my $bg_num_rand_genes;
my $has_operon;
my $biotype;
my $tf_db;
my $anchor_tf_id;
my $max_site_dist;
my $collections_str;
my $tax_groups_str;
my $tf_file;
my $min_ic;
my $threshold;
my $min_cons_level;
my $max_cons_level;
my $min_sr_level;
my $max_sr_level;
my $num_results;
my $zscore_cutoff;
my $fisher_cutoff;
my $sort_by;
my $out_dir;
my $results_dir;
my $compute_cluster_dir;
my $prefix;
my $hits;
GetOptions(
    's=s'           => \$species,
    'g=s'           => \$t_gene_file,
    'b=s'           => \$bg_gene_file,
    'gt=s'          => \$t_gene_id_type,
    'bt=s'          => \$bg_gene_id_type,
    'bnr=i'         => \$bg_num_rand_genes,
    'has_operon'    => \$has_operon,
    'biotype=s'     => \$biotype,
    'db=s'          => \$tf_db,
    'id=s'          => \$anchor_tf_id,
    'dist=i'        => \$max_site_dist,
    'co=s'          => \$collections_str,
    'tf=s'          => \$tf_file,
    'tax=s'         => \$tax_groups_str,
    'ic=s'          => \$min_ic,
    'th=s'          => \$threshold,
    'min_cl=i'      => \$min_cons_level,
    'max_cl=i'      => \$max_cons_level,
    'min_sr=i'      => \$min_sr_level,
    'max_sr=i'      => \$max_sr_level,
    'n=s'           => \$num_results,   # integer or string 'All'
    'zcutoff=f'     => \$zscore_cutoff,
    'fcutoff=f'     => \$fisher_cutoff,
    'sr=s'          => \$sort_by,
    'd=s'           => \$results_dir,
    'od=s'          => \$out_dir,
    'cc=s'          => \$compute_cluster_dir,
    'p=s'           => \$prefix,
    'hits'          => \$hits
);

$min_sr_level = DFLT_MIN_SR if !$min_sr_level;
$max_sr_level = DFLT_MAX_SR if !$max_sr_level;
$min_cons_level = DFLT_MIN_CL if !$min_cons_level;
$max_cons_level = DFLT_MAX_CL if !$max_cons_level;

if (!$out_dir and $results_dir) {
    $out_dir = $results_dir;
} elsif (!$out_dir) {
    $out_dir = ".";
}

$sort_by = 'z-score' if !$sort_by;

if (!$species) {
    pod2usage(
        -msg     => "\nPlease specify the species name",
        -verbose => 1
    );
}

if (!$t_gene_file) {
    pod2usage(
        -msg     => "\nPlease specify an input gene ID file",
        -verbose => 1
    );
}


if ($bg_gene_file && $bg_num_rand_genes) {
    pod2usage(
        -msg => "\nPlease specify EITHER a background gene file"
            . " OR a number of random background genes OR"
            . " nothing (default = all genes in oPOSSUM DB)",
        -verbose => 1
    );
}

if (!$t_gene_id_type) {
    $t_gene_id_type = DFLT_GENE_ID_TYPE;
}

if (!$bg_gene_id_type) {
    $bg_gene_id_type = DFLT_GENE_ID_TYPE;
}

#
# command arguments
#
my $parameters = " -s $species -g $t_gene_file -gt $t_gene_id_type";

if ($bg_gene_file) {
    $parameters .= " -b $bg_gene_file -bt $bg_gene_id_type";
} elsif ($bg_num_rand_genes) {
    $parameters .= " -bnr $bg_num_rand_genes";
}

$parameters .= " -has_operon" if $has_operon;
$parameters .= " -biotype $biotype" if $biotype;
$parameters .= " -db $tf_db" if $tf_db;
$parameters .= " -co $collections_str" if $collections_str;
$parameters .= " -tf $tf_file" if $tf_file;
$parameters .= " -tax $tax_groups_str" if $tax_groups_str;
$parameters .= " -ic $min_ic" if $min_ic;
$parameters .= " -n $num_results" if $num_results;
$parameters .= " -zcutoff $zscore_cutoff" if $zscore_cutoff;
$parameters .= " -fcutoff $fisher_cutoff" if $fisher_cutoff;
$parameters .= " -sr $sort_by" if $sort_by;
$parameters .= " -d $results_dir" if $results_dir;

my $ssa_command;
my $tca_command;
my $acsa_command;
my $actca_command;

$ssa_command = SSA . $parameters;

if (!$tf_file) {
    $tca_command = TCA . $parameters;
    $actca_command = ACTCA . $parameters;
}

if ($anchor_tf_id) {
    # run anchored methods
    $acsa_command = ACSA . " -id $anchor_tf_id";
    $actca_command = ACTCA . " -id $anchor_tf_id";
}

my $opdba = opossum_db_connect($species)
    || die("Could not connect to oPOSSUM DB");

my $srla = $opdba->get_SearchRegionLevelAdaptor
    || die("Could not get SearchRegionLevelAdpator");

my $sr_level_hash = $srla->fetch_search_region_level_hash();

if ($compute_cluster_dir) {
    my $qsub_file = $compute_cluster_dir . "/qsub_" . $prefix . ".sh";
    open (QSUB, ">$qsub_file") or die "Can't create $qsub_file";
    chmod(0755, $qsub_file);
}
for (my $i = $min_sr_level; $i <= $max_sr_level; $i++)
{
    for (my $j = $min_cons_level; $j <= $max_cons_level; $j++)
    {
        if ($compute_cluster_dir) {
            my $script1 = generate_sge_script($ssa_command, $out_dir, $results_dir, $compute_cluster_dir, $hits, $prefix, SSA_PREFIX, $i, $j) if $ssa_command;
            my $script2 = generate_sge_script($tca_command, $out_dir, $results_dir, $compute_cluster_dir, $hits, $prefix, TCA_PREFIX, $i, $j) if $tca_command;
            #generate_sge_script($acsa_command, $out_dir, $results_dir, $compute_cluster_dir, $hits, $prefix, ACSA_PREFIX, $i, $j) if $acsa_command;
            #generate_sge_script($actca_command, $out_dir, $results_dir, $compute_cluster_dir, $hits, $prefix, ACTCA_PREFIX, $i, $j) if $actca_command;
            print QSUB "$script1\n$script2\n";
        } else {
            run_analysis($ssa_command, $out_dir, $results_dir, $hits, $prefix, SSA_PREFIX, $i, $j) if $ssa_command;
            run_analysis($tca_command, $out_dir, $results_dir, $hits, $prefix, TCA_PREFIX, $i, $j) if $tca_command;
            #run_anchored_analysis($acsa_command, $out_dir, $results_dir, $hits, $prefix, ACSA_PREFIX, $i, $j) if $acsa_command;
            #run_anchored_analysis($actca_command, $out_dir, $results_dir, $hits, $prefix, ACTCA_PREFIX, $i, $j) if $actca_command;
        }
    }
}
close(QSUB);

########################

sub generate_sge_script
{
    my ($command, $out_dir, $results_dir, $compute_cluster_dir, $hits, $prefix, $command_prefix,
        $sr_level, $cons_level) = @_;

    my $file_prefix = $prefix . "_$command_prefix" . "_sr$sr_level" . "_cl$cons_level";
    my $job_name = $file_prefix;
    my $log_file = $results_dir . "/$file_prefix.log";
    
    my $out_file = "$out_dir/$file_prefix.out";
    my $err_file = "$out_dir/$file_prefix.err";    

    my $results_file =  $results_dir . "/$file_prefix.results.txt";
    my $hits_file = $results_dir . "/$file_prefix.hits.txt" if $hits;
    my $full_command = BASE_DIR . "/$command -o $results_file -l $log_file";
    $full_command .= " -h $hits_file" if $hits;
    $full_command .= " -srl $sr_level -cl $cons_level";
    
    my $sh_file = $compute_cluster_dir . "/run_$file_prefix.sh";
    open (SH, ">$sh_file") or die "Can't create $sh_file\n";
            
    my $text  = qq{
#!/bin/bash

### Name of the job
#\$ -N $job_name
### Declare job is non-rerunable
#\$ -r n
### Export all environment variables to batch job
#\$ -V
#\$ -o $out_file
#\$ -e $err_file
### E-mail notification on job abort
#\$ -m a
#\$ -M tjkwon\@cmmt.ubc.ca

echo \$HOSTNAME

$full_command

exit 0
};
    
    print SH $text;
    
    close(SH);
    chmod(0755, $sh_file);
    
    return $sh_file;
}

sub run_analysis
{
    my ($command, $out_dir, $results_dir, $hits, $prefix, $command_prefix, $sr_level, $cons_level) = @_;
    
    my $file_prefix = $prefix . "_$command_prefix" . "_sr$sr_level" . "_cl$cons_level";
    my $results_file = "$results_dir/$file_prefix.results.txt";
    my $hits_file = "$results_dir/$file_prefix.hits.txt" if $hits;
    my $log_file = "$out_dir/$file_prefix.log";
    
    my $full_command = $command . " -o $results_file -l $log_file";
    $full_command .= " -h $hits_file" if $hits;
    $full_command .= " -srl $sr_level -cl $cons_level";
    
    #print "$full_command\n\n";
    system $full_command;
    
    return;
}

sub run_anchored_analysis
{
    my ($command, $out_dir, $results_dir, $hits, $prefix, $command_prefix, $sr_level, $cons_level) = @_;
    
    my $sr_up = $sr_level_hash->{$sr_level}->upstream_bp();
    my $sr_down = $sr_level_hash->{$sr_level}->downstream_bp();
    
    my $file_prefix = $prefix . "$command_prefix" . "_sr$sr_level" . "_cl$cons_level";
    my $results_file = "$results_dir/$file_prefix.results.txt";
    my $hits_file = "$results_dir/$file_prefix.hits.txt";
    my $log_file = "$out_dir/$file_prefix.log";
    
    my $full_command = $command . " -up $sr_up -down $sr_down"
        . " -id $anchor_tf_id -dist $max_site_dist"
        . " -o $results_file -l $log_file";
    $full_command .= " -h $hits_file" if $hits;
    $full_command .= " -srl $sr_level -cl $cons_level";
    
    system $full_command;
    
    return;
}
