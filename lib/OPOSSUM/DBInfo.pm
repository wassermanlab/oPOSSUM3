=head1 NAME

OPOSSUM::DBInfo - DBInfo object (db_info DB record)

=head1 DESCRIPTION

A DBInfo object models the (single) record contained in the db_info
table of the oPOSSUM DB. The DBInfo object contains information about
how the oPOSSUM database was built, including the databases and software
versions which were used.

=head1 MODIFICATIONS

 2010/09/25 AK
 - added has_operon method
 
=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBInfo;

use strict;
use Carp;
use OPOSSUM::DBObject;

use vars qw(@ISA);

@ISA = qw(OPOSSUM::DBObject);

=head2 new

 Title   : new
 Usage   : $db_info = OPOSSUM::DBInfo->new(
                -build_date		        => '2010/01/01 17:23:06',
                -species	            => 'mouse',
                -latin_name		        => 'mus musculus'
                -assembly		        => 'NCBI37'
                -ensembl_db	            => 'mus_musculus_core_56_37i',
                -ucsc_db	            => 'mm9',
                -ucsc_cons_table        => 'phastCons30wayPlacental',
                -min_threshold          => 0.75,
                -max_upstream_bp	    => 10000,
                -min_cr_length	        => 20,
                -tax_group              => 'vertebrates',
                -min_ic                 => 8
				-has_operon				=> 0
           );

 Function: Construct a new DBInfo object
 Returns : a new OPOSSUM::DBInfo object

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {%args}, ref $class || $class;

    return $self;
}

=head2 build_date

 Title   : build_date
 Usage   : $date = $dbbi->build_date() or $dbbi->build_date($date);

 Function: Get/set the DB build date.
 Returns : A string.
 Args    : None or a new DB build date.

=cut

sub build_date
{
    my ($self, $build_date) = @_;

    if ($build_date) {
        $self->{-build_date} = $build_date;
    }

    return $self->{-build_date};
}

=head2 species

 Title   : species
 Usage   : $species = $dbbi->species() or $dbbi->species($species);

 Function: Get/set the name of the species
 Returns : A string
 Args    : None or a new species name

=cut

sub species
{
    my ($self, $species) = @_;

    if ($species) {
        $self->{-species} = $species;
    }

    return $self->{-species};
}

=head2 latin_name

 Title   : latin_name
 Usage   : $latin_name = $dbbi->latin_name() or $dbbi->latin_name($latin_name);

 Function: Get/set the species latin name
 Returns : A string
 Args    : None or a new species latin name

=cut

sub latin_name
{
    my ($self, $latin_name) = @_;

    if ($latin_name) {
        $self->{-latin_name} = $latin_name;
    }

    return $self->{-latin_name};
}

=head2 assembly

 Title   : assembly
 Usage   : $assembly = $dbbi->assembly() or $dbbi->assembly($assembly);

 Function: Get/set the name of the genome assembly.
 Returns : A string.
 Args    : None or a new genome assembly name.

=cut

sub assembly
{
    my ($self, $assembly) = @_;

    if ($assembly) {
        $self->{-assembly} = $assembly;
    }

    return $self->{-assembly};
}

=head2 ensembl_db

 Title   : ensembl_db
 Usage   : $db_name = $dbbi->ensembl_db()
           or $dbbi->ensembl_db($db_name);

 Function: Get/set the name of the Ensembl species core database.
 Returns : A string.
 Args    : None or a new Ensembl species core database name.

=cut

sub ensembl_db
{
    my ($self, $db_name) = @_;

    if ($db_name) {
        $self->{-ensembl_db} = $db_name;
    }

    return $self->{-ensembl_db};
}

=head2 ucsc_db

 Title   : ucsc_db
 Usage   : $db_name = $dbbi->ucsc_db()
           or $dbbi->ucsc_db($db_name);

 Function: Get/set the name of the UCSC database.
 Returns : A string.
 Args    : None or a new UCSC database name.

=cut

sub ucsc_db
{
    my ($self, $db_name) = @_;

    if ($db_name) {
        $self->{-ucsc_db} = $db_name;
    }

    return $self->{-ucsc_db};
}

=head2 ucsc_cons_table

 Title   : ucsc_cons_table
 Usage   : $table_name = $dbbi->ucsc_cons_table()
           or $dbbi->ucsc_cons_table($table_name);

 Function: Get/set the UCSC conservation (phastCons) table name.
 Returns : A string.
 Args    : None or a new UCSC conservation (phastCons) table name.

=cut

sub ucsc_cons_table
{
    my ($self, $table_name) = @_;

    if ($table_name) {
        $self->{-ucsc_cons_table} = $table_name;
    }

    return $self->{-ucsc_cons_table};
}

=head2 min_threshold

 Title   : min_threshold
 Usage   : $min_score = $dbbi->min_threshold()
           or $dbbi->min_threshold($min_score);

 Function: Get/set the minimum TFBS matrix score threshold used when building
           the database.
 Returns : A float.
 Args    : None or a new minimum score threshold.

=cut

sub min_threshold
{
    my ($self, $min_threshold) = @_;

    if (defined $min_threshold) {
        $self->{-min_threshold} = $min_threshold;
    }

    return $self->{-min_threshold};
}

sub min_tfbs_score
{
    my ($self, $min_score) = @_;

    carp "deprecated method min_tfb_score(); please use min_threshold()\n";

    return $self->min_threshold($min_score);
}

=head2 max_upstream_bp

 Title   : max_upstream_bp
 Usage   : $upstream_bp = $dbbi->max_upstream_bp()
 	       or $dbbi->max_upstream_bp($upstream_bp);

 Function: Get/set the maximum amount of upstream bp used in the DB.
 Returns : A string.
 Args    : None or a new max upstream bp amount.

=cut

sub max_upstream_bp
{
    my ($self, $upstream_bp) = @_;

    if ($upstream_bp) {
        $self->{-max_upstream_bp} = $upstream_bp;
    }

    return $self->{-max_upstream_bp};
}

# deprecated
sub upstream_bp
{
    my ($self, $upstream_bp) = @_;

    carp "deprecated method upstream_bp(); please use max_upstream_bp()\n";

    return $self->max_upstream_bp($upstream_bp);
}


#=head2 max_downstream_bp
#
# Title   : max_downstream_bp
# Usage   : $downstream_bp = $dbbi->max_downstream_bp()
# 	       or $dbbi->max_downstream_bp($downstream_bp);
#
# Function: Get/set the maximum amount of downstream bp used in the DB.
# Returns : A string.
# Args    : None or a new max downstream bp amount.
#
#=cut

#sub max_downstream_bp
#{
#    my ($self, $downstream_bp) = @_;
#
#    if ($downstream_bp) {
#        $self->{-max_downstream_bp} = $downstream_bp;
#    }
#
#    return $self->{-max_downstream_bp};
#}

=head2 min_cr_length

 Title   : min_cr_length
 Usage   : $min_cr_length = $dbbi->min_cr_length()
           or $dbbi->min_cr_length($min_cr_length);

 Function: Get/set the minimum lenght of a conserved region to include
           in the conserved_regions table.
 Returns : A string.
 Args    : None or a new upstream bp amount.

=cut

sub min_cr_length
{
    my ($self, $min_cr_length) = @_;

    if ($min_cr_length) {
        $self->{-min_cr_length} = $min_cr_length;
    }

    return $self->{-min_cr_length};
}

=head2 tax_group

 Title   : tax_group
 Usage   : $tax_group = $dbbi->tax_group()
           or $dbbi->tax_group($tax_group);

 Function: Get/set the taxonomic supergroup(s) that this species belongs to.
           This could be a single tax group or a comma separated string of
           tax groups from the JASPAR CORE collection ('vertebrates',
           'plants', 'nematodes', 'fungi' etc). Note the pluralization of
           the names. 
           NOTE: For worms, this is not just one, but all metazoans.
 Returns : A tax group string.
 Args    : None or a new tax group string.

=cut

sub tax_group
{
    my ($self, $tax_group) = @_;

    if (defined $tax_group) {
        $self->{-tax_group} = $tax_group;
    }

    return $self->{-tax_group};
}

=head2 min_ic

 Title   : min_ic
 Usage   : $min_ic = $dbbi->min_ic()
           or $dbbi->min_ic($min_ic);

 Function: Get/set the TFBS profile matrix minimum information content
           used when building the database.
 Returns : An integer.
 Args    : None or a new matrix minimum information content.

=cut

sub min_ic
{
    my ($self, $min_ic) = @_;

    if (defined $min_ic) {
        $self->{-min_ic} = $min_ic;
    }

    return $self->{-min_ic};
}

sub min_pwm_ic
{
    my ($self, $min_ic) = @_;

    carp "deprecated method min_pwm_ic(); please use min_ic()\n";

    return $self->min_ic($min_ic);
}


=head2 has_operon 

 Title   : has_operon
 Usage   : $has_operon = $dbbi->has_operon()
           or $dbbi->has_operon($has_operon);

 Function: Get/set the has_operon flag, which is set to 1 if the species
           has operons in its genome.
 Returns : Integer.
 Args    : None or a new integer (0 or 1).

=cut

sub has_operon
{
    my ($self, $has_operon) = @_;

    if (defined $has_operon) {
        $self->{-has_operon} = $has_operon;
    }

    return $self->{-has_operon};
}

1;
