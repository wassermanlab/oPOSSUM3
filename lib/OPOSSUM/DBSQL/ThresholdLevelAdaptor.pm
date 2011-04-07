
=head1 NAME

OPOSSUM::DBSQL::ThresholdLevelAdaptor - Adaptor for MySQL queries to retrieve
and store PWM (PSSM) threshold levels.

=head1 SYNOPSIS

$tla = $db_adaptor->get_ThresholdLevelAdaptor();

=head1 DESCRIPTION

In order to facilitate fast retrieval of TFBS counts from the oPOSSUM database
several count sets were pre-computed using discrete values for PWM (PSSM)
thresholds, conservation levels, and upstream/downstream search regions.
The threshold_levels table of the oPOSSUM database stores information about
the discrete PWM (PSSM) threshold scores which were used, i.e. the minimum
matrix score used to determine a TFBS 'hit'. The records are stored with a level
number and the associated matrix score threshold used in the TFBS search.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ThresholdLevelAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::ThresholdLevel;
use OPOSSUM::ThresholdLevelSet;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_levels

 Title    : fetch_levels
 Usage    : $levels = $tla->fetch_levels();
 Function : Fetch a list of all the PWM threshold levels in the DB.
 Returns  : A reference to a list of integer PWM threshold levels.
 Args	  : None.

=cut

sub fetch_levels
{
    my ($self) = @_;

    my $sql = qq{select level from threshold_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching levels\n" . $self->errstr;
        return;
    }

    my @levels;
    while (my ($level) = $sth->fetchrow_array) {
        push @levels, $level;
    }

    return @levels ? \@levels : undef;
}

=head2 fetch_threshold_levels

 Title    : fetch_threshold_levels
 Usage    : $cls = $cla->fetch_threshold_levels();
 Function : Alternate name for fetch_threshold_level_set (should really
            call fetch_threshold_level_list but this is for backward
            compatibility.

=cut

sub fetch_threshold_levels
{
    my ($self) = @_;

    return $self->fetch_threshold_level_set();
}

=head2 fetch_threshold_level_set

 Title    : fetch_threshold_level_set
 Usage    : $tls = $tla->fetch_threshold_level_set();
 Function : Fetch a set of all the threshold level objects in the DB.
 Returns  : A reference to an OPOSSUM::ThresholdLevelSet object.
 Args	  : None.

=cut

sub fetch_threshold_level_set
{
    my ($self) = @_;

    my $levels = OPOSSUM::ThresholdLevelSet->new();
    if (!$levels) {
        carp "error creating new ThresholdLevelSet object\n";
        return;
    }

    my $sql = qq{select level, threshold from threshold_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching threshold levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching threshold levels\n" . $self->errstr;
        return;
    }

    while (my @row = $sth->fetchrow_array) {
        $levels->add_threshold_level(
            OPOSSUM::ThresholdLevel->new(
                -level     => $row[0],
                -threshold => $row[1]
            )
        );
    }
    $sth->finish;

    return $levels;
}

=head2 fetch_threshold_level_list

 Title    : fetch_threshold_level_list
 Usage    : $thll = $cla->fetch_threshold_level_list();
 Function : Fetch all the threshold level objects from the DB.
 Returns  : A reference to a list of OPOSSUM::ThresholdLevel objects.
 Args	  : None.

=cut

sub fetch_threshold_level_list
{
    my ($self) = @_;

    my $sql = qq{select level, threshold from threshold_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching threshold levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching threshold levels\n" . $self->errstr;
        return;
    }

    my @thresh_levels;
    while (my @row = $sth->fetchrow_array) {
        push @thresh_levels, OPOSSUM::ThresholdLevel->new(
            -level          => $row[0],
            -threshold      => $row[1]
        );
    }
    $sth->finish;

    return @thresh_levels ? \@thresh_levels : undef;
}

=head2 fetch_threshold_level_hash

 Title    : fetch_threshold_level_hash
 Usage    : $thlh = $cla->fetch_threshold_level_hash();
 Function : Fetch all the threshervation level objects from the DB.
 Returns  : A reference to a hash of OPOSSUM::ThresholdLevel objects.
 Args	  : None.

=cut

sub fetch_threshold_level_hash
{
    my ($self) = @_;

    my $sql = qq{select level, threshold from threshold_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching threshold levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching threshold levels\n" . $self->errstr;
        return;
    }

    my %thresh_levels;
    while (my @row = $sth->fetchrow_array) {
        $thresh_levels{$row[0]} = OPOSSUM::ThresholdLevel->new(
            -level      => $row[0],
            -threshold  => $row[1]
        );
    }
    $sth->finish;

    return %thresh_levels ? \%thresh_levels : undef;
}

=head2 fetch_by_level

 Title    : fetch_by_level
 Usage    : $tl = $tla->fetch_by_level($level);
 Function : Fetch a threshold level object by it's level number.
 Returns  : An OPOSSUM::ThresholdLevel object.
 Args	  : Integer level.

=cut

sub fetch_by_level
{
    my ($self, $level) = @_;

    my $sql = qq{select threshold from threshold_levels where level = $level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching threshold by level $level\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching threshold by level $level\n" . $self->errstr;
        return;
    }

    my $threshold_level;
    if (my @row = $sth->fetchrow_array) {
        $threshold_level = OPOSSUM::ThresholdLevel->new(
            -level      => $level,
            -threshold  => $row[0]
        );
    }
    $sth->finish;

    return $threshold_level;
}

1;
