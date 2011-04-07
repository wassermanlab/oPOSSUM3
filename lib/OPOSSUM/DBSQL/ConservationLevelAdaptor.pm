
=head1 NAME

OPOSSUM::DBSQL::ConservationLevelAdaptor - Adaptor for MySQL queries to retrieve
and store conservation levels.

=head1 SYNOPSIS

$cla = $db_adaptor->get_ConservationLevelAdaptor();

=head1 DESCRIPTION

In order to facilitate fast retrieval of TFBS counts from the oPOSSUM database
several count sets were pre-computed using discrete values for PWM (PSSM)
thresholds, conservation levels, and upstream/downstream search regions.
The conservation_levels table of the oPOSSUM database stores information about
the conserved regions which were used to restrict the search for TFBSs. The
records are stored with a level number and the associated minimum conservation
cutoff at that level.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::ConservationLevelAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::ConservationLevel;
use OPOSSUM::ConservationLevelSet;

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
 Usage    : $levels = $cla->fetch_levels();
 Function : Fetch a list of all the conservation levels in the DB.
 Returns  : A reference to a list of integer conservation levels.
 Args	  : None.

=cut

sub fetch_levels
{
    my ($self) = @_;

    my $sql = qq{select level from conservation_levels order by level};

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

=head2 fetch_conservation_levels

 Title    : fetch_conservation_levels
 Usage    : $cls = $cla->fetch_conservation_levels();
 Function : Alternate name for fetch_conservation_level_set (should really
            call fetch_conservation_level_list but this is for backward
            compatibility.

=cut

sub fetch_conservation_levels
{
    my ($self) = @_;

    return $self->fetch_conservation_level_set();
}

=head2 fetch_conservation_level_set

 Title    : fetch_conservation_level_set
 Usage    : $cls = $cla->fetch_conservation_level_set();
 Function : Fetch all the conservation level objects from the DB.
 Returns  : A an OPOSSUM::ConservationLevelSet object.
 Args	  : None.

=cut

sub fetch_conservation_level_set
{
    my ($self) = @_;

    my $levels = OPOSSUM::ConservationLevelSet->new;
    if (!$levels) {
        carp "error creating OPOSSUM::ConservationLevelSet";
        return;
    }

    my $sql = qq{select level, min_conservation
        from conservation_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching conservation levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching conservation levels\n" . $self->errstr;
        return;
    }

    while (my @row = $sth->fetchrow_array) {
        $levels->add_conservation_level(
            OPOSSUM::ConservationLevel->new(
                -level              => $row[0],
                -min_conservation   => $row[1]
            )
        );
    }
    $sth->finish;

    return $levels;
}

=head2 fetch_conservation_level_list

 Title    : fetch_conservation_level_list
 Usage    : $cll = $cla->fetch_conservation_level_list();
 Function : Fetch all the conservation level objects from the DB.
 Returns  : A reference to a list of OPOSSUM::ConservationLevel objects.
 Args	  : None.

=cut

sub fetch_conservation_level_list
{
    my ($self) = @_;

    my $sql = qq{select level, min_conservation
        from conservation_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching conservation levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching conservation levels\n" . $self->errstr;
        return;
    }

    my @cons_levels;
    while (my @row = $sth->fetchrow_array) {
        push @cons_levels, OPOSSUM::ConservationLevel->new(
            -level          => $row[0],
            -min_conservation => $row[1]
        );
    }
    $sth->finish;

    return @cons_levels ? \@cons_levels : undef;
}

=head2 fetch_conservation_level_hash

 Title    : fetch_conservation_level_hash
 Usage    : $clh = $cla->fetch_conservation_level_hash();
 Function : Fetch all the conservation level objects from the DB.
 Returns  : A reference to a hash of OPOSSUM::ConservationLevel objects.
 Args	  : None.

=cut

sub fetch_conservation_level_hash
{
    my ($self) = @_;

    my $sql = qq{select level, min_conservation
		    from conservation_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching conservation levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching conservation levels\n" . $self->errstr;
        return;
    }

    my %cons_levels;
    while (my @row = $sth->fetchrow_array) {
        $cons_levels{$row[0]} = OPOSSUM::ConservationLevel->new(
            -level          => $row[0],
            -min_conservation => $row[1]
        );
    }
    $sth->finish;

    return %cons_levels ? \%cons_levels : undef;
}

=head2 fetch_by_level

 Title    : fetch_by_level
 Usage    : $cons_level = $cla->fetch_by_level($level);
 Function : Fetch a conservation level object from the DB by it's
            level number.
 Returns  : An OPOSSUM::ConservationLevel object.
 Args	  : Integer level.

=cut

sub fetch_by_level
{
    my ($self, $level) = @_;

    my $sql = qq{select min_conservation from conservation_levels
		where level = $level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching conservation level $level\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching conservation level $level\n" . $self->errstr;
        return;
    }

    my $cl;
    if (my @row = $sth->fetchrow_array) {
        $cl = OPOSSUM::ConservationLevel->new(
            -level          => $level,
            -min_conservation => $row[0]
        );
    }
    $sth->finish;

    return $cl;
}

=head2 store

 Title    : store
 Usage    : $ok = $cla->store($cl);
 Function : Store a conservation level object in the DB.
 Returns  : True on success, false otherwise.
 Args	  : An OPOSSUM::ConservationLevel object.

=cut

sub store
{
    my ($self, $cl) = @_;

    if (!ref $cl || !$cl->isa('OPOSSUM::ConservationLevel')) {
        carp "Not an OPOSSUM::ConservationLevel object";
        return;
    }

    my $sql = qq{insert into conservation_levels
        (level, min_conservation) values (?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert conservation level statement\n"
            . $self->errstr;
        return;
    }

    if (!$sth->execute($cl->level, $cl->min_conservation)) {
        carp "Error inserting conservation level\n" . $self->errstr;
        return;
    }
    $sth->finish;

    return 1;
}

1;
