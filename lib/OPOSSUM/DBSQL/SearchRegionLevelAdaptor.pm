=head1 NAME

OPOSSUM::DBSQL::SearchRegionLevelAdaptor - Adaptor for MySQL queries to retrieve
and store search region levels.

=head1 SYNOPSIS

$srla = $db_adaptor->get_SearchRegionLevelAdaptor();

=head1 DESCRIPTION

In order to facilitate fast retrieval of TFBS counts from the oPOSSUM database
several count sets were pre-computed using discrete values for PWM (PSSM)
thresholds, conservation levels, and upstream/downstream search regions.
The search_region_levels table of the oPOSSUM database stores information about
the discrete search regions which were used. The records are stored with a level
number and the associated amount of upstream/downstream sequence around the
TSS used in the TFBS search.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package OPOSSUM::DBSQL::SearchRegionLevelAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::SearchRegionLevel;
use OPOSSUM::SearchRegionLevelSet;

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
 Usage    : $levels = $srla->fetch_levels();
 Function : Fetch a list of all the search region levels in the DB.
 Returns  : A reference to a list of integer search region levels.
 Args	  : None.

=cut

sub fetch_levels
{
    my ($self) = @_;

    my $sql = qq{select level from search_region_levels order by level};

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

=head2 fetch_search_region_levels

 Title    : fetch_search_region_levels
 Usage    : $srls = $srla->fetch_search_region_levels();
 Function : Alternate name for fetch_search_region_level_set.
 Returns  : A reference to an OPOSSUM::SearchRegionLevelSet object.
 Args	  : None.

=cut

sub fetch_search_region_levels
{
    my ($self) = @_;
    
    return $self->fetch_search_region_level_set();
}


=head2 fetch_search_region_level_set

 Title    : fetch_search_region_level_set
 Usage    : $srls = $srla->fetch_search_region_level_set();
 Function : Fetch a list of all the search region level objects in the DB.
 Returns  : A reference to an OPOSSUM::SearchRegionLevelSet object.
 Args	  : None.

=cut

sub fetch_search_region_level_set
{
    my ($self) = @_;

    my $levels = OPOSSUM::SearchRegionLevelSet->new();
    if (!$levels) {
    	carp "error creating new OPOSSUM::SearchRegionLevelSet\n";
        return;
    }

    my $sql = qq{select level, upstream_bp, downstream_bp
		    from search_region_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching search region levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching search region levels\n" . $self->errstr;
        return;
    }

    while (my @row = $sth->fetchrow_array) {
        $levels->add_search_region_level(
            OPOSSUM::SearchRegionLevel->new(
                -level          => $row[0],
                -upstream_bp    => $row[1],
                -downstream_bp  => $row[2]
            )
        );
    }
    $sth->finish;

    return $levels;
}

=head2 fetch_search_region_level_list

 Title    : fetch_search_region_level_list
 Usage    : $srll = $cla->fetch_search_region_level_list();
 Function : Fetch all the search region level objects from the DB.
 Returns  : A reference to a list of OPOSSUM::SearchRegionLevel objects.
 Args	  : None.

=cut

sub fetch_search_region_level_list
{
    my ($self) = @_;

    my $sql = qq{select level, upstream_bp, downstream_bp
		    from search_region_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching search region levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching search region levels\n" . $self->errstr;
        return;
    }

    my @sr_levels;
    while (my @row = $sth->fetchrow_array) {
        push @sr_levels, OPOSSUM::SearchRegionLevel->new(
            -level          => $row[0],
            -upstream_bp    => $row[1],
            -downstream_bp  => $row[2]
        );
    }
    $sth->finish;

    return @sr_levels ? \@sr_levels : undef;
}

=head2 fetch_search_region_level_hash

 Title    : fetch_search_region_level_hash
 Usage    : $srlh = $cla->fetch_search_region_level_hash();
 Function : Fetch all the search region level objects from the DB.
 Returns  : A reference to a hash of OPOSSUM::SearchRegionLevel objects.
 Args	  : None.

=cut

sub fetch_search_region_level_hash
{
    my ($self) = @_;

    my $sql = qq{select level, upstream_bp, downstream_bp
		    from search_region_levels order by level};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching search region levels\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching search region levels\n" . $self->errstr;
        return;
    }

    my %sr_levels;
    while (my @row = $sth->fetchrow_array) {
        $sr_levels{$row[0]} = OPOSSUM::SearchRegionLevel->new(
            -level          => $row[0],
            -upstream_bp    => $row[1],
            -downstream_bp  => $row[2]
        );
    }
    $sth->finish;

    return %sr_levels ? \%sr_levels : undef;
}

=head2 fetch_by_level

 Title    : fetch_by_level
 Usage    : $srl = $srla->fetch_by_level($level);
 Function : Fetch a search region level object by it's level number.
 Returns  : An OPOSSUM::SearchRegionLevel object.
 Args	  : Integer level.

=cut

sub fetch_by_level
{
    my ($self, $level) = @_;

    my $sql = qq{
        select upstream_bp, downstream_bp
		from search_region_levels
		where level = $level
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error fetching search region level $level\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error fetching search region level $level\n" . $self->errstr;
        return;
    }

    my $srl;
    if (my @row = $sth->fetchrow_array) {
        $srl = OPOSSUM::SearchRegionLevel->new(
            -level          => $level,
            -upstream_bp    => $row[0],
            -downstream_bp  => $row[1]
        );
    }
    $sth->finish;

    return $srl;
}

1;
