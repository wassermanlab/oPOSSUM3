=head1 NAME

TFBSCluster::DBSQL::DBInfoAdaptor - Adaptor for MySQL queries to retrieve and
store global DB information.

=head1 SYNOPSIS

$cla = $db_adaptor->get_DBInfoAdaptor();

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 METHODS

=cut

package TFBSCluster::DBSQL::DBInfoAdaptor;

use strict;

use Carp;

use TFBSCluster::DBSQL::BaseAdaptor;
use TFBSCluster::DBInfo;

use vars qw(@ISA);
@ISA = qw(TFBSCluster::DBSQL::BaseAdaptor);

=head2 new

 Title    : new
 Usage    : $dbbia = TFBSCluster::DBSQL::DBInfoAdaptor->new(@args);
 Function : Create a new DBInfoAdaptor.
 Returns  : A new TFBSCluster::DBSQL::DBInfoAdaptor object.
 Args	  : An TFBSCluster::DBSQL::DBConnection object.

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_db_info

 Title    : fetch_db_info
 Usage    : $db_info = $cla->fetch_db_info();
 Function : Fetch the global database information.
 Returns  : An TFBSCluster::DBInfo object.
 Args	  : None.

=cut

sub fetch_db_info
{
    my $self = shift;

    my $sql = qq{
        select date_format(build_date, '%Y/%m/%d %T'),
			jaspar_db,
			collections,
			tax_groups,
			min_ic,
			cluster_threshold,
			radius_margin,
			comparison_method,
			cluster_method
        from db_info
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "error reading TFBSCluster DB info\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "error reading TFBSCluster DB info\n" . $self->errstr;
        return;
    }
    my @row = $sth->fetchrow_array;
    $sth->finish;
    if (!@row) {
        carp "error reading TFBSCluster DB info\n" . $self->errstr;
        return;
    }

    my $db_info = TFBSCluster::DBInfo->new(
        -build_date         => $row[0],
		-jaspar_db			=> $row[1],
		-collections		=> $row[2],
		-tax_groups			=> $row[3],
        -min_ic             => $row[4],
		-cluster_threshold	=> $row[5],
		-radius_margin		=> $row[6],
		-comparison_method	=> $row[7],
		-cluster_method		=> $row[8]
    );

    return $db_info;
}

=head2 store

 Title    : store
 Usage    : $cla->store($db_info);
 Function : Store the global database information.
 Returns  : TRUE on success, FALSE otherwise
 Args     : An TFBSCluster::DBInfo object.

=cut

sub store
{
    my ($self, $db_info) = @_;

    if (!ref $db_info || !$db_info->isa('TFBSCluster::DBInfo')) {
        carp "Not an TFBSCluster::DBInfo object";
        return 0;
    }

    my $sql = qq{
        insert into db_info (build_date, jaspar_db, collections, tax_groups,
			min_ic, cluster_threshold, radius_margin, comparison_method,
			cluster_method)
        values (now(),?,?,?,?,?,?,?,?)
    };

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert db_info statement\n" . $self->errstr;
        return 0;
    }

    #
    # Using built-in MySQL now() function in SQL statement above.
    #
    #my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    #
    #$year += 1900;  # localtime returns year as years since 1900
    #$mon += 1;      # localtime returns month in range 0-11

    # Standard MySQL datetime string format
    #my $datestr = "$year-$mon-$mday $hour:$min:$sec";

    my $ok = $sth->execute(
        #$datestr,
		$db_info->jaspar_db(),
		$db_info->collections(),
		$db_info->tax_groups(),
        $db_info->min_ic(),
		$db_info->cluster_threshold(),
		$db_info->radius_margin(),
		$db_info->comparison_method(),
		$db_info->cluster_method()
    );
    
    unless ($ok) {
        carp "Error inserting db_info\n" . $self->errstr();
    }

    $sth->finish;

    return $ok;
}

1;
