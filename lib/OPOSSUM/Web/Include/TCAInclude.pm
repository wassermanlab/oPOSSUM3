#package oPossumTCAWebInclude2;


# This module should be included in all the oPossum*TCAWeb.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to all the oPossum TCA variants.
#

use OPOSSUM::Web::Opt::BaseOpt;
use OPOSSUM::Web::Opt::SeqOpt;
use OPOSSUM::Opt::TCAOpt;

use lib TFBS_CLUSTER_LIB_PATH;

use TFBSCluster::DBSQL::DBAdaptor;

use strict;


sub cldba
{
    my $self = shift;

    if (@_) {
        $self->{-cldba} = shift;
    }

    return $self->{-cldba};
}


sub opossum_cluster_db_connect
{
    my $self = shift;

    my $dbh = TFBSCluster::DBSQL::DBAdaptor->new(
        -host     => TFBS_CLUSTER_DB_HOST,
        -dbname   => TFBS_CLUSTER_DB_NAME,
        -user     => TFBS_CLUSTER_DB_USER,
        -password => TFBS_CLUSTER_DB_PASS
    );

    if (!$dbh) {
        $self->_error("Could not connect to oPOSSUM_cluster database "
                      . TFBS_CLUSTER_DB_NAME);
    }

    $self->cldba($dbh);
}

sub fetch_tf_cluster_set
{
    my ($self, %args) = @_;

    my %matrix_args = %args;

    unless ($matrix_args{-matrixtype}) {
        $matrix_args{-matrixtype} = 'PFM';
    }

    #printf STDERR "fetch_tf_cluster_set: matrix_args = \n"
    #    . Data::Dumper::Dumper(%matrix_args) . "\n";

    my $cldba = $self->cldba();
    unless ($cldba) {
        $cldba = $self->opossum_cluster_db_connect();
    }
    
    my $tfca = $cldba->get_TFClusterAdaptor;

    my $cluster_set = TFBSCluster::TFClusterSet->new();

    my $clusters;
    #print STDERR "Getting cluster set\n";
    if ($matrix_args{-families}) {
        $clusters = $tfca->fetch_by_tf_families($matrix_args{-families});
    } elsif ($matrix_args{-family}) {
		$clusters = $tfca->fetch_by_tf_families($matrix_args{-family});
	} else {
        $clusters = $tfca->fetch_all();    
    }
    #print STDERR "# clusters = " . scalar(@$clusters) . "\n";
    $cluster_set->add_tf_cluster_list($clusters);

    die "Could not fetch TFBS clusters\n"
        if !$cluster_set || $cluster_set->size == 0;

    return $cluster_set;
}

sub get_tf_family_file
{    
    my ($self, $fam_input_method, $tempdir) = @_;

    my $q = $self->query;
    
    my $filename;
    if ($fam_input_method eq "specific") {
        my $sl = $q->param('tf_families');
        if (!$sl) {
            $self->_error("TF families not specified");
            return;
        }

        $filename = "$tempdir/tf_family.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create TF family list file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } elsif ($fam_input_method eq "upload") {
        my $file = $q->param('tf_family_upload_file');
        my $fh   = $q->upload('tf_family_upload_file');

        my $sl;
        while (my $line = <$fh>) {
            $sl .= $line;
        }

        if (!$sl) {
            $self->_error("File $file is empty\n");
            return;
        }

        $filename = "$tempdir/tf_family.txt";
        unless (open(FH, ">$filename")) {
            $self->_error("Unable to create TF family list file $filename\n");
            return;
        }
        print FH $sl;
        close(FH);
    } else {
        $self->_error("Unknown TF family list input method");
        return;
    }

    return $filename;
}

1;
