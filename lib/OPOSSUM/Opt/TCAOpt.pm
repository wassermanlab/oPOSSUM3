#
# Specific options for oPOSSUM TCA analysis
#

#
# TFBS_Cluster system paths
#
use constant TFBS_CLUSTER_HOME       => '/apps/oPOSSUM3/lib/TFBSCluster';
use constant TFBS_CLUSTER_LIB_PATH   => TFBS_CLUSTER_HOME . '/lib';


#
# oPOSSUM_cluster DB Access
#
use constant TFBS_CLUSTER_DB_HOST    => 'opossum.cmmt.ubc.ca';
use constant TFBS_CLUSTER_DB_NAME    => 'TFBS_cluster';
use constant TFBS_CLUSTER_DB_USER    => 'opossum_r';
use constant TFBS_CLUSTER_DB_PASS    => '';

1;
