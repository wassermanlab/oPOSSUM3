#
# Specific options for gene-based analysis
#

# The maximum number of target genes a user is allowed to paste/upload
use constant MAX_TARGET_GENES           => 2000;

# Form selection values and default settings
use constant DFLT_BG_NUM_RAND_GENES     => 5000;

#
# oPOSSUM DB Access
#
# The species name is dynamically appended to OPOSSUM_DB_NAME for full
# oPOSSUM DB name, e.g. oPOSSUM3_human
#
use constant OPOSSUM_DB_HOST    => 'opossum.cmmt.ubc.ca';
use constant OPOSSUM_DB_NAME    => 'oPOSSUM3';
use constant OPOSSUM_DB_USER    => 'opossum_r';
use constant OPOSSUM_DB_PASS    => '';

#
# Form selection values and default settings. These may need to be overriden 
# by specific oPOSSUM variants.
#
use constant DFLT_GENE_ID_TYPE          => 0;
use constant DFLT_CONSERVATION_LEVEL    => 1;
use constant DFLT_SEARCH_REGION_LEVEL   => 3;
use constant DFLT_THRESHOLD_LEVEL       => 3;

# Leave as undef or set to 'all' for all biotypes (not tested).
use constant DFLT_BIOTYPE               => undef;

1;
