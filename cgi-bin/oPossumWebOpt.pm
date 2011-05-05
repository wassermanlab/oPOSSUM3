#
# General options file for all oPOSSUM web application variants
# (SSA, CSA, Seq, etc.).
#
# Specific options and overriden options should be stored in files
# specific to the oPOSSUM variant, e.g. oPossumSSAWebOpt.pm,
# oPossumCSAWebOpt.pm, oPossumSeqWebOpt.pm, etc.
#
# The oPOSSUM application variant web module (e.g. oPossumCSAWeb.pm) should
# include this options file and then any options file specific to that
# variant, e.g.:
#
# use oPossumWebOpt.pm;
# use oPossumCSAWebOpt.pm;
#

#
# General constants
#
use constant VERSION                    => '3.0';
use constant DEVEL_VERSION              => 1;
use constant ADMIN_EMAIL                => 'opossum@cmmt.ubc.ca';
use constant RESULTS_TEXT_FILENAME      => 'results.txt';
use constant RESULTS_HTDOCS_FILENAME    => 'results.html';

#
# oPOSSUM system paths
#
use constant OPOSSUM_HOME           => '/usr/local/apps/oPOSSUM3';
use constant OPOSSUM_LIB_PATH       => OPOSSUM_HOME . '/lib';
use constant OPOSSUM_HTDOCS_PATH    => OPOSSUM_HOME . '/htdocs';
use constant OPOSSUM_CGI_BIN_PATH   => OPOSSUM_HOME . '/cgi-bin';
use constant OPOSSUM_TMP_PATH       => OPOSSUM_HOME . '/tmp';
use constant OPOSSUM_LOG_PATH       => OPOSSUM_HOME . '/logs';

#
# oPOSSUM web server paths and URLs
#
use constant WEB_SERVER_URL             => 'http://opossum.cisreg.ca';
use constant WEB_SERVER_HOME            => '/var/www';
use constant ABS_HTDOCS_PATH            => WEB_SERVER_HOME
                                           . '/htdocs/oPOSSUM3';
use constant ABS_CGI_BIN_PATH           => WEB_SERVER_HOME
                                           . '/cgi-bin/oPOSSUM3';
use constant ABS_HTDOCS_TEMPLATE_PATH   => ABS_HTDOCS_PATH . '/templates';
use constant ABS_HTDOCS_TMP_PATH        => ABS_HTDOCS_PATH . '/tmp';
use constant ABS_HTDOCS_RESULTS_PATH    => ABS_HTDOCS_PATH . '/results';
use constant ABS_HTDOCS_DATA_PATH       => ABS_HTDOCS_PATH . '/data';
use constant REL_HTDOCS_PATH            => '/oPOSSUM3';
use constant REL_CGI_BIN_PATH           => '/cgi-bin/oPOSSUM3';
use constant REL_HTDOCS_TMP_PATH        => REL_HTDOCS_PATH . '/tmp';
use constant REL_HTDOCS_RESULTS_PATH    => REL_HTDOCS_PATH . '/results';
use constant REL_HTDOCS_DATA_PATH       => REL_HTDOCS_PATH . '/data';

#
# External URLs
#
use constant JASPAR_URL => 'http://jaspar.genereg.net/cgi-bin/jaspar_db.pl';
use constant PAZAR_URL  => 'http://www.pazar.info/cgi-bin/display_JASPAR_profile.cgi';

#
# JASPAR DB Access
#
use constant JASPAR_DB_HOST     => 'vm5.cmmt.ubc.ca';
use constant JASPAR_DB_NAME     => 'JASPAR_2010';
use constant JASPAR_DB_USER     => 'jaspar_r';
use constant JASPAR_DB_PASS     => '';

use constant LOW_MATRIX_GC  => 0.33;
use constant HIGH_MATRIX_GC => 0.66;
use constant LOW_MATRIX_IC  => 9;
use constant HIGH_MATRIX_IC => 19;
use constant LOW_SEQ_GC     => 0.45;
use constant HIGH_SEQ_GC    => 0.55;

#
# Form selection values and default settings. These may need to be overriden 
# by specific oPOSSUM variants.
#
use constant NUM_RESULTS                => [5, 10, 20, 'All'];
use constant ZSCORE_CUTOFFS             => [5, 10, 15];
use constant DFLT_ZSCORE_CUTOFF         => 10;
use constant FISHER_CUTOFFS             => [5, 7, 9];
use constant DFLT_FISHER_CUTOFF         => 7;
use constant DFLT_NUM_RESULTS           => 'All';
use constant DFLT_THRESHOLD             => 80;

# temp. file cleanup no. of days
use constant REMOVE_TEMPFILES_OLDER_THAN    => 3;
use constant REMOVE_RESULTFILES_OLDER_THAN  => 14;

1;
