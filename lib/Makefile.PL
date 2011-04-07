require 5.008;
use ExtUtils::MakeMaker;

my $NAME     = 'POSSUM';
my $DISTNAME = 'POSSUM';
my $VERSION  = '0.1.0';

WriteMakefile(
    NAME        => $NAME,
    DISTNAME    => $DISTNAME,
    VERSION     => $VERSION,
    'dist'      => {
		    COMPRESS     => 'gzip -9f',
		    SUFFIX       => '.gz',
		    DIST_DEFAULT => 'all tardist'
		}
);
