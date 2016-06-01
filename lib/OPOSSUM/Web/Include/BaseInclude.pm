#
# This module should be included in all the oPossum*Web.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to all the oPossum variants. This
# includes utility functions as well as the common template routines like
# 'errors' and 'warnings'.
#

use OPOSSUM::Web::Opt::BaseOpt;

use lib OPOSSUM_LIB_PATH;

use Data::Dumper;    # for debugging only
use File::Temp qw/ tempfile tempdir /;
use File::Path qw/ rmtree /;

use OPOSSUM::TFSet;
use TFBS::DB::JASPAR5;
use TFBS::MatrixSet;

use Template;
use CGI::Carp qw(carpout);    # fatalsToBrowser;

use strict;

#
# High-level error routine. Call low level _error routine with current error
# and output all current errors to HTML error template.
#
sub error
{
    my ($self, $error) = @_;

    $self->_error($error) if $error;

    my $errors = $self->errors();

    my $err_str = join "\n", @$errors;
    carp "\nERROR:\n$err_str\n";

    my $error_html;
    foreach my $err (@$errors) {
        chomp $err;
        $err =~ s/\n/<br>/g;
        $error_html .= "$err<br>";
    }

    my $state = $self->state;

    my $vars = {
        abs_htdocs_path  => ABS_HTDOCS_PATH,
        rel_htdocs_path  => REL_HTDOCS_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        bg_color_class   => $state->bg_color_class(),
        title            => $state->title(),
        heading          => $state->heading(),
        section          => 'Error',
        version          => VERSION,
        devel_version    => DEVEL_VERSION,
        error            => $error_html,
        var_template     => "error.html"
    };

    my $output = $self->process_template('master.html', $vars);

    $self->errors(undef);

    return $output;
}

#
# High-level warning routine. Call low level _warning routine with current
# warning and output all current warnings to HTML warning template.
#
sub warning
{
    my ($self, $warning) = @_;

    $self->_warning($warning) if $warning;

    my $warnings = $self->warnings();

    my $warning_html;
    foreach my $warn (@$warnings) {
        chomp $warn;
        $warn =~ s/\n/<br>/g;
        $warning_html .= "$warn<br>";
    }

    my $state = $self->state;

    my $vars = {
        abs_htdocs_path  => ABS_HTDOCS_PATH,
        rel_htdocs_path  => REL_HTDOCS_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        bg_color_class   => $state->bg_color_class(),
        title            => $state->title(),
        heading          => $state->heading(),
        section          => 'Warning',
        version          => VERSION,
        devel_version    => DEVEL_VERSION,
        warning          => $warning_html,
        var_template     => "warning.html"
    };

    my $output = $self->process_template('master.html', $vars);

    $self->warnings(undef);

    return $output;
}

sub process_template
{
    my ($self, $template_name, $vars) = @_;

    my $config = {
        ABSOLUTE     => 1,
        INCLUDE_PATH => ABS_HTDOCS_TEMPLATE_PATH . "/", # or list ref
        INTERPOLATE  => 1,                       # expand "$var" in plain text
        POST_CHOMP   => 1,                       # cleanup whitespace
        #PRE_PROCESS  => 'header',                # prefix each template
        EVAL_PERL    => 1    # evaluate Perl code blocks
    };

    my $string   = '';
    my $template = Template->new($config);
    my $input    = ABS_HTDOCS_TEMPLATE_PATH . "/$template_name";
    $template->process($input, $vars, \$string) || die $template->error();

    return $string;
}

#
# Andrew - March 4, 2012
# Was wondering if I could use this to display the results.html created by
# an analysis script by oPossumWeb/oPossumGeneSSAWeb, but no go
# says the file is not found--some kind of permission problem?
#
=head3
sub process_html
{
	my ($self, $html_file) = @_;

	my $config = {
		ABSOLUTE	=> 1,
		POST_CHOMP	=> 1,
		EVAL_PERL	=> 1
	};

	my $string = '';
	my $vars;
	my $template = Template->new($config);
	$template->process($html_file, $vars, \$string) || die $template->error();

	return $string;
}
=cut

sub jdbh
{
    my $self = shift;

    if (@_) {
        $self->{-jdbh} = shift;
    }

    return $self->{-jdbh};
}

sub state
{
    my $self = shift;

    if (@_) {
        $self->{-state} = shift;
    }

    return $self->{-state};
}

sub errors
{
    my $self = shift;

    if (@_) {
        my @errors = @_;
        if (defined $errors[0]) {
            $self->param('errors', \@errors);
        } else {
            $self->param('errors', undef);
        }
    }

    return $self->param('errors');
}

sub warnings
{
    my $self = shift;

    if (@_) {
        my @warnings = @_;
        if (defined $warnings[0]) {
            $self->param('warnings', \@warnings);
        } else {
            $self->param('warnings', undef);
        }
    }

    return $self->param('warnings');
}

sub jaspar_db_connect
{
    my $self = shift;

    my $dbh = TFBS::DB::JASPAR5->connect(
        "dbi:mysql:" . JASPAR_DB_NAME . ":" . JASPAR_DB_HOST,
        JASPAR_DB_USER,
        JASPAR_DB_PASS
    );

    if (!$dbh) {
        $self->_error("Could not connect to JASPAR database " . JASPAR_DB_NAME);
    }

    $self->jdbh($dbh);
}

sub fetch_tf_set
{
    my ($self, %args) = @_;

    my %matrix_args = %args;

    unless ($matrix_args{-matrixtype}) {
        $matrix_args{-matrixtype} = 'PFM';
    }

    #printf STDERR "fetch_tf_set: matrix_args = \n"
    #    . Data::Dumper::Dumper(%matrix_args) . "\n";

    my $jdbh = $self->jdbh();

    my $matrix_set = $jdbh->get_MatrixSet(%matrix_args);

    my $tf_set = OPOSSUM::TFSet->new(-matrix_set => $matrix_set);

    return $tf_set;
}

#
# Create a local working file. Used to store input gene and TF IDs.
#
sub create_local_working_file
{
    my ($self, $dir, $name, $lines) = @_;

    my $filename = "$dir/$name";

    unless ($filename =~ /\.txt$/) {
        $filename .= '.txt';
    }

    unless (open(FH, ">$filename")) {
        $self->_error("Could not create local working file $filename - $!");
        return undef;
    }

    foreach my $line (@$lines) {
        print FH "$line\n";
    }

    close(FH);

    return $filename;
}

#
# Low-level error routine. Add latest error to internal error list and
# output to stderr (log file).
#
# Errors are now stored in the CGI::Application params rather than the
# state object.
# DJA 2012/10/17
#
#
sub _error
{
    my ($self, $error) = @_;

    return unless $error;

    #
    # Don't carp errors to log file yet. Do it in high level error routine
    # so that related errors are written in the correct order, see comments
    # above.
    # DJA 2012/10/17
    #
    #carp "\nERROR: $error\n";

    #
    # Now put new errors on the front of the list so errors will be written
    # in the correct order. More general, higher level routine's errors are
    # added latter but should be written earlier.
    # DJA 2012/10/17
    #
    my @errors;
    push @errors, $error;

    my $cur_errors = $self->errors();
    if ($cur_errors) {
        push @errors, @$cur_errors;
    }

    $self->param('errors', \@errors);

    return @errors ? \@errors : undef;
}

#
# Low-level warning routine. Add latest warning to internal warning list and
# output to stderr (log file).
#
sub _warning
{
    my ($self, $warning) = @_;

    return unless $warning;

    carp "\nWarning: $warning\n";

    my @warnings;
    push @warnings, $warning;

    my $cur_warnings = $self->warnings();
    if ($cur_warnings) {
        push @warnings, @$cur_warnings;
    }

    $self->param('warnings', \@warnings);

    return @warnings ? \@warnings : undef;
}

sub _clean_tempfiles
{
    my $self = shift;

    my @tempfiles = glob(ABS_HTDOCS_TMP_PATH . "/*");
    foreach my $file (@tempfiles) {
        unlink $file if -M $file > REMOVE_TEMPFILES_OLDER_THAN;
    }
}

sub _clean_resultfiles
{
    my $self = shift;

    my @files = glob(ABS_HTDOCS_RESULTS_PATH . "/*");

    foreach my $file (@files) {
        if (-M $file > REMOVE_RESULTFILES_OLDER_THAN) {
            if (-d $file) {
                # remove entire tree if directory
                rmtree($file, 0, 0);
            } elsif (-f $file) {
                # unlink if file
                unlink($file);
            }
        }
    }
}

sub _time
{
    my @months   = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second,     $minute,    $hour,
        $dayOfMonth, $month,     $yearOffset,
        $dayOfWeek,  $dayOfYear, $daylightSavings
    ) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime =
        "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    #print STDERR $theTime;
    return $theTime;
}

sub _session_tmp_file
{
    my $sid = shift;

    #return sprintf("%s/$sid", ABS_HTDOCS_TMP_PATH);
    return sprintf("%s/$sid", OPOSSUM_TMP_PATH);
}

1;
