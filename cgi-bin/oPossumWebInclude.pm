#
# This module should be included in all the oPossum*Web.pm modules and
# possibly the background perl scripts called by those modules. It
# contains all routines that are common to all the oPossum variants. This
# includes utility functions as well as the common template routines like
# 'errors' and 'warnings'.
#

use oPossumWebOpt;

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

    my $state = $self->state();

    $self->_error($error) if $error;

    my $errors = $state->errors();

    my $error_html;
    foreach my $err (@$errors) {
        chomp $err;
        $err =~ s/\n/<br>/g;
        $error_html .= "$err<br>";
    }

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
        errors           => $error_html,
        var_template     => "error.html"
    };

    my $output = $self->process_template('master.html', $vars);

    $self->state->errors(undef);

    return $output;
}

#
# High-level warning routine. Call low level _warning routine with current
# warning and output all current warnings to HTML warning template.
#
sub warning
{
    my ($self, $warning) = @_;

    my $state = $self->state();

    $self->_warning($warning) if $warning;

    my $warnings = $state->warnings();

    my $warning_html;
    foreach my $warn (@$warnings) {
        chomp $warn;
        $warn =~ s/\n/<br>/g;
        $warning_html .= "$warn<br>";
    }

    my $warnings = $self->warnings;

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
        warnings         => $warning_html,
        var_template     => "warning.html"
    };

    my $output = $self->process_template('master.html', $vars);

    $self->state->errors(undef);

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
    my ($self, $errors) = @_;

    if (@_) {
        $self->state->errors(shift);
    }

    return $self->state->errors();
}

sub warnings
{
    my ($self, $warnings) = @_;

    if (@_) {
        $self->state->warnings(shift);
    }

    return $self->state->warnings();
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
sub _error
{
    my ($self, $error) = @_;

    carp "\nERROR: $error\n";

    my $errors = $self->state->errors();

    push @$errors, "$error";

    $self->state->errors($errors);
}

#
# Low-level warning routine. Add latest warning to internal warning list and
# output to stderr (log file).
#
sub _warning
{
    my ($self, $warning) = @_;

    carp "\nWarning: $warning\n";

    my $warnings = $self->state->warnings();

    push @$warnings, "$warning";

    $self->state->warnings($warnings);
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

    return sprintf("%s/$sid", ABS_HTDOCS_TMP_PATH);
}

1;
