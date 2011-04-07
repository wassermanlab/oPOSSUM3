=head1 NAME

OPOSSUM::Web::State - Object for keeping track of the state information of
the OPOSSUM web application.

=head1 SYNOPSIS

    use OPOSSUM::Web::State;

    my $state;
    my $sid = $q->param('sid');
    if ($sid) {
        #
        # If session ID is already defined, load the state from file
        #
        my $filename = ABS_TMP_PATH . "/" . $sid;
        $state = OPOSSUM::Web::State->load(__Fn => $filename);
    } else {
        #
        # Create a new session with a unique session ID (in this example the
        # process ID plus current date/time
        #
        $sid = $$ . time;
        my $filename = ABS_TMP_PATH . "/" . $sid;
        $state = OPOSSUM::Web::State->new(__Fn => $filename, -sid => $sid);

        # Set some parameters
        $state->sid($sid);
        $state->species1('human');
        $state->species2('mouse');
    }

    # Save the state to a file
    $state->dumper->Purity(1);
    $state->dumper->Deepcopy(1);
    $state->commit();

=head1 DESCRIPTION

OPOSSUM::Web::State is an object used for the purpose of keeping track of
the web application's state between pages. Inherits from
Persistence::Object::Simple to save and load state information to/from a flat
file between succesive CGI calls based on a unique session ID. It uses the
AUTOLOAD facility to define it's members/methods dynamically.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package OPOSSUM::Web::State;

use strict;

use vars qw/@ISA $AUTOLOAD/;

use Carp;
use Persistence::Object::Simple;

@ISA = qw/Persistence::Object::Simple/;


sub AUTOLOAD
{
    my $self = shift;

    my $type = ref($self) || croak "$self is not an object";

    my $name = $AUTOLOAD;
    $name =~ s/.*://;    # strip fully-qualified portion

    if (@_) {
        $self->{-$name} = shift;
    }

    return $self->{-$name};
}

sub tf_info
{
    my ($self, $source, $collection, $tax_group) = @_;

    if (@_) {
        $self->{-tf_info}->{$source}->{$collection}->{$tax_group} = shift;
    }

    return $self->{-tf_info}->{$source}->{$collection}->{$tax_group};
}

#sub get_tf_info_by_src_col_tax_group
#{
#    my ($self, $src, $col, $tax_group) = @_;
#
#    #
#    # No source, collection or tax group provided.
#    # Return master TF list.
#    #
#    if (!$src && !$col && !$tax_group) {
#        return $self->tf_info();
#    }
#
#    #
#    # Create source/collection/tax group TF info if not already
#    #
#    if (!$self->{-src_col_tax_group_tf_info}) {
#        $self->_create_src_col_tax_group_tf_info();
#    }
#
#    #
#    # Source, collection and tax group all provided.
#    # Return source/collection/tax group TF info.
#    #
#    if ($src && $col && $tax_group) {
#        return $self->{-src_col_tax_group_tf_info}->{$src}->{$col}
#            ->{$tax_group};
#    }
#
#    #
#    # One or more of source, collection and tax group provided. Return
#    # corresponding subset of TF list.
#    #
#    my @sources;
#    if ($src) {
#        # Source provided
#        push @sources, $src;
#    } else {
#        # Use all sources
#        @sources = @{$self->state->tf_sources};
#    }
#
#    my @collections;
#    if ($col) {
#        # Collection provided
#        push @collections, $col;
#    } else {
#        # Use all collections
#        @collections = @{$self->state->tf_collections()};
#    }
#
#    my @tax_groups;
#    if ($tax_group) {
#        # Tax group provided
#        push @tax_groups, $tax_group;
#    } else {
#        # Use all tax groups
#        @tax_groups = @{$self->state->tf_tax_groups()};
#    }
#
#    my @tf_info;
#    foreach my $src (@sources) {
#        foreach my $col (@collections) {
#            foreach my $tax_group (@tax_groups) {
#                if ($self->{src_col_tax_group_tf_info}->{$src}) {
#                    if ($self->{src_col_tax_group_tf_info}->{$src}->{$col}) {
#                        if ($self->{src_col_tax_group_tf_info}->{$src}->{$col}
#                            ->{$tax_group}
#                        )
#                        {
#                            push @tf_info,
#                                $self->{src_col_tax_group_tf_info}->{$src}
#                                ->{$col}->{$tax_group};
#                        }
#                    }
#                }
#            }
#        }
#    }
#
#    return @tf_info ? \@tf_info : undef;
#}
#
#sub _create_src_col_tax_group_tf_info
#{
#    my $self = shift;
#
#    my $tfs = $self->{-tf_info};
#
#    return if !$tfs;
#
#    my %sources;
#    my %collections;
#    my %tax_groups;
#    my $src_col_tax_group_tfs;
#    foreach my $tf (@$tfs) {
#        my $src         = $tf->external_db;
#        my $col         = $tf->collection;
#        my $tax_group   = $tf->tax_group;
#
#        $sources{$src}          = 1;
#        $collections{$col}      = 1;
#        $tax_groups{$tax_group} = 1;
#
#        push @{$src_col_tax_group_tfs->{$src || 'unknown'}
#            ->{$col || 'unknown'}
#            ->{$tax_group || 'unknown'}},
#            $tf;
#    }
#
#    $self->{-src_col_tax_group_tf_info} = $src_col_tax_group_tfs;
#
#    $self->{-tf_sources}     = [keys %sources];
#    $self->{-tf_collections} = [keys %collections];
#    $self->{-tf_tax_groups}  = [keys %tax_groups];
#}

1;
