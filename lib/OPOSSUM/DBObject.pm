=head1 NAME

OPOSSUM::DBObject

=head1 SYNOPSIS

  my $db_id     = $storable_object->db_id();
  my $adaptor   = $storable_object->adaptor();

=head1 DESCRIPTION

This is an oPOSSUM DBObject class. It has been adapted from the
Bio::EnsEMBL::Storable class. All objects which are storable in the
oPOSSUM database should inherit from this class. It provides two
getter/setters: db_id() adaptor().

=head1 METHODS

=cut

use strict;
use warnings;

package OPOSSUM::DBObject;

use OPOSSUM::Utils::Argument qw(rearrange);

=head2 new

  Arg [-ADAPTOR] : OPOSSUM::DBSQL::BaseAdaptor
  Arg [-db_id]    : database internal id
  Example        : none 
  Caller         : internal calls
  Description    : create a new DBObject object 
  Returntype     : OPOSSUM::DBObject
  Exceptions     : Adaptor not a OPOSSUM::DBSQL::BaseAdaptor
  Status         : Stable

=cut

sub new
{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    my ($adaptor, $db_id) = rearrange(['ADAPTOR', 'db_id'], @_);

    if ($adaptor) {
        if (   !ref($adaptor)
            || !$adaptor->isa('OPOSSUM::DBSQL::BaseAdaptor'))
        {
            throw(
                '-ADAPTOR argument must be a OPOSSUM::DBSQL::BaseAdaptor'
            );
        }
    }

    return bless({'db_id' => $db_id, 'adaptor' => $adaptor}, $class);
}

=head2 db_id

  Arg [1]    : int $db_id
  Example    : none
  Description: getter/setter for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store
  Status     : Stable

=cut

sub db_id
{
    my $self = shift;

    $self->{'-db_id'} = shift if (@_);

    return $self->{'-db_id'};
}

=head2 adaptor

  Arg [1]    : OPOSSUM::DBSQL::BaseAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : OPOSSUM::DBSQL::BaseAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store
  Status     : Stable

=cut

sub adaptor
{
    my $self = shift;

    if (@_) {
        my $ad = shift;
        if ($ad
            && (!ref($ad) || !$ad->isa('OPOSSUM::DBSQL::BaseAdaptor')))
        {
            throw(
                'Adaptor argument must be a OPOSSUM::DBSQL::BaseAdaptor'
            );
        }
        $self->{'-adaptor'} = $ad;
    }

    return $self->{'-adaptor'};
}

1;
