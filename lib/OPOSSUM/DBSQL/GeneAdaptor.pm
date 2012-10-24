=head1 NAME

OPOSSUM::DBSQL::GeneAdaptor - Adaptor for MySQL queries to retrieve and
store Gene objects.

=head1 SYNOPSIS

$ga = $db_adaptor->get_GeneAdaptor();

=head1 DESCRIPTION

The genes table of the oPOSSUM database stores general gene information.

=head1 AUTHOR

 David Arenillas
 Wasserman Lab
 Centre for Molecular Medicine and Therapeutics
 University of British Columbia

 E-mail: dave@cmmt.ubc.ca

=head1 MODIFICATIONS

AK 2010-09-07
- handling of biotype parameter
- fetch calls using external gene ids now able to return multiple
  opossum genes/ids as there can be many to many relationships btw
  external ids and opossum gene ids
AK 2010-09-27
- added synonyms ending with *_list for *_genes and *_ids
- fetch_biotypes added
- fetch_random_genes added
- reorganization of functions to reduce code redundancies
DJA 2010-09-29
- added _fetch_gene_ids_by_external_id_helper based on Andrew's
  _fetch_by_external_id_helper as it is more efficient to do this then
  fetch all the genes and then get their IDs
- changed methods which call _fetch_by_external_id_helper just to retrieve
  gene IDs to call the _fetch_gene_ids_by_external_id_helper method instead
- added filter for distinct genes to _fetch_by_external_id_helper
DJA 2012-06-12
- added fetch_gene_ids_by_any_id and _fetch_gene_ids_by_any_helper to fetch
  gene IDs for any/all types of input gene ID(s) including symbols, Ensembl IDs
  or external IDs without having to know the gene ID type. The input gene IDs
  can even be of mixed types.

=cut

package OPOSSUM::DBSQL::GeneAdaptor;

use strict;

use Carp;

use OPOSSUM::DBSQL::BaseAdaptor;
use OPOSSUM::Gene;

use vars '@ISA';
@ISA = qw(OPOSSUM::DBSQL::BaseAdaptor);

=head2 new

 Title   : new
 Usage   : $ga = OPOSSUM::DBSQL::GeneAdaptor->new($db_adaptor);
 Function: Construct a new GeneAdaptor object
 Args    : An OPOSSUM::DBSQL::DBAdaptor object
 Returns : a new OPOSSUM::DBSQL::GeneAdaptor object

=cut

sub new
{
    my ($class, @args) = @_;

    $class = ref $class || $class;

    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 fetch_gene_count

 Title   : fetch_gene_count
 Usage   : $genes = $ga->fetch_gene_count($where);
 Function: Generic fetch method. Fetch gene object(s) from the DB
           with the given where clause.
 Args    : Optionally, a where clause.
 Returns : Either an OPOSSUM::Gene object or a reference to an array of
           OPOSSUM::Gene objects depending on whether query returns 1 or
           more rows.

=cut

sub fetch_gene_count
{
    my ($self, $where) = @_;

    my $sql = qq{select count(*) from genes};
	if ($where) {
        unless ($where =~ /^\s*where /) {
            $sql .= " where";
        }
		$sql .= " $where";
	}
	
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch gene count:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch gene count:\n$sql\n" . $self->errstr;
        return;
    }

    my $count;
    unless (($count) = $sth->fetchrow_array()) {
        carp "Error fetching gene count:\n$sql\n" . $self->errstr;
        return;
    }

    $sth->finish();

    return $count;
}

=head2 fetch_where

 Title   : fetch_where
 Usage   : $genes = $ga->fetch_where($where);
 Function: Generic fetch method. Fetch gene object(s) from the DB
           with the given where clause.
 Args    : Optionally, a where clause.
 Returns : Either an OPOSSUM::Gene object or a reference to an array of
           OPOSSUM::Gene objects depending on whether query returns 1 or
           more rows.

=cut

sub fetch_where
{
    my ($self, $where) = @_;

    my $sql =
        qq{select gene_id,
            ensembl_id,
			symbol,
            biotype,
            chr,
            start,
            end,
            tss,
            strand
        from genes};

	if ($where) {
        unless ($where =~ /^\s*where / or $where =~ /^\s*order /) {
            $sql .= " where";
        }
		$sql .= " $where";
	}
	
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch genes:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch genes:\n$sql\n" . $self->errstr;
        return;
    }

    my $row_count = 0;
    my @genes;
    while (my @row = $sth->fetchrow_array) {
        $row_count++;

        my $gene = OPOSSUM::Gene->new(
            -adaptor        => $self,
            -id             => $row[0],
            -ensembl_id     => $row[1],
            -symbol         => $row[2],
            -biotype        => $row[3],
            -chr            => $row[4],
            -start          => $row[5],
            -end            => $row[6],
            -tss            => $row[7],
            -strand         => $row[8]
        );

        push @genes, $gene;
    }
    $sth->finish;

    return @genes if wantarray();

    #
    # This is actually a really bad idea as it puts the onus on the calling
    # program to determine if it got an array ref or a single gene. Change
    # the code in the individual callers (e.g. fetch_by_id to return the
    # expected result.
    # DJA 2012/10/23
    #
    #if ($row_count) {
    #    if ($row_count > 1) {
    #        return \@genes;
    #    }
    #    return $genes[0];
    #}

    return @genes ? \@genes : undef;
}

=head2 fetch_gene_ids

 Title   : fetch_gene_ids
 Usage   : $gids = $ga->fetch_gene_ids($where);
 Function: Fetch list of gene IDs from the DB.
 Args    : Optionally a where clause.
 Returns : Reference to a list of internal gene IDs. If no where clause
           is provided, returns all gene IDs in the database.

=cut

sub fetch_gene_ids
{
    my ($self, $where) = @_;

    my $sql = "select gene_id from genes";

	if ($where) {
        unless ($where =~ /^\s*where / or $where =~ /^\s*order /) {
            $sql .= " where";
        }
		$sql .= " $where";
	}

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error fetching gene IDs\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error fetching gene IDs\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_gene_id_list

 Title   : fetch_gene_id_list
 Usage   : $gids = $ga->fetch_gene_id_list($where);
 Function: Fetch list of gene IDs from the DB. Synonym of fetch_gene_ids.
 Args    : Optionally a where clause.
 Returns : Reference to a list of internal gene IDs. If no where clause
           is provided, returns all gene IDs in the database.

=cut

sub fetch_gene_id_list
{
	my ($self, $where) = @_;
	
	return $self->fetch_gene_ids($where);
}

=head2 fetch_gene_id_fields

 Title   : fetch_gene_id_fields
 Usage   : $ids = $ga->fetch_gene_id_fields($field, $where);
 Function: Fetch list of IDs of type specified by $field from the DB.
 Args    : Optionally a field type and where clause.
 Returns : Reference to a list of IDs. If no field type is specified,
           returns internal gene IDs, otherwise the IDs stored in the
           given field. If no where clause is provided, returns all IDs in
           the database.

=cut

sub fetch_gene_id_fields
{
    my ($self, $id_field, $where) = @_;

    if (!$id_field) {
        return $self->fetch_gene_ids;
    }

    my $sql = "select distinct $id_field from genes";
	if ($where) {
        unless ($where =~ /^\s*where /) {
            $sql .= " where";
        }
		$sql .= " $where";
	}

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error fetching gene $id_field fields\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error fetching gene $id_field fields\n" . $self->errstr;
        return;
    }

    my @ids;
    while (my ($id) = $sth->fetchrow_array) {
        push @ids, $id;
    }
    $sth->finish;

    return @ids ? \@ids : undef;
}

=head2 fetch_by_id

 Title   : fetch_by_id
 Usage   : $gene = $ga->fetch_by_id($id);
 Function: Fetch a gene object from the DB using its ID.
 Args    : The unique internal gene ID.
 Returns : An OPOSSUM::Gene object.

=cut

sub fetch_by_id
{
    my ($self, $id) = @_;

    my $where = "where gene_id = $id";

    my $genes = $self->fetch_where($where);

    if ($genes && $genes->[0]) {
        return $genes->[0];
    }

    return undef;
}

=head2 fetch_by_gene_id

 Title   : fetch_by_gene_id
 Usage   : $gene = $ga->fetch_by_gene_id($id);
 Function: Fetch a gene object from the DB using its ID. Synonym of fetch_by_id.
 Args    : The unique internal gene ID.
 Returns : An OPOSSUM::Gene object.

=cut

sub fetch_by_gene_id
{
    my ($self, $id) = @_;

    return $self->fetch_by_id($id);
}

=head2 fetch_by_ensembl_id

 Title   : fetch_by_ensembl_id
 Usage   : $gene = $ga->fetch_by_ensembl_id($id);
 Function: Fetch a gene object from the DB using its Ensembl ID.
 Args    : The unique gene Ensembl ID.
 Returns : An OPOSSUM::Gene object.

=cut

sub fetch_by_ensembl_id
{
    my ($self, $ensembl_id) = @_;

    return if !$ensembl_id;

    my $where = "where ensembl_id = '$ensembl_id'";

    my $genes = $self->fetch_where($where);

    if ($genes && $genes->[0]) {
        return $genes->[0];
    }

    return undef;
}

=head2 fetch_by_symbol

 Title   : fetch_by_symbol
 Usage   : $gene = $ga->fetch_by_symbol($symbol);
 Function: Fetch a gene object from the DB using its gene symbol.
 Args    : The unique gene symbol.
 Returns : An OPOSSUM::Gene object.

=cut
sub fetch_by_symbol
{
    my ($self, $symbol) = @_;

    return if !$symbol;

    my $where = "where symbol = '$symbol'";

    my $genes = $self->fetch_where($where);

    if ($genes && $genes->[0]) {
        return $genes->[0];
    }

    return undef;
}

=head2 fetch_by_gene_ids

 Title   : fetch_by_gene_ids
 Usage   : $genes = $ga->fetch_by_gene_ids($id_list);
 Function: Fetch a list of gene objects from the DB according to a list
           of gene IDs.
 Args    : A reference to a list of unique internal gene IDs.
 Returns : A reference to a list of OPOSSUM::Gene objects.

=cut

sub fetch_by_gene_ids
{
	my ($self, $gene_ids) = @_;
	
	my $where = "gene_id in (" . join(",", @$gene_ids) . ")";

	return $self->fetch_where($where);
}

=head2 fetch_by_gene_id_list

 Title   : fetch_by_gene_id_list
 Usage   : $genes = $ga->fetch_by_gene_id_list($id_list);
 Function: Fetch a list of gene objects from the DB according to a list
           of gene IDs. Synonym of fetch_by_gene_ids.
 Args    : A reference to a list of unique internal gene IDs.
 Returns : A reference to a list of OPOSSUM::Gene objects.

=cut

sub fetch_by_gene_id_list
{
	my ($self, $gene_ids) = @_;
	
	return $self->fetch_by_gene_ids($gene_ids);
}

=head2 fetch_by_ensembl_ids

 Title   : fetch_by_ensembl_ids
 Usage   : $genes = $ga->fetch_by_ensembl_ids($id_list);
 Function: Fetch a list of gene objects from the DB according to a list
           of Ensembl IDs.
 Args    : A reference to a list of unique Ensembl IDs.
 Returns : A reference to a list of OPOSSUM::Gene objects.

=cut

sub fetch_by_ensembl_ids
{
	my ($self, $ensembl_ids) = @_;
	
	my $where = "ensembl_id in ('" . join("','", @$ensembl_ids) . "')";

	return $self->fetch_where($where);
}

=head2 fetch_by_ensembl_id_list

 Title   : fetch_by_ensembl_id_list
 Usage   : $genes = $ga->fetch_by_ensembl_id_list($id_list);
 Function: Fetch a list of gene objects from the DB according to a list
           of Ensembl IDs. Synonym of fetch_by_ensembl_ids.
 Args    : A reference to a list of unique Ensembl IDs.
 Returns : A reference to a list of OPOSSUM::Gene objects.

=cut

sub fetch_by_ensembl_id_list
{
	my ($self, $ensembl_ids) = @_;
	
	return $self->fetch_by_ensembl_ids($ensembl_ids);
}

=head2 fetch_by_symbols

 Title   : fetch_by_symbols
 Usage   : $genes = $ga->fetch_by_symbols($symbols);
 Function: Fetch a list of gene objects from the DB according to a list
           of gene symbols.
 Args    : A reference to a list of unique gene symbols.
 Returns : A reference to a list of OPOSSUM::Gene objects.

=cut

sub fetch_by_symbols
{
	my ($self, $symbols) = @_;
	
	my $where = "symbol in ('" . join("','", @$symbols) . "')";

	return $self->fetch_where($where);
}

=head2 fetch_by_symbol_list

 Title   : fetch_by_symbol_list
 Usage   : $genes = $ga->fetch_by_symbol_list($id_list);
 Function: Fetch a list of gene objects from the DB according to a list
           of gene symbols. Synonym of fetch_by_symbols.
 Args    : A reference to a list of unique gene symbols.
 Returns : A reference to a list of OPOSSUM::Gene objects.

=cut

sub fetch_by_symbol_list
{
	my ($self, $symbols) = @_;
	
	return $self->fetch_by_symbols($symbols);
}


#
# External ID methods
#

=head2 _fetch_by_external_id_helper

 Title   : _fetch_by_external_id_helper
 Usage   : $genes = $ga->_fetch_by_external_id_helper(
               -ext_id              => $ext_id,
               -ext_ids             => $ext_ids,
               -ext_id_type         => $ext_id_type,
               -mapped_ext_ids      => $mapped_ext_ids,
               -unmapped_ext_ids    => $unmapped_ext_ids,
               -gene_id_ensembl_ids => $gene_id_ensembl_ids,
               -gene_id_ext_ids     => $gene_id_ext_ids,
               -ext_id_gene_ids     => $ext_id_gene_ids
           );
 Function: Fetch gene objects from the DB by external gene ID(s).
           Internal method only.
 Args    : A single external ID or listref of external IDs,
           OPTIONAL external ID type,
           OPTIONAL output arguments:
               mapped_ext_ids - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               unmapped_ext_IDs - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               gene_id_ensembl_ids - a hashref which is set to the mapping
                    between opossum gene IDs and Ensembl IDs (1-to-1)
               gene_id_ext_ids - a hashref which is set to the mapping from
                    opossum gene IDs to external IDs (1-to-many),
               ext_id_gene_ids mapping - a hashref which is set to the
                    mapping from external IDs to opossum gene IDs (1-to-many)
 Returns : A listref of OPOSSUM::Gene objects.

=cut

sub _fetch_by_external_id_helper
{
    my ($self, %args) = @_;
	
	my $ext_id              = $args{-ext_id};
    my $ext_ids             = $args{-ext_ids};
	my $ext_id_type         = $args{-ext_id_type};
	
	unless ($ext_id or $ext_ids) {
		carp "Must specify external gene id(s)\n";
		return;
	}
	
	if ($ext_id and $ext_ids) {
		carp "Both ext_id and ext_ids specified\n";
	}

    #
    # The following are optional return variables which are updated by the
    # routine.
    #
    my $mapped_ext_ids      = $args{-mapped_ext_ids};       # listref
    my $unmapped_ext_ids    = $args{-unmapped_ext_ids};     # listref
    my $gene_id_ensembl_ids = $args{-gene_id_ensembl_ids};  # hashref
    my $gene_id_ext_ids     = $args{-gene_id_ext_ids};      # hashref
    my $ext_id_gene_ids     = $args{-ext_id_gene_ids};      # hashref

    push @$ext_ids, $ext_id if $ext_id;
	
    my $sql =
        qq{select g.gene_id,
			g.ensembl_id,
			g.symbol,
			g.biotype,
			g.chr,
			g.start,
			g.end,
			g.tss,
			g.strand
		from genes g, external_gene_ids xgi
		where g.gene_id = xgi.gene_id
            and xgi.external_id = ?};
			
    if ($ext_id_type) {
        $sql .= " and xgi.id_type = $ext_id_type";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch gene:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @genes;
    my %gene_included;
    my %ext_id_included;
    foreach my $ext_id (@$ext_ids) {
        if (!$sth->execute($ext_ids)) {
            carp "Error executing fetch gene:\n$sql\nwith external_id $ext_id\n"
                . $self->errstr;
            next;
        }

        #
        # There can be a many-to-many mapping of oPOSSUM gene IDs and
        # external gene IDs.
        #
        while (my @row = $sth->fetchrow_array()) {
            my $gene_id = $row[0];
            my $ensid   = $row[1];

            #
            # There is a 1-to-1 mapping of oPOSSUM gene IDs
            # and Ensembl IDs.
            #
            unless ($gene_included{$gene_id}) {
                $gene_included{$gene_id} = 1;

                push @genes, OPOSSUM::Gene->new(
                    -adaptor        => $self,
                    -id             => $row[0],
                    -ensembl_id     => $row[1],
                    -symbol         => $row[2],
                    -biotype        => $row[3],
                    -chr            => $row[4],
                    -start          => $row[5],
                    -end            => $row[6],
                    -tss            => $row[7],
                    -strand         => $row[8]
                );

                $gene_id_ensembl_ids->{$gene_id} = $ensid;
            }

            #
            # Each external gene ID can have multiple oPOSSUM gene IDs
            # mapped to it.
            #
            unless ($ext_id_included{$ext_id}) {
                push @$mapped_ext_ids, $ext_id;

                $ext_id_included{$ext_id} = 1;
            }

            push @{$ext_id_gene_ids->{$ext_id}}, $gene_id;
            push @{$gene_id_ext_ids->{$gene_id}}, $ext_id;
        }
    }
    $sth->finish();

    #
    # Determine which of the input Ensembl / external gene IDs were not mapped
    # to an oPOSSUM gene ID (missing).
    #
    foreach my $ext_id (@$ext_ids) {
        unless ($ext_id_included{$ext_id}) {
            push @$unmapped_ext_ids, $ext_id;
        }
    }

    return @genes ? \@genes : undef;
}

=head2 _fetch_gene_ids_by_external_id_helper

 Title   : _fetch_gene_ids_by_external_id_helper
 Usage   : $gene_ids = $ga->_fetch_gene_ids_by_external_id_helper(
               -ext_id              => $ext_id,
               -ext_ids             => $ext_ids,
               -ext_id_type         => $ext_id_type,
               -mapped_ext_ids      => $mapped_ext_ids,
               -unmapped_ext_ids    => $unmapped_ext_ids,
               -gene_id_ensembl_ids => $gene_id_ensembl_ids,
               -gene_id_ext_ids     => $gene_id_ext_ids,
               -ext_id_gene_ids     => $ext_id_gene_ids
           );
 Function: Fetch gene IDs from the DB by external gene ID(s).
           Internal method only.
 Args    : A single external ID or listref of external IDs,
           OPTIONAL external ID type,
           OPTIONAL output arguments:
               mapped_ext_ids - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               unmapped_ext_IDs - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               gene_id_ensembl_ids - a hashref which is set to the mapping
                    between opossum gene IDs and Ensembl IDs (1-to-1)
               gene_id_ext_ids - a hashref which is set to the mapping from
                    opossum gene IDs to external IDs (1-to-many),
               ext_id_gene_ids mapping - a hashref which is set to the
                    mapping from external IDs to opossum gene IDs (1-to-many)
 Returns : A listref of gene IDs.

=cut

sub _fetch_gene_ids_by_external_id_helper
{
    my ($self, %args) = @_;
	
	my $ext_id          = $args{-ext_id};
    my $ext_ids         = $args{-ext_ids};
	my $ext_id_type     = $args{-ext_id_type};
	
	unless ($ext_id or $ext_ids) {
		carp "Must specify external gene id(s)\n";
		return;
	}
	
	if ($ext_id and $ext_ids) {
		carp "Both ext_id and ext_ids specified\n";
	}
	
    #
    # The following are optional return variables which are updated by the
    # routine.
    #
    my $mapped_ext_ids      = $args{-mapped_ext_ids};       # listref
    my $unmapped_ext_ids    = $args{-unmapped_ext_ids};     # listref
    my $gene_id_ensembl_ids = $args{-gene_id_ensembl_ids};  # hashref
    my $gene_id_ext_ids     = $args{-gene_id_ext_ids};      # hashref
    my $ext_id_gene_ids     = $args{-ext_id_gene_ids};      # hashref

    push @$ext_ids, $ext_id if $ext_id;
	
    my $sql =
        qq{select g.gene_id, g.ensembl_id
           from genes g, external_gene_ids xgi
           where g.gene_id = xgi.gene_id
           and xgi.external_id = ?};
			
    if ($ext_id_type) {
        $sql .= " and xgi.id_type = $ext_id_type";
    }

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch gene ID:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @gene_ids;
    my %gene_included;
    my %ext_id_included;
    foreach my $ext_id (@$ext_ids) {
        if (!$sth->execute($ext_id)) {
            carp "Error executing fetch gene ID:\n$sql\n"
                . "with external_id = $ext_id"
                . $self->errstr;
            next;
        }

        #
        # There can be a many-to-many mapping of oPOSSUM gene IDs and
        # external gene IDs.
        #
        while (my @row = $sth->fetchrow_array) {
            my $gene_id = $row[0];
            my $ensid   = $row[1];

            #
            # There is a 1-to-1 mapping of oPOSSUM gene IDs
            # and Ensembl IDs.
            #
            unless ($gene_included{$gene_id}) {
                $gene_included{$gene_id} = 1;

                push @gene_ids, $gene_id;

                $gene_id_ensembl_ids->{$gene_id} = $ensid;
            }

            #
            # Each external gene ID can have multiple oPOSSUM gene IDs
            # mapped to it.
            #
            unless ($ext_id_included{$ext_id}) {
                push @$mapped_ext_ids, $ext_id;

                $ext_id_included{$ext_id} = 1;
            }

            push @{$ext_id_gene_ids->{$ext_id}}, $gene_id;
            push @{$gene_id_ext_ids->{$gene_id}}, $ext_id;
        }
    }
    $sth->finish;

    #
    # Determine which of the input Ensembl / external gene IDs were not mapped
    # to an oPOSSUM gene ID (missing).
    #
    foreach my $ext_id (@$ext_ids) {
        unless ($ext_id_included{$ext_id}) {
            push @$unmapped_ext_ids, $ext_id;
        }
    }

    return @gene_ids ? \@gene_ids : undef;
}

=head2 _fetch_gene_ids_by_any_id_helper

 Title   : _fetch_gene_ids_by_any_id_helper
 Usage   : $gene_ids = $ga->_fetch_gene_ids_by_any_id_helper(
               -ext_id              => $ext_id,
               -ext_ids             => $ext_ids,
               -mapped_ext_ids      => $mapped_ext_ids,
               -unmapped_ext_ids    => $unmapped_ext_ids,
               -gene_id_ensembl_ids => $gene_id_ensembl_ids,
               -gene_id_ext_ids     => $gene_id_ext_ids,
               -ext_id_gene_ids     => $ext_id_gene_ids
           );
 Function: Fetch oPOSSUM gene IDs from the DB by all possible gene
           ID/name/symbol mappings, i.e. the ensembl_id column or symbol
           column of the gene record or an external id mapped to the
           gene record that matches the input gene ID/name/symbol(s).
           Internal method only.
 Args    : A single value or listref of gene ID/name/symbol(s),
           OPTIONAL output arguments:
               mapped_ext_ids - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               unmapped_ext_IDs - a listref which is set to the list of
                    external IDs which are NOT mapped to opossum genes,
               gene_id_ensembl_ids - a hashref which is set to the mapping
                    between opossum gene IDs and Ensembl IDs (1-to-1)
               gene_id_ext_ids - a hashref which is set to the mapping from
                    opossum gene IDs to external IDs (1-to-many),
               ext_id_gene_ids mapping - a hashref which is set to the
                    mapping from external IDs to opossum gene IDs (1-to-many)
 Returns : A listref of gene IDs.

=cut

sub _fetch_gene_ids_by_any_id_helper
{
    my ($self, %args) = @_;
	
	my $ext_id          = $args{-ext_id};
    my $ext_ids         = $args{-ext_ids};
	
	unless ($ext_id or $ext_ids) {
		carp "Must specify external gene id(s)\n";
		return;
	}
	
	if ($ext_id and $ext_ids) {
		carp "Both ext_id and ext_ids specified\n";
	}
	
    #
    # The following are optional return variables which are updated by the
    # routine.
    #
    my $mapped_ext_ids      = $args{-mapped_ext_ids};       # listref
    my $unmapped_ext_ids    = $args{-unmapped_ext_ids};     # listref
    my $gene_id_ensembl_ids = $args{-gene_id_ensembl_ids};  # hashref
    my $gene_id_ext_ids     = $args{-gene_id_ext_ids};      # hashref
    my $ext_id_gene_ids     = $args{-ext_id_gene_ids};      # hashref

    push @$ext_ids, $ext_id if $ext_id;
	
    my $sql =
        qq{select gene_id, ensembl_id from genes
           where ensembl_id = ? or symbol = ? or gene_id = (
               select gene_id from external_gene_ids where external_id = ?)};
			
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch gene ID:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @gene_ids;
    my %gene_included;
    my %ext_id_included;
    foreach my $ext_id (@$ext_ids) {
        if (!$sth->execute($ext_id)) {
            carp "Error executing fetch gene ID:\n$sql\n"
                . "with external_id = $ext_id"
                . $self->errstr;
            next;
        }

        #
        # There can be a many-to-many mapping of oPOSSUM gene IDs and
        # external gene IDs.
        #
        while (my @row = $sth->fetchrow_array) {
            my $gene_id = $row[0];
            my $ensid   = $row[1];

            #
            # There is a 1-to-1 mapping of oPOSSUM gene IDs
            # and Ensembl IDs.
            #
            unless ($gene_included{$gene_id}) {
                $gene_included{$gene_id} = 1;

                push @gene_ids, $gene_id;

                $gene_id_ensembl_ids->{$gene_id} = $ensid;
            }

            #
            # Each external gene ID can have multiple oPOSSUM gene IDs
            # mapped to it.
            #
            unless ($ext_id_included{$ext_id}) {
                push @$mapped_ext_ids, $ext_id;

                $ext_id_included{$ext_id} = 1;
            }

            #
            # Note in the case that the input external gene IDs were themselves
            # Ensembl IDs, they will be mapped to the oPOSSUM gene IDs below
            # and in the $gene_id_ensembl_ids hash.
            #
            push @{$ext_id_gene_ids->{$ext_id}}, $gene_id;
            push @{$gene_id_ext_ids->{$gene_id}}, $ext_id;
        }
    }
    $sth->finish;

    #
    # Determine which of the input Ensembl / external gene IDs were not mapped
    # to an oPOSSUM gene ID (missing).
    #
    foreach my $ext_id (@$ext_ids) {
        unless ($ext_id_included{$ext_id}) {
            push @$unmapped_ext_ids, $ext_id;
        }
    }

    return @gene_ids ? \@gene_ids : undef;
}

=head2 fetch_by_external_id

 Title   : fetch_by_external_id
 Usage   : $genes = $ga->fetch_by_external_id($ext_id, $ext_id_type);
 Function: Fetch gene objects from the DB by an external gene ID.
 Args    : The external ID and external ID type.
           OPTIONAL return arguments (see *_helper method for details)
 Returns : A listref of OPOSSUM::Gene objects.

=cut

sub fetch_by_external_id
{
    my ($self, $ext_id, $ext_id_type, %return_args) = @_;

    return $self->_fetch_by_external_id_helper(
		-ext_id      => $ext_id,
		-ext_id_type => $ext_id_type,
        %return_args
	);
}

=head2 fetch_by_external_ids

 Title   : fetch_by_external_ids
 Usage   : $genes = $ga->fetch_by_external_ids($ext_ids, $ext_id_type);
 Function: Fetch gene objects from the DB by external gene IDs.
 Args    : The external IDs and external ID type.
           OPTIONAL return arguments (see *_helper method for details)
 Returns : A listref of OPOSSUM::Gene objects.

=cut

sub fetch_by_external_ids
{
    my ($self, $ext_ids, $ext_id_type, %return_args) = @_;

    return $self->_fetch_by_external_id_helper(
		-ext_ids     => $ext_ids,
		-ext_id_type => $ext_id_type,
        %return_args
	);
}

=head2 fetch_by_external_id_list

 Title   : fetch_by_external_id_list
 Usage   : $genes = $ga->fetch_by_external_id_list($ext_ids, $ext_id_type);
 Function: Fetch gene objects from the DB by external gene IDs.
			Synonym of fetch_by_external_ids.
 Args    : The external ID and external ID type.
           OPTIONAL return arguments (see *_helper method for details)
 Returns : A listref of OPOSSUM::Gene objects.

=cut

sub fetch_by_external_id_list
{
    my ($self, $ext_ids, $ext_id_type, %return_args) = @_;

    return $self->_fetch_by_external_id_helper(
		-ext_ids     => $ext_ids,
		-ext_id_type => $ext_id_type,
        %return_args
	);
}

#
# The method name implies fetching a single ID but the method returns a
# listref of IDs which could cause problems. Best to comment it out and
# force users to use fetch_gene_ids_by_external_id method instead.
#
#=head2 fetch_gene_id_by_external_id
#
# Title   : fetch_gene_id_by_external_id
# Usage   : $geneid = $ga->fetch_gene_id_by_external_id(
#               $ext_id, $ext_id_type
#           );
# Function: Fetch gene ID from the DB corresponding to an external
#           gene ID.
# Args    : The external ID and external ID type.
# Returns : A listref of oPOSSUM gene IDs
#
#=cut
#
#sub fetch_gene_id_by_external_id
#{
#    my ($self, $ext_id, $ext_id_type) = @_;
#
#    my $genes = $self->_fetch_by_external_id_helper(
#        -ext_id => $ext_id,
#        -ext_id_type => $ext_id_type
#    );
#
#    return if !$genes;
#
#    my @geneids;
#    foreach my $gene (@$genes) {
#        push @geneids, $gene->id;
#    }
#
#    return \@geneids;
#}

=head2 fetch_gene_ids_by_external_id

 Title   : fetch_gene_ids_by_external_id
 Usage   : $gene_ids = $ga->fetch_gene_ids_by_external_id(
               $ext_id, $ext_id_type
           );
 Function: Fetch gene IDs from the DB corresponding to an external
           gene ID.
 Args    : An external gene ID,
           OPTIONAL external ID type.
           OPTIONAL return arguments (see *_helper method for details)
 Returns : A listref of oPOSSUM gene IDs

=cut

sub fetch_gene_ids_by_external_id
{
	my ($self, $ext_id, $ext_id_type, %return_args) = @_;
	
    return $self->_fetch_gene_ids_by_external_id_helper(
        -ext_id      => $ext_id,
        -ext_id_type => $ext_id_type,
        %return_args
    );
}

=head2 fetch_gene_ids_by_any_id

 Title   : fetch_gene_ids_by_any_id
 Usage   : $gene_ids = $ga->fetch_gene_ids_by_any_id($id);
 Function: Fetch gene IDs from the DB corresponding to any type of gene
           ID/name/symbol.
 Args    : Some kind of gene ID/name/symbol,
           OPTIONAL return arguments (see *_helper method for details)
 Returns : A listref of oPOSSUM gene IDs

=cut

sub fetch_gene_ids_by_any_id
{
	my ($self, $id, %return_args) = @_;
	
    return $self->_fetch_gene_ids_by_any_id_helper(
        -ext_id => $id,
        %return_args
    );
}

=head2 fetch_gene_ids_by_external_ids

 Title   : fetch_gene_ids_by_external_ids
 Usage   : $gene_ids = $ga->fetch_gene_ids_by_external_ids(
               $ext_ids, $ext_id_type
           );
 Function: Fetch a list of gene IDs from the DB corresponding to a
           list of external gene IDs. Duplicate gene IDs are
           filtered (included only once in the returned list).
 Args    : The external ID type and a listref of external gene IDs.
           OPTIONAL return arguments (see *_helper method for details)
 Returns : A listref of gene IDs.

=cut

sub fetch_gene_ids_by_external_ids
{
    my ($self, $ext_ids, $ext_id_type, %return_args) = @_;

    return if !$ext_ids || !$ext_ids->[0];

	return $self->_fetch_gene_ids_by_external_id_helper(
		-ext_ids     => $ext_ids,
		-ext_id_type => $ext_id_type,
        %return_args
	);
}

=head2 fetch_gene_ids_by_ensembl_ids

 Title   : fetch_gene_ids_by_ensembl_ids
 Usage   : $gene_ids = $ga->fetch_gene_ids_by_ensembl_ids(
               $ensembl_ids,
               -mapped_ensembl_ids      => $mapped_ensids,
               -unmapped_ensembl_ids    => $unmapped_ensids,
               -gene_id_ensembl_ids     => $gene_id_ensids,
               -ensembl_id_gene_ids     => $ensid_gene_ids
           );
 Function: Fetch gene IDs from the DB by Ensembl ID(s).
 Args    : A listref of Ensembl IDs,
           OPTIONAL output arguments:
               mapped_ensids - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               unmapped_ensids - a listref which is set to the list of
                    external IDs which are actually mapped to opossum genes,
               gene_id_ensids - a hashref which is set to the mapping
                    between opossum gene IDs and Ensembl IDs (1-to-1)
               ensid_gene_ids mapping - a hashref which is set to the
                    mapping from external IDs to opossum gene IDs (1-to-1)
 Returns : A listref of gene IDs.

=cut

sub fetch_gene_ids_by_ensembl_ids
{
    my ($self, $ensids, %return_args) = @_;
	
	unless ($ensids) {
		carp "Must specify Ensembl id(s)\n";
		return;
	}

    my $mapped_ensids   = $return_args{-mapped_ensembl_ids};
    my $unmapped_ensids = $return_args{-unmapped_ensembl_ids};
    my $gene_id_ensids  = $return_args{-gene_id_ensembl_ids};
    my $ensid_gene_ids  = $return_args{-ensembl_id_gene_ids};
	
    my $sql = qq{select gene_id from genes where ensembl_id = ?};
			
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch gene ID:\n$sql\n"
            . $self->errstr;
        return;
    }

    my @gene_ids;
    my %ensid_included;
    foreach my $ensid (@$ensids) {
        if (!$sth->execute($ensid)) {
            carp "Error executing fetch gene ID:\n$sql\n"
                . "with ensembl_id = $ensid"
                . $self->errstr;
            next;
        }

        #
        # There is a 1-to-1 mapping of oPOSSUM gene IDs
        # and Ensembl IDs.
        #
        if (my @row = $sth->fetchrow_array) {
            my $gene_id = $row[0];

            push @gene_ids, $gene_id;

            $gene_id_ensids->{$gene_id} = $ensid;
            $ensid_gene_ids->{$ensid} = $gene_id;

            push @$mapped_ensids, $ensid;

            $ensid_included{$ensid} = 1;
        }
    }
    $sth->finish;

    #
    # Determine which of the input Ensembl IDs were not mapped
    # to an oPOSSUM gene ID (missing).
    #
    foreach my $ensid (@$ensids) {
        unless ($ensid_included{$ensid}) {
            push @$unmapped_ensids, $ensid;
        }
    }

    return @gene_ids ? \@gene_ids : undef;
}

=head2 fetch_gene_id_list_by_external_id_list

 Title   : fetch_gene_id_list_by_external_id_list
 Usage   : $geneids = $ga->fetch_gene_id_list_by_external_id_list(
               $ext_ids, $ext_id_type
           );
 Function: Fetch a list of gene IDs from the DB corresponding to a
           list of external gene IDs. Duplicate gene IDs are
           filtered (included only once in the returned list).
		   Synonym of fetch_gene_ids_by_external_ids
 Args    : The external ID type and a listref of external gene IDs.
 Returns : A listref of gene IDs.

=cut

sub fetch_gene_id_list_by_external_id_list
{
    my ($self, $ext_ids, $ext_id_type) = @_;

    return $self->fetch_gene_ids_by_external_ids($ext_ids, $ext_id_type);
}


=head2 fetch_random_genes

 Title   : fetch_random_random_genes
 Usage   : $gids = $ga->fetch_random_genes(
               -num_genes => 5,
               -biotype => 'protein_coding' | -biotypes => $biotypes
           );
 Function: Fetches a random list of OPOSSUM::GENE objects.
 Args    : The number of gene ids to retrieve (integer) and the
           biotype (string) or biotypes (listref of string) of genes to be
           fetched.
 Returns : A listref of OPOSSUM::Gene objects.

=cut

sub fetch_random_genes
{
	my ($self, %args) = @_;
	
	my $num_genes   = $args{-num_genes};
	my $biotype     = $args{-biotype};
	my $biotypes    = $args{-biotypes};
	
	if ($biotype and $biotypes) {
		carp "Both biotype and biotypes specified\n";
		push @$biotypes, $biotype;
		$biotype = undef;
	}
	
    my $where;
	if ($biotypes) {
		$where = "where biotype in ('" . join("','", @$biotypes) . "')";
	} elsif ($biotype) {
        unless ($biotype =~ /^all$/i) {
            $where = "where biotype = '$biotype'";
        }
	}

	$where .= " order by rand() limit $num_genes";
	
	return $self->fetch_where($where);
}

=head2 fetch_random_gene_list

 Title   : fetch_random_random_gene_list
 Usage   : $gids = $ga->fetch_random_gene_list(
               -num_genes => 5,
               -biotype => 'protein_coding' | -biotypes => $biotypes
           );
 Function: Fetches a random list of OPOSSUM::Gene objects. Synonym of
           fetch_random_genes.
 Args    : The number of gene ids to retrieve (integer) and the
           biotype (string) or biotypes (listref of string) of genes to be
           fetched.
 Returns : A listref of OPOSSUM::Gene objects.

=cut

sub fetch_random_gene_list
{
	my ($self, %args) = @_;

	return $self->fetch_random_genes(%args);
}

=head2 fetch_random_gene_ids

 Title   : fetch_random_random_gene_ids
 Usage   : $gids = $ga->fetch_random_gene_ids($num);
 Function: Fetches a random list of gene IDs.
 Args    : An integer - number of gene ids to retrieve.
 Returns : A listref of gene IDs.

=cut

sub fetch_random_gene_ids
{
	my ($self, %args) = @_;
	
	my $num_genes   = $args{-num_genes};
	my $biotype     = $args{-biotype};
	my $biotypes    = $args{-biotypes};
	
	if ($biotype and $biotypes) {
		carp "Both biotype and biotypes specified\n";
		push @$biotypes, $biotype;
		$biotype = undef;
	}
	
    my $where;
	if ($biotypes) {
		$where = "where biotype in ('" . join("','", @$biotypes) . "')";
	} elsif ($biotype) {
        unless ($biotype =~ /^all$/i) {
            $where = "where biotype = '$biotype'";
        }
	}

	$where .= " order by rand() limit $num_genes";
	
	return $self->fetch_gene_ids($where);
}

=head2 fetch_random_gene_id_list

 Title   : fetch_random_random_gene_id_list
 Usage   : $gids = $ga->fetch_random_gene_id_list(%args);
 Function: Fetches a random list of gene IDs. Synonym of fetch_random_gene_ids.
 Args    : Number of gene ids and gene biotypes.
 Returns : A listref of gene IDs.

=cut

sub fetch_random_gene_id_list
{
	my ($self, %args) = @_;

	return $self->fetch_random_gene_ids(%args);
}

=head2 fetch_biotypes

 Title   : fetch_biotypes
 Usage   : $biotypes = $ga->fetch_biotypes($where);
 Function: Fetches the list of gene biotypes available.
 Args    : [Optional] where clause string.
 Returns : A listref of gene biotype strings

=cut

sub fetch_biotypes
{
	my ($self, $where) = @_;
	
    my $sql = qq{select distinct(biotype) from genes};
	if ($where) {
        unless ($where =~ /^\s*where /) {
            $where = "where $where";
        }
		$sql .= " $where";
	}
	
    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing fetch biotypes:\n$sql\n" . $self->errstr;
        return;
    }

    if (!$sth->execute) {
        carp "Error executing fetch biotypes:\n$sql\n" . $self->errstr;
        return;
    }

    my @biotypes;
    while (my ($biotype) = $sth->fetchrow_array()) {
        push @biotypes, $biotype;
    }
    $sth->finish();
	
	return @biotypes ? \@biotypes : undef;
}

=head2 fetch_biotype_list

 Title   : fetch_biotype_list
 Usage   : $biotypes = $ga->fetch_biotype_list($where);
 Function: Fetches the list of gene biotypes available. Synonym of
           fetch_biotypes.
 Args    : [Optional] where clause string.
 Returns : A listref of gene biotype strings

=cut

sub fetch_biotype_list
{
	my ($self, $where) = @_;

	return $self->fetch_biotypes($where);
}

=head2 store

 Title   : store
 Usage   : $id = $ga->store($gene);
 Function: Store gene in the database.
 Args    : The gene (OPOSSUM::Gene) to store.
 Returns : A database ID of the newly stored gene.

=cut

sub store
{
    my ($self, $gene) = @_;

    if (!ref $gene || !$gene->isa('OPOSSUM::Gene')) {
        carp "Not an OPOSSUM::Gene object";
        return;
    }

    my $sql
        = qq{insert into genes
		    (ensembl_id, symbol, biotype, chr, start, end, tss, strand)
		    values (?,?,?,?,?,?,?,?)};

    my $sth = $self->prepare($sql);
    if (!$sth) {
        carp "Error preparing insert gene statement\n" . $self->errstr;
        return;
    }

    if (!$sth->execute(
            $gene->ensembl_id, $gene->symbol, $gene->biotype,
            $gene->chr, $gene->start, $gene->end, $gene->tss, $gene->strand)
        )
    {
        carp "Error inserting gene\n" . $self->errstr;
        return;
    }
    $sth->finish;

    return $sth->{'mysql_insertid'};
}

1;
