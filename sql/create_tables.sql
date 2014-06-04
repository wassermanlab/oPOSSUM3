--
-- Create table statements for the oPOSSUM DB
--

--
-- The db_info table contains meta information about the database itself
--
drop table if exists db_info;
create table db_info (
    build_date          datetime not NULL,
    species             varchar(32) not NULL,
    latin_name          varchar(32) not NULL,
    assembly            varchar(16) not NULL,
    ensembl_db          varchar(64) not NULL,
    ucsc_db             varchar(64) not NULL,
    ucsc_cons_table     varchar(32) not NULL,
    min_threshold       real(4,3) unsigned not NULL,
    max_upstream_bp     int(5) not NULL,
    min_cr_length       int(4) not NULL,
    tax_group           varchar(64) NOT NULL,
    min_ic              tinyint(2) unsigned NOT NULL,
    has_operon          tinyint(1) unsigned default 0
) engine = MyISAM;

--
-- The genes table contains information about genes
--
-- NOTE: the start and end define the entire gene region including any
-- upstream portion. In a +ve strand gene the start coordinate defines the
-- 5' end of the upstream part of the gene region and in a -ve strand gene,
-- the end defines the 5' end of the upstream part of the gene region. In
-- both cases, the tss defines the 5' most transcription start position.
-- Therefore it is always true that start < tss < end regardless of strand.
--
-- NOTE: The chr field seems excessively large but when processing fly,
-- a group of genes were located on a chromosome called
-- 'dmel_mitochondrion_genome'.
--
drop table if exists genes;
create table genes (
    gene_id             int(10) unsigned not NULL,
    ensembl_id          varchar(32) not NULL,
    symbol              varchar(32) default NULL,
--  description         varchar(1024) default NULL,
    biotype             varchar(40) default NULL,
    chr                 varchar(32) not NULL,
    start               int(10) unsigned not NULL,
    end                 int(10) unsigned not NULL,
    tss                 int(10) unsigned not NULL,
    strand              tinyint(1) not NULL,
    primary key (gene_id),
    index (ensembl_id),
    index (symbol)
) engine = MyISAM;

drop table if exists operons;
create table operons (
    operon_id               int(10) unsigned not NULL,
    symbol                  varchar(32) not NULL,
    gene_id                 int(10) not NULL,
    primary key (operon_id,gene_id),
    index (operon_id),
    index (symbol)
) engine = MyISAM;

drop table if exists promoters;
create table promoters (
    gene_id                 int(10) unsigned not NULL,
    tss                     int(10) unsigned not NULL,
    ensembl_transcript_id   varchar(25),
    index (gene_id)
) engine = MyISAM;

drop table if exists exons;
create table exons (
    gene_id       int(10) unsigned not NULL,
    start         int(10) unsigned not NULL,
    end           int(10) unsigned not NULL,
    index (gene_id)
) engine = MyISAM;

drop table if exists conserved_regions;
create table conserved_regions (
    gene_id                 int(10) unsigned not NULL,
    conservation_level      tinyint(1) unsigned not NULL,
    start                   int(10) unsigned not NULL,
    end                     int(10) unsigned not NULL,
    conservation            real(4,3) unsigned not NULL,
    gc_content              real(4,3) unsigned not NULL,
    index (gene_id, conservation_level)
) engine = MyISAM;

drop table if exists conserved_region_lengths;
create table conserved_region_lengths (
    gene_id                 int(10) unsigned not NULL,
    conservation_level      tinyint(1) unsigned not NULL,
    search_region_level     tinyint(1) unsigned not NULL,
    length                  int(6) unsigned not NULL,
    primary key (conservation_level, search_region_level, gene_id)
) engine = MyISAM;

drop table if exists conserved_tfbss;
create table conserved_tfbss (
    gene_id             int(10) unsigned not NULL,
    tf_id               varchar(16) not NULL,
    start               int(10) unsigned not NULL,
    end                 int(10) unsigned not NULL,
    strand              tinyint(1) not NULL,
    score               real(5,3) not NULL,
    rel_score           real(4,3) not NULL,
    seq                 varchar(40),
    conservation_level  tinyint(1) unsigned not NULL,
    conservation        real(4,3) unsigned not NULL,
    index (gene_id, tf_id),
    index (gene_id, conservation_level)
) engine = MyISAM;

drop table if exists tfbs_counts;
create table tfbs_counts (
    gene_id             int(10) unsigned not NULL,
    tf_id               varchar(16) not NULL,
    conservation_level  tinyint(1) unsigned not NULL,
    threshold_level     tinyint(1) unsigned not NULL,
    search_region_level tinyint(1) unsigned not NULL,
    count               int(6) unsigned,
    primary key (conservation_level, threshold_level, search_region_level, tf_id, gene_id)
) engine = MyISAM;

drop table if exists tfbs_cluster_counts;
create table tfbs_cluster_counts (
    gene_id             int(10) unsigned not NULL,
    cluster_id          int(10) unsigned not NULL,
    conservation_level  tinyint(1) unsigned not NULL,
    threshold_level     tinyint(1) unsigned not NULL,
    search_region_level tinyint(1) unsigned not NULL,
    count               int(6) unsigned,
    sum_length              int(10) unsigned,
    primary key (gene_id, cluster_id, conservation_level, search_region_level, threshold_level),
    index (cluster_id)
) engine = MyISAM;

drop table if exists conservation_levels;
create table conservation_levels (
    level               tinyint(1) unsigned not NULL,
    min_conservation    real(4,3) unsigned not NULL
) engine = MyISAM;

drop table if exists search_region_levels;
create table search_region_levels (
    level               tinyint(1) unsigned not NULL,
    upstream_bp         smallint(5) unsigned,
    downstream_bp       smallint(5) unsigned 
) engine = MyISAM;

drop table if exists threshold_levels;
create table threshold_levels (
    level               tinyint(1) unsigned not NULL,
    threshold           real(4,3) unsigned not NULL
) engine = MyISAM;

drop table if exists external_gene_id_types;
create table external_gene_id_types (
    id_type     tinyint(2) unsigned not NULL,
    name        varchar(40) not NULL,
    dblink_name varchar(40) not NULL
) engine = MyISAM;

drop table if exists external_gene_ids;
create table external_gene_ids (
    gene_id         int(10) unsigned not NULL,
    id_type         tinyint(2) unsigned not NULL,
    external_id     varchar(40) not NULL,
    primary key (gene_id, id_type, external_id),
    index (id_type, external_id)
) engine = MyISAM;
