# coding: utf-8

# SQLAlchemy Mappings generated  from Rfam 14.0 using sqlacodegen (https://pypi.org/project/sqlacodegen/)
# sqlacodegen mysql+mysqlconnector://rfamro@mysql-rfam-public.ebi.ac.uk:4497/Rfam --outfile rfam_db.py

import sqlalchemy.orm
import sqlalchemy.engine.url
from sqlalchemy import (
    CHAR,
    Column,
    DECIMAL,
    Date,
    DateTime,
    Enum,
    Float,
    ForeignKey,
    Index,
    String,
    TIMESTAMP,
    Table,
    Text,
    text,
)
from sqlalchemy.dialects.mysql import (
    BIGINT,
    INTEGER,
    LONGBLOB,
    LONGTEXT,
    MEDIUMINT,
    MEDIUMTEXT,
    TINYINT,
    TINYTEXT,
    VARCHAR,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from dotenv import load_dotenv, find_dotenv
import os

Base = declarative_base()
metadata = Base.metadata


def rfam_session():
    """
    Create a new connection to the Rfam database.
    :return:
    """

    dotenv_path = find_dotenv()
    load_dotenv(dotenv_path)

    # Set up a connection to a local docker copy of Rfam if appropriate environment variables have been set, otherwise
    # connect to the public Rfam database

    database = os.environ.get("MYSQL_DATABASE", default="Rfam")
    host = os.environ.get("MYSQL_HOST", default="mysql-rfam-public.ebi.ac.uk")
    port = os.environ.get("MYSQL_PORT", default=4497)
    username = os.environ.get("MYSQL_USER", default="rfamro")
    password = os.environ.get("MYSQL_PASSWORD")

    # Set up the connection to the MySQL server and create a session
    database_url = sqlalchemy.engine.url.URL(
        "mysql+mysqlconnector",
        username=username,
        password=password,
        host=host,
        port=port,
        database=database,
    )
    engine = sqlalchemy.create_engine(database_url, encoding="utf8")
    SessionMaker = sqlalchemy.orm.sessionmaker(bind=engine)
    session = SessionMaker()

    return session


class Author(Base):
    __tablename__ = "author"

    author_id = Column(INTEGER(11), primary_key=True, unique=True)
    name = Column(String(20), nullable=False)
    last_name = Column(String(50))
    initials = Column(String(4))
    orcid = Column(String(19))
    synonyms = Column(String(100))


class Clan(Base):
    __tablename__ = "clan"

    clan_acc = Column(String(7), primary_key=True, unique=True)
    id = Column(String(40))
    previous_id = Column(TINYTEXT)
    description = Column(String(100))
    author = Column(TINYTEXT)
    comment = Column(LONGTEXT)
    created = Column(DateTime, nullable=False)
    updated = Column(TIMESTAMP, nullable=False)

    family = relationship("Family", secondary="clan_membership")


class DbVersion(Base):
    __tablename__ = "db_version"

    rfam_release = Column(Float(4, True), primary_key=True)
    rfam_release_date = Column(DateTime, nullable=False)
    number_families = Column(INTEGER(10), nullable=False)
    embl_release = Column(TINYTEXT, nullable=False)
    genome_collection_date = Column(DateTime)
    refseq_version = Column(INTEGER(11))
    pdb_date = Column(DateTime)
    infernal_version = Column(String(45))


t_dead_clan = Table(
    "dead_clan",
    metadata,
    Column("clan_acc", String(7), nullable=False, unique=True, server_default=text("''")),
    Column("clan_id", String(40), nullable=False),
    Column("comment", MEDIUMTEXT),
    Column("forward_to", String(7)),
    Column("user", TINYTEXT, nullable=False),
)


t_dead_family = Table(
    "dead_family",
    metadata,
    Column("rfam_acc", String(7), nullable=False, unique=True, server_default=text("''")),
    Column("rfam_id", String(40), nullable=False),
    Column("comment", MEDIUMTEXT),
    Column("forward_to", String(7)),
    Column("title", String(150)),
    Column("user", TINYTEXT, nullable=False),
)


class EnsemblName(Base):
    __tablename__ = "ensembl_names"

    insdc = Column(String(50), primary_key=True, server_default=text("''"))
    ensembl = Column(String(50))


t_family_author = Table(
    "family_author",
    metadata,
    Column("rfam_acc", String(7), nullable=False),
    Column("author_id", INTEGER(11), nullable=False),
    Column("desc_order", INTEGER(4), nullable=False),
)


class Genome(Base):
    __tablename__ = "genome"

    upid = Column(String(20), primary_key=True, unique=True, server_default=text("''"))
    assembly_acc = Column(String(20))
    assembly_version = Column(INTEGER(4))
    wgs_acc = Column(String(20))
    wgs_version = Column(INTEGER(4))
    assembly_name = Column(String(100))
    assembly_level = Column(Enum("contig", "chromosome", "scaffold", "complete-genome"))
    study_ref = Column(String(20))
    description = Column(MEDIUMTEXT)
    total_length = Column(BIGINT(20))
    ungapped_length = Column(BIGINT(20))
    circular = Column(TINYINT(3))
    ncbi_id = Column(INTEGER(10), nullable=False, index=True)
    scientific_name = Column(String(300))
    common_name = Column(String(200))
    kingdom = Column(String(50))
    num_rfam_regions = Column(INTEGER(10))
    num_families = Column(INTEGER(10))
    is_reference = Column(TINYINT(1), nullable=False, server_default=text("'1'"))
    is_representative = Column(TINYINT(1), nullable=False, server_default=text("'0'"))
    created = Column(DateTime, nullable=False)
    updated = Column(TIMESTAMP, nullable=False, server_default=text("CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP"))


t_genseq = Table(
    "genseq",
    metadata,
    Column("upid", String(20), nullable=False, index=True, server_default=text("''")),
    Column("rfamseq_acc", String(20), nullable=False, index=True, server_default=text("''")),
    Column("chromosome_name", String(100)),
    Column("chromosome_type", String(100)),
    Column("version", String(6)),
)


class Keyword(Base):
    __tablename__ = "keywords"
    __table_args__ = (
        Index("rfam_kw_idx", "description", "rfam_general", "literature", "wiki", "pdb_mappings", "clan_info"),
    )

    rfam_acc = Column(String(7), primary_key=True, server_default=text("''"))
    rfam_id = Column(String(40))
    description = Column(String(100), server_default=text("'NULL'"))
    rfam_general = Column(LONGTEXT)
    literature = Column(LONGTEXT)
    wiki = Column(LONGTEXT)
    pdb_mappings = Column(LONGTEXT)
    clan_info = Column(LONGTEXT)


class LiteratureReference(Base):
    __tablename__ = "literature_reference"

    pmid = Column(INTEGER(10), primary_key=True)
    title = Column(TINYTEXT)
    author = Column(Text)
    journal = Column(TINYTEXT)


class Motif(Base):
    __tablename__ = "motif"

    motif_acc = Column(String(7), primary_key=True)
    motif_id = Column(String(40), index=True)
    description = Column(String(75))
    author = Column(TINYTEXT)
    seed_source = Column(TINYTEXT)
    gathering_cutoff = Column(Float(5, True))
    trusted_cutoff = Column(Float(5, True))
    noise_cutoff = Column(Float(5, True))
    cmbuild = Column(TINYTEXT)
    cmcalibrate = Column(TINYTEXT)
    type = Column(String(50))
    num_seed = Column(BIGINT(20))
    average_id = Column(Float(5, True))
    average_sqlen = Column(Float(7, True))
    ecmli_lambda = Column(Float(10, True))
    ecmli_mu = Column(Float(10, True))
    ecmli_cal_db = Column(MEDIUMINT(9), server_default=text("'0'"))
    ecmli_cal_hits = Column(MEDIUMINT(9), server_default=text("'0'"))
    maxl = Column(MEDIUMINT(9), server_default=text("'0'"))
    clen = Column(MEDIUMINT(9), server_default=text("'0'"))
    match_pair_node = Column(TINYINT(1), server_default=text("'0'"))
    hmm_tau = Column(Float(10, True))
    hmm_lambda = Column(Float(10, True))
    wiki = Column(String(80))
    created = Column(DateTime, nullable=False)
    updated = Column(TIMESTAMP, nullable=False, server_default=text("CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP"))


class Pdb(Base):
    __tablename__ = "pdb"

    pdb_id = Column(String(4), primary_key=True)
    keywords = Column(TINYTEXT)
    title = Column(MEDIUMTEXT)
    date = Column(TINYTEXT)
    resolution = Column(DECIMAL(5, 2), server_default=text("'0.00'"))
    method = Column(TINYTEXT)
    author = Column(MEDIUMTEXT)


class Refseq(Base):
    __tablename__ = "refseq"

    refseq_acc = Column(String(14), primary_key=True)
    description = Column(MEDIUMTEXT)
    species = Column(MEDIUMTEXT)
    ncbi_taxid = Column(INTEGER(10))


t_rnacentral_matches = Table(
    "rnacentral_matches",
    metadata,
    Column("rfamseq_acc", String(20), nullable=False, index=True, server_default=text("''")),
    Column("seq_start", BIGINT(19), nullable=False, index=True, server_default=text("'0'")),
    Column("seq_end", BIGINT(19), nullable=False, index=True),
    Column("md5", String(32), nullable=False),
    Column("rnacentral_id", String(13), index=True),
)


class Taxonomy(Base):
    __tablename__ = "taxonomy"

    ncbi_id = Column(INTEGER(10), primary_key=True, server_default=text("'0'"))
    species = Column(String(100), nullable=False, index=True, server_default=text("''"))
    tax_string = Column(MEDIUMTEXT)
    tree_display_name = Column(String(100))
    align_display_name = Column(String(50))


t_taxonomy_websearch = Table(
    "taxonomy_websearch",
    metadata,
    Column("ncbi_id", INTEGER(10), index=True, server_default=text("'0'")),
    Column("species", String(100), index=True),
    Column("rgt", INTEGER(10), index=True),
    Column("taxonomy", MEDIUMTEXT),
    Column("lft", INTEGER(10), index=True),
    Column("parent", INTEGER(10), index=True),
    Column("level", String(200), index=True),
    Column("minimal", TINYINT(1), nullable=False, index=True, server_default=text("'0'")),
    Column("rank", String(100)),
)


class Version(Base):
    __tablename__ = "version"

    rfam_release = Column(Float(4, True), primary_key=True)
    rfam_release_date = Column(Date, nullable=False)
    number_families = Column(INTEGER(10), nullable=False)
    embl_release = Column(TINYTEXT, nullable=False)


class Wikitext(Base):
    __tablename__ = "wikitext"

    auto_wiki = Column(INTEGER(10), primary_key=True)
    title = Column(VARCHAR(150), nullable=False, unique=True)


t_clan_database_link = Table(
    "clan_database_link",
    metadata,
    Column("clan_acc", ForeignKey("clan.clan_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("db_id", TINYTEXT, nullable=False),
    Column("comment", TINYTEXT),
    Column("db_link", TINYTEXT, nullable=False),
    Column("other_params", TINYTEXT),
)


t_clan_literature_reference = Table(
    "clan_literature_reference",
    metadata,
    Column("clan_acc", ForeignKey("clan.clan_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("pmid", ForeignKey("literature_reference.pmid", ondelete="CASCADE"), nullable=False, index=True),
    Column("comment", TINYTEXT),
    Column("order_added", TINYINT(3)),
)


class Family(Base):
    __tablename__ = "family"

    rfam_acc = Column(String(7), primary_key=True, unique=True)
    rfam_id = Column(String(40), nullable=False, index=True)
    auto_wiki = Column(ForeignKey("wikitext.auto_wiki"), nullable=False, index=True)
    description = Column(String(75))
    author = Column(TINYTEXT)
    seed_source = Column(TINYTEXT)
    gathering_cutoff = Column(Float(5, True))
    trusted_cutoff = Column(Float(5, True))
    noise_cutoff = Column(Float(5, True))
    comment = Column(LONGTEXT)
    previous_id = Column(TINYTEXT)
    cmbuild = Column(TINYTEXT)
    cmcalibrate = Column(TINYTEXT)
    cmsearch = Column(TINYTEXT)
    num_seed = Column(BIGINT(20))
    num_full = Column(BIGINT(20))
    num_genome_seq = Column(BIGINT(20))
    num_refseq = Column(BIGINT(20))
    type = Column(String(50))
    structure_source = Column(TINYTEXT)
    number_of_species = Column(BIGINT(20))
    number_3d_structures = Column(INTEGER(11))
    tax_seed = Column(MEDIUMTEXT)
    ecmli_lambda = Column(Float(10, True))
    ecmli_mu = Column(Float(10, True))
    ecmli_cal_db = Column(MEDIUMINT(9), server_default=text("'0'"))
    ecmli_cal_hits = Column(MEDIUMINT(9), server_default=text("'0'"))
    maxl = Column(MEDIUMINT(9), server_default=text("'0'"))
    clen = Column(MEDIUMINT(9), server_default=text("'0'"))
    match_pair_node = Column(TINYINT(1), server_default=text("'0'"))
    hmm_tau = Column(Float(10, True))
    hmm_lambda = Column(Float(10, True))
    created = Column(DateTime, nullable=False)
    updated = Column(TIMESTAMP, nullable=False, server_default=text("CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP"))

    wikitext = relationship("Wikitext")


t_motif_database_link = Table(
    "motif_database_link",
    metadata,
    Column(
        "motif_acc", ForeignKey("motif.motif_acc", ondelete="CASCADE", onupdate="CASCADE"), nullable=False, index=True
    ),
    Column("db_id", TINYTEXT, nullable=False),
    Column("comment", TINYTEXT),
    Column("db_link", TINYTEXT, nullable=False),
    Column("other_params", TINYTEXT),
)


t_motif_file = Table(
    "motif_file",
    metadata,
    Column("motif_acc", ForeignKey("motif.motif_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("seed", LONGBLOB, nullable=False),
    Column("cm", LONGBLOB, nullable=False),
)


t_motif_literature = Table(
    "motif_literature",
    metadata,
    Column(
        "motif_acc",
        ForeignKey("motif_old.motif_acc", ondelete="CASCADE", onupdate="CASCADE"),
        nullable=False,
        index=True,
    ),
    Column(
        "pmid",
        ForeignKey("literature_reference.pmid", ondelete="CASCADE", onupdate="CASCADE"),
        nullable=False,
        index=True,
    ),
    Column("comment", TINYTEXT),
    Column("order_added", TINYINT(3)),
)


t_motif_pdb = Table(
    "motif_pdb",
    metadata,
    Column(
        "motif_acc",
        ForeignKey("motif_old.motif_acc", ondelete="CASCADE", onupdate="CASCADE"),
        nullable=False,
        index=True,
    ),
    Column("pdb_id", String(4), nullable=False, index=True),
    Column("chain", String(4)),
    Column("pdb_start", MEDIUMINT(9)),
    Column("pdb_end", MEDIUMINT(9)),
)


class Rfamseq(Base):
    __tablename__ = "rfamseq"

    rfamseq_acc = Column(String(20), primary_key=True, unique=True, server_default=text("''"))
    accession = Column(String(15), nullable=False)
    version = Column(INTEGER(6), nullable=False, index=True)
    ncbi_id = Column(ForeignKey("taxonomy.ncbi_id", ondelete="CASCADE"), nullable=False, index=True)
    mol_type = Column(
        Enum(
            "protein",
            "genomic DNA",
            "DNA",
            "ss-DNA",
            "RNA",
            "genomic RNA",
            "ds-RNA",
            "ss-cRNA",
            "ss-RNA",
            "mRNA",
            "tRNA",
            "rRNA",
            "snoRNA",
            "snRNA",
            "scRNA",
            "pre-RNA",
            "other RNA",
            "other DNA",
            "unassigned DNA",
            "unassigned RNA",
            "viral cRNA",
            "cRNA",
            "transcribed RNA",
        ),
        nullable=False,
    )
    length = Column(INTEGER(10), server_default=text("'0'"))
    description = Column(String(250), nullable=False, server_default=text("''"))
    previous_acc = Column(MEDIUMTEXT)
    source = Column(CHAR(20), nullable=False)

    ncbi = relationship("Taxonomy")


t_alignment_and_tree = Table(
    "alignment_and_tree",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("type", Enum("seed", "seedTax", "genome", "genomeTax"), nullable=False),
    Column("alignment", LONGBLOB),
    Column("tree", LONGBLOB),
    Column("treemethod", TINYTEXT),
    Column("average_length", Float(7, True)),
    Column("percent_id", Float(5, True)),
    Column("number_of_sequences", BIGINT(20)),
)


t_clan_membership = Table(
    "clan_membership",
    metadata,
    Column("clan_acc", ForeignKey("clan.clan_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, unique=True),
)


t_database_link = Table(
    "database_link",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("db_id", TINYTEXT, nullable=False),
    Column("comment", TINYTEXT),
    Column("db_link", TINYTEXT, nullable=False),
    Column("other_params", TINYTEXT),
)


t_family_literature_reference = Table(
    "family_literature_reference",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("pmid", ForeignKey("literature_reference.pmid", ondelete="CASCADE"), nullable=False, index=True),
    Column("comment", TINYTEXT),
    Column("order_added", TINYINT(3)),
)


t_family_ncbi = Table(
    "family_ncbi",
    metadata,
    Column("ncbi_id", ForeignKey("taxonomy.ncbi_id", ondelete="CASCADE"), nullable=False, index=True),
    Column("rfam_id", String(40)),
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
)


t_features = Table(
    "features",
    metadata,
    Column("rfamseq_acc", ForeignKey("rfamseq.rfamseq_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("database_id", String(50), nullable=False),
    Column("primary_id", String(100), nullable=False),
    Column("secondary_id", String(255)),
    Column("feat_orient", TINYINT(3), nullable=False, server_default=text("'0'")),
    Column("feat_start", BIGINT(19), nullable=False, server_default=text("'0'")),
    Column("feat_end", BIGINT(19), nullable=False, server_default=text("'0'")),
    Column("quaternary_id", String(150)),
)


t_full_region = Table(
    "full_region",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False),
    Column("rfamseq_acc", ForeignKey("rfamseq.rfamseq_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("seq_start", BIGINT(19), nullable=False, server_default=text("'0'")),
    Column("seq_end", BIGINT(19), nullable=False),
    Column("bit_score", Float(7, True), nullable=False, server_default=text("'0.00'")),
    Column("evalue_score", String(15), nullable=False, server_default=text("'0'")),
    Column("cm_start", MEDIUMINT(8), nullable=False),
    Column("cm_end", MEDIUMINT(8), nullable=False),
    Column("truncated", Enum("0", "5", "3", "53"), nullable=False),
    Column("type", Enum("seed", "full"), nullable=False, server_default=text("'full'")),
    Column("is_significant", TINYINT(1), nullable=False, index=True),
    Index("full_region_acc_sign", "rfam_acc", "is_significant"),
)


t_html_alignment = Table(
    "html_alignment",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("type", Enum("seed", "genome", "seedColorstock", "genomeColorstock"), nullable=False, index=True),
    Column("html", LONGBLOB),
    Column("block", INTEGER(6), nullable=False, index=True),
    Column("html_alignmentscol", VARCHAR(45)),
)


t_matches_and_fasta = Table(
    "matches_and_fasta",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc"), nullable=False, index=True),
    Column("match_list", LONGBLOB),
    Column("fasta", LONGBLOB),
    Column("type", Enum("rfamseq", "genome", "refseq"), nullable=False),
)


t_motif_family_stats = Table(
    "motif_family_stats",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc"), nullable=False, index=True),
    Column("motif_acc", ForeignKey("motif_old.motif_acc"), nullable=False, index=True),
    Column("num_hits", INTEGER(11)),
    Column("frac_hits", DECIMAL(4, 3)),
    Column("sum_bits", DECIMAL(12, 3)),
    Column("avg_weight_bits", DECIMAL(12, 3)),
)


t_motif_matches = Table(
    "motif_matches",
    metadata,
    Column(
        "motif_acc",
        ForeignKey("motif_old.motif_acc", ondelete="CASCADE", onupdate="CASCADE"),
        nullable=False,
        index=True,
    ),
    Column(
        "rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE", onupdate="CASCADE"), nullable=False, index=True
    ),
    Column(
        "rfamseq_acc",
        ForeignKey("rfamseq.rfamseq_acc", ondelete="CASCADE", onupdate="CASCADE"),
        nullable=False,
        index=True,
    ),
    Column("rfamseq_start", BIGINT(19)),
    Column("rfamseq_stop", BIGINT(19)),
    Column("query_start", INTEGER(11)),
    Column("query_stop", INTEGER(11)),
    Column("motif_start", INTEGER(11)),
    Column("motif_stop", INTEGER(11)),
    Column("e_value", String(15)),
    Column("bit_score", Float(7, True)),
)


t_motif_ss_image = Table(
    "motif_ss_image",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("motif_acc", ForeignKey("motif_old.motif_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("image", LONGBLOB),
)


t_pdb_full_region = Table(
    "pdb_full_region",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("pdb_id", String(4), nullable=False, index=True),
    Column("chain", VARCHAR(4)),
    Column("pdb_start", MEDIUMINT(8), nullable=False),
    Column("pdb_end", MEDIUMINT(8), nullable=False),
    Column("bit_score", Float(7, True), nullable=False, server_default=text("'0.00'")),
    Column("evalue_score", String(15), nullable=False, server_default=text("'0'")),
    Column("cm_start", MEDIUMINT(8), nullable=False),
    Column("cm_end", MEDIUMINT(8), nullable=False),
    Column("hex_colour", String(6), server_default=text("'NULL'")),
    Column("is_significant", TINYINT(1), nullable=False, index=True, server_default=text("'1'")),
)


t_secondary_structure_image = Table(
    "secondary_structure_image",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column(
        "type",
        Enum(
            "cons",
            "dist",
            "ent",
            "fcbp",
            "cov",
            "disttruc",
            "maxcm",
            "norm",
            "rchie",
            "species",
            "ss",
            "rscape",
            "rscape-cyk",
        ),
        index=True,
    ),
    Column("image", LONGBLOB),
)


t_seed_region = Table(
    "seed_region",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("rfamseq_acc", ForeignKey("rfamseq.rfamseq_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("seq_start", BIGINT(19), nullable=False, server_default=text("'0'")),
    Column("seq_end", BIGINT(19), nullable=False),
)


t_sunburst = Table(
    "sunburst",
    metadata,
    Column("rfam_acc", ForeignKey("family.rfam_acc", ondelete="CASCADE"), nullable=False, index=True),
    Column("data", LONGBLOB, nullable=False),
    Column("type", Enum("rfamseq", "genome", "refseq"), nullable=False),
)
