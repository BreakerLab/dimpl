# -*- coding: utf-8 -*-
import logging
import os
import shutil
import gzip
import numpy as np
import pandas as pd
import time
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
from urllib.request import urlopen
from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC
from bx.intervals.intersection import Interval, IntervalTree
from src.data.rfam_db import rfam_session, Genome, t_full_region, t_genseq, Family
from itertools import compress
from sqlalchemy import or_
from urllib.error import HTTPError
import sys
import plotly as py
from Bio import SeqIO, SeqRecord


def main():
    """Runs data processing scripts to turn raw data from (../raw) into
    cleaned data ready to be analyzed (saved in ../processed).
    """

    logger = logging.getLogger(__name__)
    logger.info("downloading genomes from NCBI")
    session = rfam_session()
    genomes_query = (
        session.query(Genome)
        .filter(or_(Genome.kingdom == "archaea", Genome.kingdom == "bacteria"))
        .filter(Genome.assembly_acc != "")
    )
    session.close()

    for genome in genomes_query:
        download_genome(genome)

    for genome in genomes_query:
        pickle_filename = "/home/jovyan/work/data/interim/igr_df_pickles/" + genome.assembly_acc + ".bz2"

        if not os.path.isfile(pickle_filename):
            igr_df = extract_igrs(genome)
            annotated_igr_df = annotate_igrs(genome, igr_df)

            annotated_igr_df.to_pickle(pickle_filename)


def download_genome(genome):
    """
    Given a genome object, download the appropriate Genbank flat file from NCBI
    :param genome: A sqlalchemy ORM object as defined in rfam_db
    """

    Entrez.email = os.environ.get("ENTREZ_EMAIL")
    Entrez.api_key = os.environ.get("ENTREZ_APIKEY")

    genbank_filename = "data/raw/download/{}_{}_genomic.gbff.gz".format(genome.assembly_acc, genome.assembly_name)

    if not os.path.isfile(genbank_filename):
        logger = logging.getLogger(__name__)
        logger.info("downloading genome assembly {} ({})".format(genome.assembly_acc, genome.scientific_name))
        # Search NCBI to get the record ID for a given Assembly ID
        search_term = "{} [Assembly Accession]".format(genome.assembly_acc)
        try:
            assembly_query_handle = Entrez.esearch(db="assembly", term=search_term, field="ASAC")

        except HTTPError as e:
            if e.code == 429:
                time.sleep(5)
                assembly_query_handle = Entrez.esearch(db="assembly", term=search_term, field="ASAC")
            else:
                raise

        assembly_query_result = Entrez.read(assembly_query_handle)
        assembly_query_handle.close()
        assembly_record_ids = assembly_query_result["IdList"]

        # Check to make sure a single record was returned, not multiple
        if len(assembly_record_ids) == 1:
            record_id = assembly_record_ids[0]
        else:
            # Get summaries of the duplicate files
            result = Entrez.read(
                Entrez.esummary(
                    id=",".join(assembly_query_result["IdList"]),
                    db="assembly",
                )
            )
            # Get a list of the assembly accessions
            accession_list = [
                document_summary["AssemblyAccession"]
                for document_summary in result["DocumentSummarySet"]["DocumentSummary"]
            ]
            # Check if the assembly accession matches search term (disregard GCA vs GCF)
            accession_match = [accession[3:] in search_term for accession in accession_list]
            # Extract the index of the matching accession
            index = list(compress(range(len(accession_list)), accession_match))

            if len(index) != 1:
                raise ValueError(
                    "{} records were returned by Entrez when searching for assembly {}".format(
                        len(assembly_record_ids), genome.assembly_acc
                    )
                )
            else:
                record_id = result["DocumentSummarySet"]["DocumentSummary"][index[0]].attributes["uid"]

        # Download the record summary for that genome assembly
        assembly_record_handle = Entrez.esummary(db="assembly", id=record_id)
        assembly_record_summary_set = Entrez.read(assembly_record_handle, validate=False)
        assembly_record_handle.close()
        assembly_record_summary = assembly_record_summary_set["DocumentSummarySet"]["DocumentSummary"][0]

        # Pull the FTP path from the assembly record summary.
        ftp_path = assembly_record_summary["FtpPath_GenBank"] + ""
        # Extract the part of the ftp_path that comes after the final '/'
        ftp_directory = ftp_path[ftp_path.rfind("/") + 1 :]
        # Build the full path to the genomic gbff.tz file
        genbank_ftp_path = "{}/{}_genomic.gbff.gz".format(ftp_path, ftp_directory)

        with open(genbank_filename, "wb") as genbank_file:
            request = urlopen(genbank_ftp_path)
            shutil.copyfileobj(request, genbank_file)

    return


def extract_igrs(genome, igr_length_cutoff=1):
    """For a genbank file, create a SeqRecord with all features annotated, identify all of the inter-genic regions,
    and store the corresponding data in a DataFrame"""

    genbank_filename = "data/raw/download/{}_{}_genomic.gbff.gz".format(genome.assembly_acc, genome.assembly_name)

    if not os.path.exists(genbank_filename):
        download_genome(genome)

    with gzip.open(genbank_filename, "rt") as genbank_file:
        igr_list = []

        # Create a Seqrecord for each "chromosome" in the file
        for seqrecord in SeqIO.parse(genbank_file, "gb"):
            cds_list = []

            # Loop over the genome file, get the CDS features on each of the strands
            for feature in seqrecord.features:
                if feature.type == "CDS":
                    coding_start = feature.location.start
                    coding_end = feature.location.end
                    cds_list.append((coding_start, coding_end))

            for i, position_pair in enumerate(cds_list):
                # Extract the sequence from current start position to previous end position
                igr_start = cds_list[i - 1][1]
                igr_end = position_pair[0]
                igr_sequence = seqrecord.seq[igr_start:igr_end]
                igr_len = len(igr_sequence)

                if igr_len >= igr_length_cutoff:  # Check if IGR is longer than cutoff.
                    igr_gc = GC(igr_sequence)
                    igr_list.append(
                        {
                            "accession": seqrecord.name,
                            "start": igr_start,
                            "end": igr_end,
                            "length": igr_len,
                            "gc": igr_gc,
                            "sequence": igr_sequence,
                        }
                    )
        igr_df = pd.DataFrame(igr_list, columns=["accession", "start", "end", "length", "gc", "sequence"])
        igr_df.index.name = "igr_index"

    return igr_df


def annotate_igrs(genome, igr_df):
    """
    Annotate the inter-genic regions listed in a dataframe with any available annotations from Rfam

    Parameters
    ----------
    genome: src.data.rfam_db.Genome
        The genome object for the organism who's IGR's are being analyzed
    igr_df: pandas.Dataframe
        The dataframe with the columns 'accession', 'start', 'end', 'length', 'gc'
    Returns
    -------
    annotated_igr_df: pandas.Dataframe
    """

    # Initialize connection to Rfam database
    session = rfam_session()

    # Get the list of "rfamseq_acc" numbers for a given organism
    rfamseq_acc_list = session.query(t_genseq.c.rfamseq_acc).filter(t_genseq.c.upid == genome.upid).distinct().all()

    # Create a list to store all the interval trees
    annotation_tree_dict = {}

    for rfamseq_acc in rfamseq_acc_list:
        # Pull rfamseq_acc out of the list
        rfamseq_acc = rfamseq_acc[0]

        rna_query = session.query(t_full_region).filter(t_full_region.c.rfamseq_acc == rfamseq_acc)
        rna_list = rna_query.all()

        # Make an interval tree for all of the RNA annotations to allow for rapid overlap search
        annotation_tree = IntervalTree()

        # Go though and add each RNA annotation to the interval tree
        for rna in rna_list:
            start = min(rna.seq_start, rna.seq_end)
            end = max(rna.seq_start, rna.seq_end)

            annotation_interval = Interval(start=start, end=end, chrom=rna.rfamseq_acc, value=rna)
            annotation_tree.insert_interval(annotation_interval)

        rfamseq_acc_stripped = rfamseq_acc.partition(".")[0]
        annotation_tree_dict[rfamseq_acc_stripped] = annotation_tree

    # Make an empty list of all the igrs with annotations
    annotated_igr_list = []
    for accession, accession_igr_df in igr_df.groupby("accession"):
        # Lookup the RNA annotation tree for the given accession
        try:
            annotation_tree = annotation_tree_dict[accession]
        except KeyError:
            print(
                "IGR dataframe key: {} not found. Available keys are: {}".format(accession, annotation_tree_dict.keys())
            )

        # For each IGR find all of the overlaps with annotated RNAs
        for igr in accession_igr_df.itertuples():
            overlap_list = annotation_tree.find(igr.start, igr.end)
            for overlap in overlap_list:
                # Add the IGR to the annotated_igr_list
                annotated_igr_list.append({"igr_index": igr[0], "rfam_acc": overlap.value.rfam_acc})

    # Convert annotated_igr_list into dataframe and merge on the rfam_acc
    annotated_igr_df = pd.merge(
        igr_df, pd.DataFrame(annotated_igr_list, columns=["igr_index", "rfam_acc"]), on="igr_index", how="left"
    )

    # Look up the information for all of the RNA families represented in this genome
    rna_family_query = (
        session.query(Family)
        .with_entities(Family.rfam_acc, Family.rfam_id, Family.description, Family.type)
        .filter(Family.rfam_acc.in_(annotated_igr_df["rfam_acc"].dropna().unique()))
    )
    rna_families_df = pd.read_sql(rna_family_query.statement, rna_family_query.session.bind)

    merged_igr_df = pd.merge(annotated_igr_df, rna_families_df, on="rfam_acc", how="left")

    combined_descriptions = (
        merged_igr_df.dropna()
        .groupby("igr_index")
        .agg(
            dict(
                rfam_acc=lambda x: ",".join(set(x)),
                rfam_id=lambda x: ",".join(set(x)),
                type=lambda x: ",".join(set(x)),
                description=lambda x: "<br>".join(set(x)),
            )
        )
    )
    merged_igr_df.drop_duplicates(["igr_index"], inplace=True)
    merged_igr_df.reset_index(inplace=True, drop=True)
    merged_igr_df.update(combined_descriptions)

    merged_igr_df["category"] = merged_igr_df.apply(lambda row: categorize_igr(row), axis=1)

    merged_igr_df["log_length"] = np.log(merged_igr_df["length"])
    session.close()
    return merged_igr_df


def categorize_igr(row):
    """
    Categorize an igr_df row based on a pre-defined set of labels
    :param row: pandas.Series
    :return: string
    """

    rfam_id_dict = {
        "rRNA": "rRNA",
        "tRNA": "tRNA",
        "intron": "Intron",
        "RNaseP": "RNase P",
        "6S": "6S RNA",
        "tmRNA": "tmRNA",
    }
    description_dict = {"riboswitch": "Riboswitch", "ribozyme": "Ribozyme"}
    type_dict = {"Gene; sRNA;": "sRNA"}
    category = ""

    if type(row.rfam_acc) == float:
        category = "No Known RNA"
    elif len(row.rfam_acc) > 7:
        category = "Multiple"
    else:
        for key, value in rfam_id_dict.items():
            if key in row.rfam_id:
                category = value
        for key, value in description_dict.items():
            if key in row.description:
                category = value
        if row.type == "Gene; sRNA;":
            category = "sRNA"
    if not category:
        category = "Miscellaneous"

    return category


# function to write .sto files:
def write_sto(directory, igr_name, igr_seq, suffix=""):
    if suffix != "":
        suffix = "_" + suffix
    output_filename = "{}/{}{}.sto".format(directory, igr_name.replace("/", "_"), suffix)
    with open(output_filename, "w") as output_file:
        output_file.write("# STOCKHOLM 1.0\n\n")
        output_file.write("{}       {}\n".format(igr_name, igr_seq))
        output_file.write("#=GC SS_cons{}{}\n//".format((" " * (len(igr_name) - 5)), ("." * len(igr_seq))))
        output_file.close()


if __name__ == "__main__":
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
