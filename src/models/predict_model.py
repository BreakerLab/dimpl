import re
import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqUtils import GC
from src.data.make_dataset import write_sto


def extract_parameters(selection_file_name):
    regex = r"([^/]+)\_selection_(-?\d+\.\d)_(-?\d+\.\d)_(-?\d+\.\d).fasta.xml"
    result = re.search(regex, selection_file_name)
    if result:
        assembly_acc = result.groups()[0]
        class_weight_mod = float(result.groups()[1])
        c_exp = float(result.groups()[2])
        gamma_exp = float(result.groups()[3])

        return assembly_acc, class_weight_mod, c_exp, gamma_exp
    else:
        raise ValueError("Could not extract selection parameters from filename {}".format(selection_file_name))


def process_blast(
    fasta_dir,
    xml_dir,
    svm_clf,
    step_dirname="infernal_step1",
    suffix="igr",
    poor_blast_score=40,
    poor_score_weight=0.5,
    hypothetical_weight=0.5,
    orf_score_cutoff=4,
    score_increment=2,
    svm_decision_cutoff=0.0,
):
    parent_dir = "/".join(fasta_dir.split("/")[:-1])
    step_dir = "{}/{}".format(parent_dir, step_dirname)

    if not os.path.exists(step_dir):
        os.makedirs(step_dir)

    with open(step_dir + "/blast_log.txt", "w") as log_file:
        # Get list of all fasta files
        fasta_filelist = glob.glob(fasta_dir + "/*.fasta")

        for fasta_file in fasta_filelist:
            fasta_record = [record for record in SeqIO.parse(fasta_file, "fasta")]

            igr_address = ".".join(fasta_file.split("/")[-1].split(".")[:-1])
            blast_filename = xml_dir + "/" + igr_address + ".xml"
            with open(blast_filename, "r") as blast_file:
                blast_records = NCBIXML.parse(blast_file)

                # For each blast_record, make a nucleotide_orf_score array of zeros to represent the length
                # and position of nucleotide in each of the query sequences:
                for blast_record in blast_records:
                    nucleotide_orf_score = [0] * blast_record.query_length

                    # Define query_record and query_record.seq as the name and nucleotide query
                    # sequence, respectively. Also, define the genomic start and end positions.

                    query_record = fasta_record[0]
                    query_record.seq = query_record.seq[0 : len(query_record.seq) - 1]
                    parse_igr_name = re.search(r"\S+_(\d+)-(\d+)", igr_address)
                    genomic_start = int(parse_igr_name.groups()[0])
                    genomic_stop = int(parse_igr_name.groups()[1])

                    log_file.write("IGR {} is divided into:\n".format(igr_address))

                    # Now, for each alignment in the BLASTx search, perform this for each "hit":
                    for alignment in blast_record.alignments:
                        # For each hit (hsp), if the hit score is less than 40, change the
                        # nucleotide_flag to 1, otherwise change nucleotide_flag to 2.
                        if "hypothetical protein" in alignment.hit_def.lower():
                            alignment_score_increment = score_increment * hypothetical_weight
                        else:
                            alignment_score_increment = score_increment

                        for hsp in alignment.hsps:
                            # If BLASTx score for alignment is <40 or a black bar, then increment by 1,
                            # else 2.
                            if hsp.score < poor_blast_score:
                                hsp_score_increment = alignment_score_increment * poor_score_weight
                            else:
                                hsp_score_increment = alignment_score_increment

                            # For each nucleotide in the range from the start to the end of the
                            # query sequence, keep track of the nucleotide_flag.
                            for i in range(hsp.query_start - 1, hsp.query_end):
                                nucleotide_orf_score[i] = nucleotide_orf_score[i] + hsp_score_increment

                    # Initialize the program at the beginning of query sequence.
                    current_location = -1
                    start = 0

                    # Assume we start in an intergenic region.
                    in_igr_flag = 1
                    # If multiple hits overlap, then we are not starting in an intergenic region.
                    if nucleotide_orf_score[0] >= orf_score_cutoff:
                        in_igr_flag = 0

                        # Loop through each nucleotide until the end of query.
                    while not len(nucleotide_orf_score) == 0:
                        nucleotide_flag = nucleotide_orf_score.pop(
                            0
                        )  # Pop removes one nucleotide from the nucleotide_orf_score for each iteration.
                        current_location = current_location + 1

                        # If in protein coding region, in_igr_flag = 0.
                        if not in_igr_flag:
                            # If protein hit is confident continue to pop out nucleotides.
                            if nucleotide_flag >= orf_score_cutoff:
                                continue
                            # If protein hit is not confident, flag as a potential IGR (in_igr_flag = 1) and continue to else statement.
                            else:
                                start = current_location
                                in_igr_flag = 1
                        # If currently in an igr, in_igr_flag = 1.
                        else:
                            # If you reach the end of an intergenic region, create a seq_record for the IGR you just finished.
                            if nucleotide_flag >= orf_score_cutoff or len(nucleotide_orf_score) == 0:
                                stop = current_location - 1
                                # If you reach the end of the query sequence,
                                if len(nucleotide_orf_score) == 0:
                                    stop = stop + 1
                                sequence_start = genomic_start + start
                                sequence_stop = genomic_start + stop
                                accession = igr_address.split("_")[0]

                                sub_record = query_record[start : stop + 1]
                                sub_record.id = "{}_{}-{}".format(accession, sequence_start, sequence_stop)

                                # If the seq_record meets the IGR length_cutoff and GC_cutoff requirements, record the new genomic start/end sites and new corresponding nucleotide sequence.
                                # Record action as "kept" in log file.
                                if len(sub_record.seq) > 0:
                                    sub_record_series = pd.DataFrame(
                                        [[GC(sub_record.seq), np.log(len(sub_record))]], columns=["gc", "log_length"]
                                    )
                                    prediction = svm_clf.decision_function(sub_record_series)

                                    if prediction > svm_decision_cutoff:
                                        sto_dir = "{}/{}_{}".format(step_dir, sub_record.id, suffix)
                                        if not os.path.exists(sto_dir):
                                            os.makedirs(sto_dir)
                                        write_sto(sto_dir, sub_record.id, sub_record.seq, suffix=suffix)

                                        log_string = "\tKept:   {} length: {} GC: {:.2f}%\n".format(
                                            sub_record.id, len(sub_record), GC(sub_record.seq)
                                        )
                                        log_file.write(log_string)

                                    # If the seq_record does not meet the minimum cut_off requirements, toss and record action as "tossed" in log file.
                                    else:
                                        log_string = "\tTossed: {} length: {} GC: {:.2f}%\n".format(
                                            sub_record.id, len(sub_record), GC(sub_record.seq)
                                        )
                                        log_file.write(log_string)
                                in_igr_flag = 0
                            # If the IGR continues
                            else:
                                continue

    log_file.close()
