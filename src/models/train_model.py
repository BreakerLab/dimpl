import pandas as pd
from sklearn import svm
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def genome_svm_selection(annotated_df, class_weight_mod=0.5, c_exp=3, gamma_exp=-2, output_fasta=False, output_prefix=""):


    y = annotated_df['category'] != 'No Known RNA'
    total_igrs = len(y)
    total_knowns = y.sum()
    total_unknowns = total_igrs - total_knowns

    class_weight = {False: total_knowns / total_igrs, True: (total_unknowns / total_igrs * class_weight_mod)}
    svm_clf = svm.SVC(C=10 ** c_exp, class_weight=class_weight, gamma=10 ** gamma_exp)
    svm_clf.fit(annotated_df.loc[:, ["gc", "log_length"]], y)


    selection = pd.Series(svm_clf.predict(annotated_df.loc[:, ["gc", "log_length"]]))
    selected_unknowns = selection & annotated_df['rfam_acc'].isnull()



    selected_igr_list = [SeqRecord.SeqRecord(row.sequence, id=("{}/{}-{}".format(row.accession, row.start + 1, row.end)))
        for row in annotated_df.loc[selected_unknowns, ["accession", "start", "end", "sequence"]].itertuples()]

    if output_fasta:
        if output_prefix == "":
            raise ValueError("Output prefix is a required argument when outputting a fasta")
        outputfilename = "data/interim/{}_selection_{}_{}_{}.fasta".format(output_prefix, class_weight_mod, c_exp,
                                                                           gamma_exp)
        SeqIO.write(selected_igr_list, outputfilename, "fasta")