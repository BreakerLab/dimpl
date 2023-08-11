#!/usr/bin/env python3
#
#  GFF3 to BED6 format file converter.
#
#  Glenn Gaffield - Nov 2019
#
#  COMMAND LINE USAGE:  python gff2bed.py [--make-genome-file | -g] [--desc-only | -d] gff-filename [ bed-filename ]
#
#  "bed-filename" is optional.  If not specified, the bed file is written to the same directory
#  as the gff-filename argument.  To write to STDOUT use "-".
#
#  OPTIONS:
#  --make-genome-file | -g      Creates a .genome file in the same directory as the .bed file,
#                               that contains accession number and length of sequence.
#                               This file is typically used by BEDTools.
#  --desc-only | -d             Don't display IDs when something other than an ID is available.
#                               In other words, only show ID as a last resort.
#
#  Python module usage: gff2bed.convert(gff_filename, bed_filename, [make_genome_file_flag], [desc_only_flag])
#
#  BED field#	1	2	3	4	5	6
#  GFF field#	1	4	5	9**	6	7
#    ** = where the data comes from depends on feature type

LABEL_SIZE = 80

# Ignore "Broken Pipe" error when writing to stdout and piping to something the terminates output prematurely (ex. more, less, head)
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)
import os

#
# FUNCTIONS


# Open a file, and handle gzip'd files using gzip module
def openfile(filename, mode="r"):
    import sys
    import gzip

    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    elif filename.endswith("genome"):
        return open(filename, mode)
    elif not filename:
        # If not generating a genome file, dump those writes to /dev/null
        return open(os.devnull, mode)
    elif filename == "-" and mode == "w":
        return sys.stdout
    else:
        return open(filename, mode)


def convert(gff_file, bed_file, make_genome_file=False, desc_only=False):
    import gzip
    import re

    attributes = {}
    attribs_saved = {}
    previous_ids = []
    previous_data = {}
    previous_accno = ""

    if not bed_file:
        bed_file = gff_file.replace(".gff", ".bed")
        bed_file = gff_file.replace(".gff.gz", ".bed")

    genome_file = ""
    if make_genome_file and bed_file != "-":
        genome_file = bed_file.replace(".bed", ".genome")

    with openfile(gff_file, "r") as gff:
        with openfile(bed_file, "w") as bed:
            with openfile(genome_file, "w") as genome:
                EOF = False
                first_line = True

                # Read gff file line-by-line until end-of-file is detected
                while not EOF:
                    line = gff.readline().decode("utf-8")
                    if not line:
                        EOF = True
                        gff_id = ""
                    else:
                        # Ignore comment lines
                        if line[0] == "#":
                            continue

                        # Parse Line into fields
                        (
                            accno,
                            source,
                            gff_type,
                            start,
                            end,
                            score,
                            strand,
                            phase,
                            gff_attributes,
                        ) = line.rstrip().split("\t")

                        # "region" record is used to generate the genome file
                        if gff_type == "region" and start == "1" and accno != previous_accno:
                            size = end
                            # Write to genome file, accession number and genome size
                            genome.write(accno + "\t" + size + "\n")
                            previous_accno = accno

                        # Ignore "region" and "repeat" features
                        if gff_type in ["region", "direct_repeat", "inverted_repeat", "repeat_region", "tandem_repeat"]:
                            continue

                        # Parse GFF field 9 attributes and collect them together as we read lines from the file until the exon changes
                        attrib_string = gff_attributes.split(";")
                        attributes = {}
                        for attribute in attrib_string:
                            (tag, value) = attribute.split("=")
                            # print(tag+" is |"+value+"|")
                            attributes[tag] = value

                        # Group related records
                        gff_id = attributes.get("Parent", attributes.get("ID"))

                    # Create 1 bed file line per unique exon in the gff file.
                    # Multipl gff file lines are grouped together by IDs found in the "ID" and "parent" fields.
                    if not first_line and (gff_id not in previous_ids):
                        skip = False
                        p = previous_data

                        # Create a meaningful label for the feature that (ideally) consists of three parts:
                        #   Feature Type + Gene Name + Some ID

                        # Start by assigning data fields to the three parts of the label (which covers about 90% of the data).
                        # Then we'll override these values for the exceptions to these assignments.
                        feature_type = attribs_saved.get("gene_biotype", p["type"])
                        feature_gene = attribs_saved.get("gene", "")
                        feature_id = attribs_saved.get("Name", attribs_saved.get("Dbxref", ""))
                        feature_desc = attribs_saved.get(
                            "product", attribs_saved.get("Note", attribs_saved.get("function", ""))
                        )

                        # "mobile_genetic_element" feature type sometimes has info in the "mobile_element_type" field.
                        # Otherwise, just label it as "mobile_genetic_element".
                        if feature_type == "mobile_genetic_element":
                            feature_desc = attribs_saved["mobile_element_type"]
                            if ":" in feature_desc:
                                (feature_desc, feature_id) = feature_desc.split(":")
                            if feature_desc == "other":
                                feature_desc = feature_id
                                feature_id = ""

                        # The pseudogene feature type has an ID on the end of the "inference" field
                        if feature_type in ["pseudogene"]:
                            feature_infer_str = attribs_saved.get("inference", "")
                            feature_infer_arr = feature_infer_str.split(":")

                            if feature_infer_str and feature_infer_arr[-2] in ["RefSeq", "HMM", "SwissProt"]:
                                feature_id = feature_infer_arr[-1]
                                # Mark ID as "approximate" when it's commented as being "similar" to the annotated protein ID.
                                if feature_infer_str.find("similar") >= 0:
                                    feature_id = "≈" + feature_id

                        # "sequence_feature" and "binding_site" feature types have useful info in the "Note" or "function" fields.
                        if feature_type in ["sequence_feature", "binding_site"]:
                            #                    feature_desc = attribs_saved.get('Note',attribs_saved.get('function',''))

                            # Skip features like: "23S ribosomal RNA rRNA prediction is too short"
                            if feature_desc.endswith("too short"):
                                skip = True
                            feature_type = ""

                        # If gene name is also in the id, replace id with the contents of the "Dbxref" field.
                        # If gene name is not blank, append ":".
                        if feature_gene:
                            if feature_id == feature_gene:
                                feature_id = attribs_saved.get("Dbxref", "")
                                # Cut off characters after (and including) a comma.
                                n = feature_id.find(",")
                                feature_id = feature_id[0:n] if n >= 0 else feature_id

                        # If feature_type is "protein_coding", use the "product" field desciption in place of feature_type.
                        if feature_type == "protein_coding":
                            feature_desc = attribs_saved.get("product", "")
                            if feature_desc.find("protein") < 0:
                                feature_desc += " protein"
                            feature_type = ""

                        # If feature_type is "riboswitch", use the "Note" field desciption in place of feature_type.
                        if feature_type == "riboswitch":
                            feature_desc = attribs_saved.get("Note", "")
                            if feature_desc.find("riboswitch") >= 0:
                                feature_type = ""

                        # Final feature label is made up of three parts (any one or two could be missing)
                        # In description only mode, do not append ID unless the rest of the label is blank.
                        div1 = ""
                        div2 = ""
                        if feature_type and feature_gene or feature_type and feature_desc:
                            div1 = ":"
                        if feature_gene and feature_desc:
                            div2 = ":"

                        # Cleanup and trim description, if one exists.
                        if feature_desc:
                            # Remove leading numbers
                            feature_desc = re.sub(r"^[0-9]+ ", "", feature_desc)
                            # Remove Pfam & Prosite IDs, if option to show only descriptions was selected.
                            if desc_only:
                                feature_desc = re.sub(r"^(Pfam.*PF[0-9]+|PS[0-9]+) ", "", feature_desc)

                            # Special case where Note contains fasta comment format, field sep = |
                            # The last field contains useful label info.
                            n = feature_desc.rfind("|")
                            if n >= 0:
                                feature_desc = feature_desc[n + 1 :].strip()
                            # Trim description to approximately 20 chars, trim at word boundary only.
                            # Final length can be greater than 20 chars if adding a few more chars would create a full word.
                            if len(feature_desc) > LABEL_SIZE:
                                i = feature_desc.rfind(" ", 0, LABEL_SIZE)
                                if i < 0:
                                    i = LABEL_SIZE - 1
                                j = feature_desc.find(" ", LABEL_SIZE - 1)
                                if j < 0:
                                    j = len(feature_desc)
                                if LABEL_SIZE - 1 - i < j - LABEL_SIZE - 1:
                                    k = i
                                else:
                                    k = j
                                feature_desc = feature_desc[0:k]

                        label = feature_type + div1 + feature_gene + div2 + feature_desc

                        if (desc_only and not label) or not desc_only:
                            div3 = ""
                            if label and feature_id:
                                div3 = ":"
                            label += div3 + feature_id

                        # Fix encoded characters
                        label = label.replace("%2C", ",").replace("%3B", ";").replace("%25", "%")

                        # BED format uses a zero-based start position, so subtract 1 from GFF start position.
                        zero_based_start = str(int(p["start"]) - 1)

                        bed_line = "\t".join([p["accno"], zero_based_start, p["end"], label, p["score"], p["strand"]])

                        if not skip:
                            # print("BED:\t"+bed_line)
                            bed.write(bed_line + "\n")

                        attribs_saved = {}
                        previous_ids = []

                    previous_ids.append(attributes.get("ID"))

                    previous_data = {
                        "accno": accno,
                        "start": start,
                        "end": end,
                        "type": gff_type,
                        "score": score,
                        "strand": strand,
                        "phase": phase,
                    }
                    first_line = False

                    # Merge the current record's attributes with those from the previous record(s)
                    attribs_saved.update(attributes)
    return


if __name__ == "__main__":
    import sys
    import os

    arg_index = 1
    make_genome_file = False
    desc_only = False

    for arg in sys.argv:
        if arg.startswith("--make-genome-file") or arg.startswith("-g"):
            arg_index += 1
            make_genome_file = True

        if arg.startswith("--desc-only") or arg.startswith("-d"):
            arg_index += 1
            desc_only = True

    args = len(sys.argv) - arg_index
    if args == 0:
        print("Error: GFF filename argument is required.")
        sys.exit(1)

    if args >= 1:
        gff_file = sys.argv[arg_index]
        bed_file = ""

    if args == 2:
        bed_file = sys.argv[arg_index + 1]

    convert(gff_file, bed_file, make_genome_file, desc_only)
