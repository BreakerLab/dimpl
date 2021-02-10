from Bio import Entrez
import pandas as pd
import numpy as np
import os
import urllib.request
import shutil
import gzip
import sys
from Bio import SeqIO
import configparser
import subprocess
from IPython.core.display import display, HTML, SVG
from urllib.error import HTTPError
import time
from Bio import AlignIO
import string
from itertools import compress
import genomeview
from src.shell.gff2bed import convert
import re

Entrez.email = os.environ.get("ENTREZ_EMAIL")
Entrez.api_key = os.environ.get("ENTREZ_APIKEY")

def get_nuccore_id(hit_accession):
    '''
    get the nuccore id of the hit by searching the nuccore database with the hit accession
    '''
    try:
        #search the nuccore database for info on the hit_accession
        nuccore_search_handle = Entrez.esearch(term=hit_accession, field='ACCN', db='nuccore')
        result = Entrez.read(nuccore_search_handle)
        nuccore_search_handle.close()
        nuccore_id = result['IdList'][0]

        return nuccore_id
    
    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return get_nuccore_id(hit_accession)
        else:
            raise

def fetch_deprecated_record_id(nuccore_id):
    '''
    returns record id for deprecated accessions
    '''
    try:
        fetch_record_handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="gb", retmode="xml")
        result = Entrez.read(fetch_record_handle)
        fetch_record_handle.close()
        summary = result[0]

        acc = result[0]['GBSeq_xrefs'][-1]['GBXref_id']

        search_term = "{} [Assembly Accession]".format(acc)
        assembly_query_handle = Entrez.esearch(db="assembly", term=search_term, field='ASAC')
        assembly_query_result = Entrez.read(assembly_query_handle)
        assembly_query_handle.close()
        assembly_record_ids = assembly_query_result["IdList"]
        
        # Check to make sure a single record was returned, not multiple. Then save assembly record ID.
        if len(assembly_record_ids) == 1:
            record_id = assembly_record_ids[0]
   
        else:
            # Get summaries of the duplicate files
            summary_handle= Entrez.esummary(id=','.join(assembly_query_result['IdList']), db='assembly')
            result = Entrez.read(summary_handle)
            summary_handle.close()
            # Get a list of the assembly accessions
            accession_list = [document_summary['AssemblyAccession'] for document_summary in result['DocumentSummarySet']['DocumentSummary']]
            # Check if the assembly accession matches search term (disregard GCA vs GCF)
            accession_match = [accession[3:] in search_term for accession in accession_list]
            # Extract the index of the matching accession
            index = list(compress(range(len(accession_list)), accession_match))
            
            if len(index) != 1:
                raise ValueError("{} records were returned by Entrez when searching for assembly {}".format(len(assembly_record_ids), genome.assembly_acc))
            else:
                record_id = result['DocumentSummarySet']['DocumentSummary'][index[0]].attributes['uid'] 

        return record_id
        
    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return get_nuccore_id(hit_accession)
        else:
            raise

def get_assembly_link(nuccore_id):
    '''
    returns record ID from assembly link.
    '''
    try:
        #get link to genome assembly information
        assembly_link_handle =Entrez.elink(id = nuccore_id, dbfrom = 'nuccore', linkfrom = 'nuccore_assembly' , db = 'assembly')
        assembly_query_result = Entrez.read(assembly_link_handle)
        assembly_link_handle.close()
        #print(assembly_query_result)
        
        try:
            #save assembly link ID for assembly information
            assembly_record_ids = assembly_query_result[0]['LinkSetDb'][0]
            
            # Check to make sure a single record was returned, not multiple. Then save assembly record ID.
            if len(assembly_record_ids['Link']) == 1:
                record_id = assembly_record_ids['Link'][0]['Id']
                
            else:
                record_id = fetch_deprecated_record_id(nuccore_id)
            
        except IndexError as ierr:
            record_id = fetch_deprecated_record_id(nuccore_id)

        return record_id
    
    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return get_assembly_link(nuccore_id)
        else:
            raise

def get_assembly_document(record_id):
    '''
    get the assembly accession of the hit from summary information of assembly record
    if nuccoreid is deprecated - and record_id is 0 - then return 0
    '''
    try:
        #get information on assembly for the specific assembly accession
        assembly_record_summary_handle = Entrez.esummary(db="assembly", id=record_id)
        result = Entrez.read(assembly_record_summary_handle, validate = False)
        assembly_record_summary_handle.close()

        #extract assembly record summary for ftp path later
        assembly_record_summary = result['DocumentSummarySet']['DocumentSummary'][0]

        return assembly_record_summary

    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return get_assembly_document(record_id)
        else:
            raise

def find_seq_in_alignment(id, alignment):
    '''returns alignment sequence'''
    for seqrecord in alignment:
        if seqrecord.id == id:
            s = str(seqrecord.seq)
            s = s.replace('.', '')
            s = s.replace('-', '')
            return s

def get_taxonomy(nuccore_id):
    
    try:
        record_id = get_assembly_link(nuccore_id)
        tax_id = get_assembly_document(record_id)['Taxid']
        taxonomy_info_handle= Entrez.efetch(db = 'taxonomy', id = tax_id, retmode = 'xml')
        result = Entrez.read(taxonomy_info_handle)
        taxonomy_info_handle.close()
        return result[0]['Lineage']
    
    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return get_taxonomy(nuccore_id)
        else:
            raise

def download_gff_fna(hit_accession):
    '''
    downloads the gff and fasta files for the hit accession
    returns filename
    '''
    
    #get link to assembly record 
    nuccore_id = get_nuccore_id(hit_accession)
    record_id = get_assembly_link(nuccore_id)

    #get document information of assembly record
    assembly_record_summary = get_assembly_document(record_id)

    #get assembly information
    assembly_accn = assembly_record_summary['AssemblyAccession']
    assembly_name = assembly_record_summary['AssemblyName']
    
    # Pull the FTP path from the assembly record summary.
    ftp_path = assembly_record_summary['FtpPath_RefSeq'] + ''
    base_filename = ftp_path[ftp_path.rfind("/") + 1:]

    output_folder= "data/raw/download"
    #if the output folder doesn't exist, create it, in the directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    #create filenames for the gff and fasta files
    refseq_gff_zip_filename = "data/raw/download/{}_genomic.gff.gz".format(base_filename)
    refseq_fasta_zip_filename = "data/raw/download/{}_genomic.fna.gz".format(base_filename)
    
    if not os.path.isfile(refseq_gff_zip_filename):

        # Build the full path to the genomic gff/fna files
        refseq_gff_zip_ftp_path = "{}/{}_genomic.gff.gz".format(ftp_path, base_filename)
        refseq_fasta_zip_ftp_path = "{}/{}_genomic.fna.gz".format(ftp_path, base_filename)

        #create gff/fna files
        with open(refseq_gff_zip_filename, 'wb') as refseq_gff_zip_file:
            request_gff = urllib.request.urlopen(refseq_gff_zip_ftp_path)

            #copy the content of source file to destination file.
            shutil.copyfileobj(request_gff, refseq_gff_zip_file)

        with open(refseq_fasta_zip_filename, 'wb') as refseq_fasta_zip_file:
            request_fasta = urllib.request.urlopen(refseq_fasta_zip_ftp_path)
            shutil.copyfileobj(request_fasta, refseq_fasta_zip_file)


        # Build the full path to the UNZIPPED genomic gff/fna files
        refseq_gff_ftp_path_unzip = "{}/{}_genomic.gff".format(ftp_path, base_filename)
        refseq_fasta_ftp_path_unzip = "{}/{}_genomic.fasta".format(ftp_path, base_filename)

        #unzip gff.gz file and save as .gff file
        input_gff = gzip.GzipFile(refseq_gff_zip_filename, 'rb')
        s_gff = input_gff.read()
        input_gff.close()

        output_gff = open("data/raw/download/{}_genomic.gff".format(base_filename), 'wb')
        output_gff.write(s_gff)
        output_gff.close()

        #unzip fna.gz file and save as .fna file
        input_fna = gzip.GzipFile(refseq_fasta_zip_filename, 'rb')
        s_fna = input_fna.read()
        input_fna.close()

        output_fna = open("data/raw/download/{}_genomic.fna".format(base_filename), 'wb')
        output_fna.write(s_fna)
        output_fna.close()
        
    return base_filename

def build_context_image(hit_row, alignment, upstream_range = 4000, downstream_range = 4000):
    
    #extract hit accession of target from hit_row
    hit_accession= hit_row.target_name
    
    #extract start nt of target, stop nt of target and strand from hit_row
    
    start= hit_row.seq_from
    stop= hit_row.seq_to
    seq_length = stop-start+1
    strand = hit_row.strand
    
    
    #if strand is + then flip is true, if strand is - then flip is false. This is needed in the genome browser url.
    flip_val = 'false'
    if(strand == "-"):
        flip_val = 'true'

    #get the assembly accession of the hit and download the fasta and gff files
    if hit_row.assembly_accession=='nan':
        base_filename = download_gff_fna(hit_accession)
    else:
        base_filename = hit_row.assembly_accession
    
    fna_file = "data/raw/download/"+base_filename+"_genomic.fna"
    gff_file = "data/raw/download/"+base_filename+"_genomic.gff"

    #extract sequence of the part of the target that matches the query
    target_sequence = find_seq_in_alignment(hit_row.target_coords, alignment)
    target_sequence = target_sequence.replace('T','U')
    #extract E-value
    e_value = hit_row.e_value
    
    #extract %gc
    percent_gc = hit_row.gc
    
    #extract score
    score = hit_row.score   
    #extract target-name
    target_name = hit_row.target_coords
    
    #extract taxonomy
    if  hit_row.lineage=='nan':
        lineage = get_taxonomy(hit_accession)
    else:
        lineage = hit_row.lineage
    
    #set range
    down_limit = start - downstream_range
    up_limit = stop + upstream_range
    
    #print statements for variables to be shown in image

    print("Match #{}".format(int(hit_row.name)+1))
    print("E-value: "+ str(e_value))
    print("%GC: "+ str(percent_gc))
    print("Score: "+ str(score))
    print("Genome Assembly: " +str(base_filename))
    print("Target: "+ target_name)
    print("Lineage: " + lineage)
    print("Matched Sequence: "+target_sequence)

   
    #clickable link
    genome_browser_url ='https://www.ncbi.nlm.nih.gov/projects/sviewer/?id={}&v={}:{}&c=FF6600&theme=Details&flip={}&select=null&content=3&color=0&label=1&geneModel=0&decor=0&layout=0&spacing=0&alncolor=on&m={},{}&mn=5,3'.format(hit_accession, down_limit, up_limit, flip_val, start, stop)   
    try:
        display(HTML('<a href="{}")>genome browser</a>'.format(genome_browser_url)))
    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return display(HTML('<a href="{}")>genome browser</a>'.format(genome_browser_url)))
        else:
            raise

    gff_file_zip="/home/jovyan/work/data/raw/download/"+ base_filename+ "_genomic.gff.gz"
    bed_file="/home/jovyan/work/data/raw/features/"+ base_filename+ "_genomic.bed"

    
    convert(gff_file_zip,bed_file, desc_only=True)

    def prerender(renderer, element):
        # Prerenderers get run before the track is rendered.
        # Draw non-feature graphic elements.

        # IGR is in forward orientation
        if start < stop:
            x1 = element.scale.topixels(start) # converting genomic coordinates to screen coordinates
            x2 = element.scale.topixels(stop)
            # Draw vertical hit bar
            yield from renderer.rect(x1, 0, x2-x1, element.height-14, fill="lightblue", stroke="none")
            # Add "HIT" text
            yield from renderer.text(x1+(x2-x1)/2, element.height-2, "HIT", size=12, anchor="middle")
            # Draw hit direction arrow (forward/positive stand)
            yield from renderer.block_arrow(x1+(x2-x1)/2+16, element.height-11, 0, 9, arrow_width=9,
                direction="right", fill="black", stroke="none")

        # IGR is in reverse orientation
        if start > stop:
            x1 = element.scale.topixels(start) # converting genomic coordinates to screen coordinates
            x2 = element.scale.topixels(stop)
            # Draw vertical hit bar
            yield from renderer.rect(x2, 0, x1-x2, element.height-14, fill="lightblue", stroke="none")
            # Add "HIT" text
            yield from renderer.text(x1+(x2-x1)/2, element.height-2, "HIT", size=12, anchor="middle")
            # Draw hit direction arrow (reverse/negative stand)
            yield from renderer.block_arrow(x2+(x1-x2)/2-16, element.height-11, 0, 9, arrow_width=9,
                direction="left", fill="black", stroke="none")

    doc=genomeview.visualize_data({"":bed_file},hit_accession,down_limit,up_limit)
    cur_track = genomeview.get_one_track(doc, "")
    cur_track.prerenderers = [prerender]

    display(doc)

    return base_filename, lineage


def get_all_images(results_csv_filename, alignment):

    results_df = pd.read_csv(results_csv_filename, index_col=0)
    results_df['lineage'] = results_df['lineage'].astype(str)
    results_df['assembly_accession'] = results_df['assembly_accession'].astype(str)
    updated_results_df = results_df.copy(deep=True)
    for index, row in results_df.iterrows():
        assembly_accession, lineage = build_context_image(row, alignment, upstream_range = 4000, downstream_range = 4000)
        updated_results_df.loc[index, 'lineage'] = lineage
        updated_results_df.loc[index, 'assembly_accession'] = assembly_accession
    updated_results_df.to_csv(results_csv_filename)
        
def build_target_coords(target_name, seq_from, seq_to):
    return "{}/{}-{}".format(target_name, seq_from, seq_to)

def run_rscape(outdir, sto_filename, fold=True, output=True):
    truncated_filename = sto_filename[sto_filename.rfind('/')+1:sto_filename.rfind('.sto')]
    
    if fold:
        arguments = ['R-scape', '--fold', '--outdir', outdir, sto_filename]
    else:
        arguments = ['R-scape', '--r2rall', '--outdir', outdir, sto_filename]
    
    result = subprocess.run(arguments, capture_output=True)
    if output:
        print(result.stdout.decode())
    
    # List of the suffixes of excess files to delete
    deleted_file_suffix = ['cov', 'dplot.ps','dplot.svg', 'power',  'sorted.cov', 'surv', 'surv.ps', 'surv.svg', 'R2R.sto', 'R2R.sto.pdf' ]
    
    for suffix in deleted_file_suffix:
        file_to_delete = "{}/{}_1.{}".format(outdir, truncated_filename, suffix)
        if os.path.exists(file_to_delete):
            os.remove(file_to_delete)
    
    svg_filename = "{}/{}_1.R2R.sto.svg".format(outdir, truncated_filename)
    
    if fold:
        deleted_fold_suffix = ['dplot.ps', 'dplot.svg', 'power', 'R2R.sto', 'R2R.sto.pdf', 'cov']
        for suffix in deleted_fold_suffix:
            file_to_delete = "{}/{}_1.fold.{}".format(outdir, truncated_filename, suffix)
            if os.path.exists(file_to_delete):
                os.remove(file_to_delete)
        os.remove("{}/{}_1.sorted.fold.cov".format(outdir, truncated_filename, suffix))
        svg_filename = "{}/{}_1.fold.R2R.sto.svg".format(outdir, truncated_filename)

    
    display(SVG(filename=svg_filename))
    
def tar_subdir_members(tar, import_tar_name):
    tar_foldername = re.sub(r'(\.done)*\.tar\.gz','',import_tar_name.split('/')[-1]) + '/'
    tar_foldername_length = len(tar_foldername)
    for member in tar.getmembers():
        if member.path.startswith(tar_foldername):
            member.path = member.path[tar_foldername_length:]
            yield member
