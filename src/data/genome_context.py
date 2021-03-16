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
from IPython.core.display import display, HTML, SVG, Markdown
from urllib.error import HTTPError, URLError
import time
from Bio import AlignIO
import string
from itertools import compress
import genomeview
from src.shell.gff2bed import convert
import re
import base64
from dotenv import load_dotenv, find_dotenv

# Load OS environment variables from ".env" file.
dotenv_path = find_dotenv()
load_dotenv(dotenv_path,override=True)

Entrez.email = os.environ.get("ENTREZ_EMAIL")
Entrez.api_key = os.environ.get("ENTREZ_APIKEY")
Entrez.sleep_between_tries = 1

# Breaker Lab Environmental Dataset - requires Yale VPN access
BL_url_top = "http://bl.biology.yale.edu/lab/scripts"
BL_api = BL_url_top + "/" + "dimpl_api.pl"
BL_user = "Dimpl_api"
BL_api_key = os.environ.get("BL_APIKEY")

# A tidy (no traceback) way of stopping after a trapped error has been handled.
# Use "raise StopException".
class StopExecution(Exception):
    def _render_traceback_(self):
        pass

def init_environmentals_api():
    '''
    setup digest authentication for Breaker Lab server access
    '''
    global BLaccess
    
    if not BL_api_key:
        BLaccess = False
        return
    
    # Create an authentication manager
    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    
    # Parse API key
    try:
        apikey = base64.standard_b64decode(BL_api_key[::-1]+'==').rstrip().decode('utf-8')
        digestauth = base64.standard_b64decode(apikey.split()[1][::-1]+'==').rstrip().decode('utf-8')
        keyuser = base64.standard_b64decode(apikey.split()[2]+'==').rstrip().decode('utf-8')
        #print("Authenticated Breaker Lab User:",keyuser)
    except UnicodeDecodeError:
        print("BL API key error.  Check for typo or acquire a new key.")
        print("    key=",BL_api_key)
        raise StopExecution
    
    # Add digest credentials to the password manager
    password_mgr.add_password(None, BL_url_top, BL_user, digestauth)

    handler = urllib.request.HTTPDigestAuthHandler(password_mgr)

    # create "opener" (OpenerDirector instance)
    global opener
    opener = urllib.request.build_opener(handler)

    # Check if Breaker Lab server is accessible
    try:
        opener.open(BL_api)
    except HTTPError as e:
        if e.code == 443:
            #print("BL server found.")
            BLaccess = True

            # Install the opener.
            # Now all calls to urllib.request.urlopen use our opener.
            urllib.request.install_opener(opener)
        else:
            print("ERROR: BL server access problem:",e)
            raise StopExecution
    except URLError as e:
        print("BL server not found.",e)
        BLaccess = False
        raise StopExecution
        
def get_nuccore_id(hit_accession):
    '''
    get the nuccore id of the hit by searching the nuccore database with the hit accession
    '''
    try:
        #search the nuccore database for info on the hit_accession
        nuccore_search_handle = Entrez.esearch(term=hit_accession, field='ACCN', db='nuccore')
        result = Entrez.read(nuccore_search_handle)
        nuccore_search_handle.close()
    
        # Check if an empty set is returned
        if(result['IdList']):
            nuccore_id = result['IdList'][0]
            return nuccore_id
        else:
            # print(result['WarningList']['OutputMessage'])
            return "NOT FOUND"
    
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
            return fetch_deprecated_record_id(nuccore_id)
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
    except IndexError:
        time.sleep(0.5)
        return get_assembly_document(record_id)

def find_seq_in_alignment(id, alignment):
    '''returns alignment sequence'''
    for seqrecord in alignment:
        # Remove nn| prefix before doing comparison
        if re.sub('^[0-9]+\|','',seqrecord.id) == id:
            s = str(seqrecord.seq)
            s = s.replace('.', '')
            s = s.replace('-', '')
            return s

def get_taxonomy(tax_id):
    '''
    returns lineage given a tax_id
    '''
    try:
        taxonomy_info_handle= Entrez.efetch(db = 'taxonomy', id = tax_id, retmode = 'xml')
        result = Entrez.read(taxonomy_info_handle)
        taxonomy_info_handle.close()
        return result[0]['Lineage']
    
    except HTTPError as e:
        if(e.code == 429):
            time.sleep(0.5)
            return get_taxonomy(tax_id)
        else:
            raise

def download_gff(hit_row):
    '''
    downloads a gff file for the hit accession
    returns filename
    '''
    hit_accession = hit_row.target_name
    
    output_folder= "data/raw/download"
    #if the output folder doesn't exist, create it, in the directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    #get link to assembly record 
    nuccore_id = get_nuccore_id(hit_accession)
    found_on_refseq = (nuccore_id != 'NOT FOUND')

    if found_on_refseq:
        record_id = get_assembly_link(nuccore_id)
    
        #get document information of assembly record
        assembly_record_summary = get_assembly_document(record_id)

        #get assembly information
        hit_row.tax_id = assembly_record_summary['Taxid']
        
        # Pull the FTP path from the assembly record summary.
        ftp_path = assembly_record_summary['FtpPath_RefSeq'] + ''
        base_filename = ftp_path[ftp_path.rfind("/") + 1:]
    else:
        base_filename = 'ENV_'+hit_accession
        # Use Taxonomy ID for "metegenome".
        hit_row.tax_id = 256318

    # Create filename for the gff file
    refseq_gff_zip_filename = "{}/{}_genomic.gff.gz".format(output_folder,base_filename)

    if not (os.path.isfile(refseq_gff_zip_filename) and os.path.getsize(refseq_gff_zip_filename) > 0):
        if found_on_refseq:           
            # Build the full path to the genomic gff file
            refseq_gff_zip_ftp_path = "{}/{}_genomic.gff.gz".format(ftp_path, base_filename)

        elif BLaccess:
            #If not found at NCBI, and we have BL access, do an environmental lookup
            url = 'http://bl.biology.yale.edu/lab2/scripts/dimpl_dev.pl'
            values = {'accno' : hit_accession, 'gzip' : '1'}
            data = urllib.parse.urlencode(values).encode('ascii')
            refseq_gff_zip_ftp_path = urllib.request.Request(url, data)
        else:
            print("Accession number",hit_accession,"not found on NCBI Entrez")
            return 'ERROR'
                
        # Create gff file
        # Do the web request before creating a file
        try:
            request_gff = urllib.request.urlopen(refseq_gff_zip_ftp_path)
        except HTTPError as e:
            if e.code == 404 or e.code == 444:
                return base_filename + ' - NOT FOUND'
            elif e.code == 445:
                return base_filename + ' - NO FEATURES FOUND'
                
        with open(refseq_gff_zip_filename, 'wb') as refseq_gff_zip_file:
            #copy the content of source file to destination file.
            shutil.copyfileobj(request_gff, refseq_gff_zip_file)

        #unzip gff.gz file and save as .gff file
        input_gff = gzip.GzipFile(refseq_gff_zip_filename, 'rb')
        s_gff = input_gff.read()
        input_gff.close()

        output_gff = open("data/raw/download/{}_genomic.gff".format(base_filename), 'wb')
        output_gff.write(s_gff)
        output_gff.close()
       
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
        base_filename = download_gff(hit_row)
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

        lineage = get_taxonomy(hit_row.tax_id)
    else:
        lineage = hit_row.lineage
    
    #set range
    down_limit = start - downstream_range
    up_limit = stop + upstream_range
    
    #print statements for variables to be shown in image

    print("Match #{}".format(int(hit_row.name)+1))
    print("E-value:          " + str(e_value))
    print("%GC:              " + str(percent_gc))
    print("Score:            " + str(score))
    print("Target:           " + target_name)
    print("Genome Assembly:  " + str(base_filename))
    print("Lineage:          " + lineage)
    print("Matched Sequence: " + target_sequence)

    #skip browser link and context graph when there is no bed file found
    if 'NO FEATURES' in base_filename or 'NOT FOUND' in base_filename or 'ERROR' in base_filename:
        return ('nan','nan')
    
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
        # prerenderers get run before the track is rendered
        if start < stop:
            x1 = element.scale.topixels(start) # converting genomic coordinates to screen coordinates
            x2 = element.scale.topixels(stop)
            yield from renderer.rect(x1, 0, x2-x1, element.height, fill="lightblue", stroke="none")
        if start > stop:
            x1 = element.scale.topixels(start) # converting genomic coordinates to screen coordinates
            x2 = element.scale.topixels(stop)
            yield from renderer.rect(x1, 0, x1-x2, element.height, fill="lightblue", stroke="none")
   
    doc=genomeview.visualize_data({"":bed_file},hit_accession,down_limit,up_limit)
    cur_track = genomeview.get_one_track(doc, "")
    cur_track.prerenderers = [prerender]

    display(doc)

    return base_filename, lineage


def get_all_images(results_csv_filename, alignment):

    init_environmentals_api()
    results_df = pd.read_csv(results_csv_filename, index_col=0)
    results_df['lineage'] = results_df['lineage'].astype(str)
    results_df['assembly_accession'] = results_df['assembly_accession'].astype(str)
    updated_results_df = results_df.copy(deep=True)
    for index, row in results_df.iterrows():
        display(Markdown('---'))
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
