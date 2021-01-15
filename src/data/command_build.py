import os
import tarfile
import glob
import shutil

def build_infernal_commands(data_dir, step_name, no_secondary_structure=False):

    if not os.path.exists(data_dir + '/scripts'):
        os.makedirs(data_dir + '/scripts')

    with open("src/shell/infernal_run_template.sh", 'r') as infernal_run_template: 
        #Skip the first line of the template
        infernal_run_template.readline()
        # Capture the rest of the data
        data = infernal_run_template.read()
    infernal_run_script = "{}/{}_run.sh".format(data_dir, step_name)
    with open(infernal_run_script, 'w') as infernal_run: 
        infernal_run.write("#!/bin/bash\nSTEPNAME={}\n".format(step_name) + data)

    # Make the script executable.
    os.chmod(infernal_run_script, 0o775)

    # Copy the configuration file
    shutil.copy("src/shell/cluster.conf", "{}/scripts/".format(data_dir))

    # Copy the tar helper utility
    shutil.copy("src/shell/make_tar.sh", "{}/scripts/".format(data_dir))

    # Add the outdir variable to the infernal source file
    with open("src/shell/infernal_source_template.sh", 'r') as infernal_source_template:
        data = infernal_source_template.read()
    with open("{}/scripts/{}_infernal_source.sh".format(data_dir, step_name), 'w') as infernal_source: 
        # If necessary, add line to provide the no secondary structure command to cmbuild
        if no_secondary_structure:
            infernal_source.write("NOSS='--noss'\n")
        infernal_source.write("STEPNAME={}\n".format(step_name) + data)

    # Add the outdir variable to the infernal source file
    with open("src/shell/cmfinder_source_template.sh", 'r') as cmfinder_source_template: data = cmfinder_source_template.read()
    with open("{}/scripts/{}_cmfinder_source.sh".format(data_dir, step_name), 'w') as cmfinder_source: cmfinder_source.write("STEPNAME={}\n".format(step_name) + data)
    
    # Copy the fix coords scripts to the scripts directory
    shutil.copy("src/shell/fix_sto_coords.sh", "{}/scripts/".format(data_dir))
    shutil.copy("src/shell/fix_tblout_coords.sh", "{}/scripts/".format(data_dir))
    shutil.copy("src/shell/dedupe_sample.sh", "{}/scripts/".format(data_dir))

    # Create the Output directory if it doesn't already exist
    if not os.path.exists("{}/{}".format(data_dir,step_name)):
        os.makedirs("{}/{}".format(data_dir,step_name))

    igr_list = glob.glob("{}/{}/*/".format(data_dir,step_name))
    motifs_dir = "{}/{}/motifs/".format(data_dir,step_name)
    if motifs_dir in igr_list:
        igr_list.remove(motifs_dir)
        
    with open("{}/scripts/{}_infernal_jobfile.sh".format(data_dir, step_name), 'w') as jobfile:
        for igr_filename in igr_list:
            igr_name = igr_filename[:-1].split('/')[-1]
            jobfile.write("SEQNAME={}; source scripts/{}_infernal_source.sh; $CMBUILDCOMMAND && $CMCALIBRATECOMMAND && $CMSEARCHCOMMAND\n".format(igr_name, step_name))
            
    with open("{}/scripts/{}_cmfinder_jobfile.sh".format(data_dir,step_name), 'w') as jobfile:
        for igr_filename in igr_list:
            igr_name = igr_filename[:-1].split('/')[-1]
            jobfile.write("SEQNAME={}; source scripts/{}_cmfinder_source.sh; bash -c \"$CMFINDERCOMMAND\" \n".format(igr_name, step_name))
            
