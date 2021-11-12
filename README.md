# DIMPL (Discovery of Intergenic Motifs PipeLine)
==============================

### Summary

The DIMPL discovery pipeline enables rapid extraction and selection of bacterial IGRs that are enriched for structured ncRNAs. DIMPL also automates the subsequent computational steps necessary for their functional identification.

### Requirements

#### For Local Computer

* [Docker Desktop/Engine](https://hub.docker.com/search?q=&type=edition&offering=community&sort=updated_at&order=desc)

#### For Compute Cluster

* [Infernal v1.1.4](http://eddylab.org/infernal/)
* [BLAST+ v2.11.0+](https://www.ncbi.nlm.nih.gov/books/NBK279680/)
* [CMfinder v0.4.1.18+](https://sourceforge.net/projects/weinberg-cmfinder/)
* [Dead Simple Queue v1.05+](https://github.com/ycrc/dSQ)
* [Slurm v20.02+](https://slurm.schedmd.com/quickstart.html)

### Quick-start

#### Cluster Configuration

1. Download the IGR search database (filename: s50.igr.fasta) from [this link](https://app.globus.org/file-manager?origin_id=347584ae-43bf-11ec-a6bf-9b4f84e67de8&origin_path=%2F) to your cluster using [Globus FTP](https://www.globus.org/).

2. Ensure the availability of the BLAST nr database on your computational cluster. Follow [these instructions](https://www.ncbi.nlm.nih.gov/books/NBK537770/) for updating/downloading the latest version.

#### Local Configuration

1. Download the source code (into any folder on drive)

`wget https://github.com/BreakerLab/dimpl/archive/v1.0.1.tar.gz`
`tar xzvf v1.0.1.tar.gz`

2. Download the docker image. 

`docker pull breakerlab/dimpl`

3. Configure docker to grant containers access to the folder where the DIMPL repository is located

4. Modify the configuration file found at dimpl/src/shell/cluster.conf with the database locations and appropriate commands for importing utilities on your cluster. 

5. Run `./start.sh` in the main repository directory. Follow the first-time configuration instructions (asks for email and NCBI API key).

6. Follow the link generated by the `start.sh` script to access the DIMPL jupyter notebooks.

#### Data Transfer between Local Machine and Cluster

The DIMPL notebooks generate compressed .tar.gz files consisting of all the scripts and data necessary to run the more computationally demanding steps on a cluster. These .tar.gz files are placed in the directory data/export. After transferring the files to a cluster they should be unpacked using the command `tar xzvf data-dir.tar.gz`. When tasks on the cluster complete the directory should be recompressed using the command `tar czvf data-dir.tar.gz data-dir`.

### File Organization
------------

    ├── .env                    <- File generated during configuration step of start.sh
    ├── LICENSE
    ├── README.md               <- This document
    ├── start.sh                <- Script to perform initial configuration and start the docker container
    ├── data
    │   ├── export              <- Where DIMPL places data and bash script tar.gz files  
    │   ├── import              <- Where to place re-compressed tar.gz files that have been run on a compute cluster
    │   ├── interim             <- Where processed genomic data is stored during analysis
    │   └── raw                 <- The original genomic data.
    │
    ├── docs                    <- Sphinx documentation for DIMPL
    │
    ├── notebooks               <- Jupyter notebooks for the various steps of DIMPL
    │   ├── 1-Genome-IGR-Selection.ipynb    <- 
    │   ├── 2-BLAST-Processing.ipynb        <- 
    │   ├── 3-IGR-Report.ipynb              <- 
    │   └── 4-Motif-Refinement.ipynb        <- 
    │
    ├── requirements.txt        <- The requirements file for reproducing the analysis environment, e.g.
    │                              generated with `pip freeze > requirements.txt`
    │
    ├── setup.py                <- makes project pip installable (pip install -e .) so src can be imported
    └── src                     <- Source code for use in this project.
        └── shell               <- The original genomic data.
            └── cluster.conf    <- Configuration file for the compute environment.
