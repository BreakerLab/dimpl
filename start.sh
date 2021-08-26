#!/usr/bin/env bash

# ##################################################
# Script for starting docker environment for
#
version="1.0.2"               # Sets version variable
#
# Templated from https://github.com/natelandau/shell-scripts
scriptTemplateVersion="1.3.0" # Version of scriptTemplate.sh that this script is based on

#
# 2019-05-20 DATE - v1.0.0  - First Creation
# 2020-02-19 DATE - v1.0.1  - Added prompt for BL API Key
# 2021-08-25 DATE - v1.0.2  - Improved security on Windows 10.  
#                             This script now runs under WSL1&2/Ubuntu without exposing Docker port 2375.
#
# ##################################################

# Provide a variable with the location of this script.
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source Scripting Utilities
# -----------------------------------
# These shared utilities provide many functions which are needed to provide
# the functionality in this boilerplate. This script will fail if they can
# not be found.
# -----------------------------------

utilsLocation="${scriptPath}/src/lib/shell-scripts/utils.sh" #

if [ -f "${utilsLocation}" ]; then
  source "${utilsLocation}"
else
  echo "Please find the file util.sh and add a reference to it in this script. Exiting."
  exit 1
fi

# trapCleanup Function
# -----------------------------------
# Any actions that should be taken if the script is prematurely
# exited.  Always call this function at the top of your script.
# -----------------------------------
function trapCleanup() {
  echo ""
  if is_dir "${tmpDir}"; then
    rm -r "${tmpDir}"
  fi
  die "Exit trapped."  # Edit this if you like.
}

# Set Flags
# -----------------------------------
# Flags which can be overridden by user input.
# Default values are below
# -----------------------------------
quiet=0
printLog=0
verbose=0
force=0
strict=0
debug=0
configure=0
make=0
args=()

# Set Temp Directory
# -----------------------------------
# Create temp directory with three random numbers and the process ID
# in the name.  This directory is removed automatically at exit.
# -----------------------------------
tmpDir="/tmp/${scriptName}.$RANDOM.$RANDOM.$RANDOM.$$"
(umask 077 && mkdir "${tmpDir}") || {
  die "Could not create temporary directory! Exiting."
}

# Logging
# -----------------------------------
# Log is only used when the '-l' flag is set.
# -----------------------------------
logFile="$HOME/Library/Logs/${scriptBasename}.log"


function mainScript() {

# If the .env file is not present, configure
if [[ ${configure} = 1 ]]; then
    info "Overwriting any existing .env configuration file."
    configure
elif is_not_file "${scriptPath}/.env"; then
    info "No .env file found. Prompting for configuration info."
    configure
else
    info "Configured .env file found."
fi

# Make sure docker is installed then launch notebook
if type_not_exists 'docker-compose'; then
    die "docker-compose command not found. Please install docker"
else
    info "Launching/Restarting 'notebook' docker container"
    drun docker-compose pull
    drun docker-compose up -d
    sleep 2
    #info "Logging into Globus"
    #docker exec -it --user jovyan dimpl_container globus login --no-local-server
    # Generate Access URL
    token=$(drun docker exec dimpl_container sh -c '"jupyter notebook list"' |  grep -Po '\?token\=(\S+)')
    port=$(drun docker port dimpl_container | head -1 | cut -d":" -f2)
    echo "Access the GC-IGR jupyter notebook at http://localhost:$port/$token"
fi



if [[ ${make} == 1 ]]; then
    info "Running command 'make data' inside docker container"
    makeDataset
fi

safeExit
}

function configure() {

if is_file "${scriptPath}/.env"; then
  rm ./.env
fi
touch .env

input "Enter email address (Required for data requests from NCBI E-utils): "
read ENTREZ_EMAIL
if [[ ${ENTREZ_EMAIL} != "" ]]; then
    echo "ENTREZ_EMAIL=${ENTREZ_EMAIL}" >> .env
fi


input "Enter NCBI E-utils API Key (Recommended, allows faster requests from E-utils): "
read ENTREZ_APIKEY
if [[ ${ENTREZ_APIKEY} != "" ]]; then
    echo "ENTREZ_APIKEY=${ENTREZ_APIKEY}" >> .env
fi

input "Enter Breaker Lab API Key (Required for lab members, others can skip it): "
read BL_APIKEY
if [[ ${BL_APIKEY} != "" ]]; then
    echo "BL_APIKEY=${BL_APIKEY}" >> .env
fi

seek_confirmation "Customize Rfam database location? (Used when running local copies of Rfam): "
if is_confirmed; then

    input "Enter MYSQL database name: "
    read MYSQL_DATABASE;
    if [[ ${MYSQL_DATABASE} != "" ]]; then
        echo "MYSQL_DATABASE=${MYSQL_DATABASE}" >> .env
    fi

    input "Enter MYSQL host: "
    read MYSQL_HOST;
    if [[ ${MYSQL_HOST} != "" ]]; then
        echo "MYSQL_HOST=${MYSQL_HOST}" >> .env
    fi

    input "Enter MYSQL port: "
    read MYSQL_PORT;
    if [[ ${MYSQL_PORT} != "" ]]; then
        echo "MYSQL_PORT=${MYSQL_PORT}" >> .env
    fi

    input "Enter MYSQL username: "
    read MYSQL_USER;
    if [[ ${MYSQL_USER} != "" ]]; then
        echo "MYSQL_USER=${MYSQL_USER}" >> .env
    fi

    notice "Notice: Entering the MYSQL password below saves it in plaintext in the directory's .env file. \n
            For greater security, leave the line below blank and use the -p or --password option when running start.sh"
    input "Enter MYSQL password: "
    read MYSQL_PASSWORD;
    if [[ ${MYSQL_PASSWORD} != "" ]]; then
        echo "MYSQL_PASSWORD=${MYSQL_PASSWORD}" >> .env
    fi
fi
}


# Run a Docker command 
# -----------------------------------
# Determine if running on Windows (WSL), if so, run via powershell.exe
# Docker commands especially.
# -----------------------------------
function drun() {
    command=`uname -r | grep "Microsoft" | sed -r 's/.+/powershell.exe -Command /'`
    $command $@
}


function makeDataset() {

parse_yaml docker-compose.yml "DOCKER_" > "${tmpDir}/parsed_yaml.txt"
source "${tmpDir}/parsed_yaml.txt"

drun docker exec -w /home/jovyan/work --user jovyan ${DOCKER_services_notebook_container_name} make data
}

############## Begin Options and Usage ###################


# Print usage
usage() {
  echo -n "${scriptName} [OPTION]... [FILE]...

This is my script template.

 Options:
  -m, --make        Run the 'make data' command inside docker-container upon completion
  -c, --configure   Prompt for environment variables, even if already present
  -p, --password    MYSQL password
  --force           Skip all user interaction.  Implied 'Yes' to all actions.
  -q, --quiet       Quiet (no output)
  -l, --log         Print log to file
  -s, --strict      Exit script with null variables.  i.e 'set -o nounset'
  -v, --verbose     Output more information. (Items echoed to 'verbose')
  -d, --debug       Runs script in BASH debug mode (set -x)
  -h, --help        Display this help and exit
      --version     Output version information and exit
"
}

# Iterate over options breaking -ab into -a -b when needed and --foo=bar into
# --foo bar
optstring=h
unset options
while (($#)); do
  case $1 in
    # If option is of type -ab
    -[!-]?*)
      # Loop over each character starting with the second
      for ((i=1; i < ${#1}; i++)); do
        c=${1:i:1}

        # Add current char to options
        options+=("-$c")

        # If option takes a required argument, and it's not the last char make
        # the rest of the string its argument
        if [[ $optstring = *"$c:"* && ${1:i+1} ]]; then
          options+=("${1:i+1}")
          break
        fi
      done
      ;;

    # If option is of type --foo=bar
    --?*=*) options+=("${1%%=*}" "${1#*=}") ;;
    # add --endopts for --
    --) options+=(--endopts) ;;
    # Otherwise, nothing special
    *) options+=("$1") ;;
  esac
  shift
done
set -- "${options[@]}"
unset options

# Print help if no arguments were passed.
# Uncomment to force arguments when invoking the script
# [[ $# -eq 0 ]] && set -- "--help"

# Read the options and set stuff
while [[ $1 = -?* ]]; do
  case $1 in
    -h|--help) usage >&2; safeExit ;;
    --version) echo "$(basename $0) ${version}"; safeExit ;;
    -c|--configure) shift; configure=1 ;;
    -m|--make) shift; make=1 ;;
    -p|--password) shift; echo "Enter MYSQL Password: "; stty -echo; read MYSQLPASS; stty echo;
      echo ;;
    -v|--verbose) verbose=1 ;;
    -l|--log) printLog=1 ;;
    -q|--quiet) quiet=1 ;;
    -s|--strict) strict=1;;
    -d|--debug) debug=1;;
    --force) force=1 ;;
    --endopts) shift; break ;;
    *) die "invalid option: '$1'." ;;
  esac
  shift
done

# Store the remaining part as arguments.
args+=("$@")

############## End Options and Usage ###################




# ############# ############# #############
# ##       TIME TO RUN THE SCRIPT        ##
# ##                                     ##
# ## You shouldn't need to edit anything ##
# ## beneath this line                   ##
# ##                                     ##
# ############# ############# #############

# Trap bad exits with your cleanup function
trap trapCleanup EXIT INT TERM

# Exit on error. Append '||true' when you run the script if you expect an error.
set -o errexit

# Run in debug mode, if set
if [ "${debug}" == "1" ]; then
  set -x
fi

# Exit on empty variable
if [ "${strict}" == "1" ]; then
  set -o nounset
fi

# Bash will remember & return the highest exitcode in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`, for example.
set -o pipefail

# Run your script
mainScript

safeExit # Exit cleanly
