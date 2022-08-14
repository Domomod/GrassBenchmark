#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
if [[ ":${PATH}:" != *":${SCRIPT_DIR}:"* ]]; then
  echo -e "${B_RED}[ERROR] Please add \"${SCRIPT_DIR}\" to your path and run the script again.${NC}"
  exit 1
fi

. common_color.sh

Help()
{
   # Display Help
   echo -e "${LCY}This will run all the programs included in the."
   echo -e "Prouced output will be stored in pwd: ${PWD}/lumpy"
   echo -e 
   echo -e "Syntax: $0 [-g|r|R|l|h]"
   echo -e "options:"
   echo -e "g     Path to genome (fasta format)"
   echo -e "r     Path to reads1 (fastq format)"
   echo -e "R     Path to reads2 (fastq format)"
   echo -e "l     Value of read length"
   echo -e "h     Show this help message"
   echo -e ${NC}
}


. common_parse.sh

while getopts 'g:l:i:r:R:h' OPTION; do
    ParseGetOpts ${OPTION} 
done

. common_sanitize.sh

shopt -s extglob

for SCRIPT in $(cd ${SCRIPT_DIR}; ls run_!(all).sh); do
    C=${SCRIPT%.sh}
    echo -e "${LCY}[+] Calling ${SCRIPT}${NC}"

    ${SCRIPT} "$@"
done