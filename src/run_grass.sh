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
   echo -e "${LCY}This will run all steps for lumpy pipeline in a single process."
   echo -e "Prouced output will be stored in ${PWD}/lumpy"
   echo -e 
   echo -e "Syntax: $0 [-g|r|R|l|h]"
   echo -e "options:"
   echo -e "g     Path to genome (fasta format)"
   echo -e "I     Path to genome index (fai format)"
   echo -e "L     Path to genome lengths (length format)"
   echo -e "r     Path to reads1 (fastq format)"
   echo -e "R     Path to reads2 (fastq format)"
   echo -e "l     Value of read length"
   echo -e "h     Show this helpe message"
   echo -e ${NC}
}

. common_parse.sh

#[CONFUSING] getopt's ':' means "no argument" for LONG options, and "requires argument" for SHORT option 
LONG_LIST=("allow-overwrite:")
SHORT_LIST=("g:l:i:r:R:h")

opts=$(getopt \
  --longoptions "$(printf "%s:," "${LONG_LIST[@]}")" \
  --name "$(basename "$0")" \
  --options "${SHORT_LIST[@]}" \
  -- "$@"
) || exit

eval set --$opts

while [[ $# -gt 0 ]]; do
    OPTION=${1##-}
    OPTION=${OPTION##-}
    OPTARG=$2

    ParseGetOpts
    
    if [[ ! -z ${PARSE_BREAK} ]] && [[ ${PARSE_BREAK} == true ]]; then
      break
    fi
    shift 2
done

. common_sanitize.sh
CX_INDEX=${CX_REFERENCE}.fai
CX_REFERENCE_LENGTH=${CX_REFERENCE%.*}.lengths
echo -e "Genome index file (realpath): \"$CX_INDEX\""
echo -e "Genome lengths file (realpath): \"$CX_REFERENCE_LENGTH\""

#Exit on error
set -o errexit

if [[ -d grass_sv  ]] && [[ ! -z ${CX_ALLOW_OVERWITE+x} ]]; then
  echo -e "${YELLOW}[Warning] Flag \"allow-overwite\" set to true, recurisvely removing $PWD/grass_sv directory.${NC}"
  rm -rf grass_sv
fi 

echo -e "${LCY}[+] Creating grass_sv directory [PWD is \"${PWD}\"]${NC}"
mkdir grass_sv

echo -e "${LCY}[+] Going into grass_sv directory [PWD is \"${PWD}\"]${NC}"
cd grass_sv

echo -e "${LCY}[+] Run GrassSV.py run_standalone [PWD is \"${PWD}\"]${NC}"
GrassSV.py run_standalone -g ${CX_REFERENCE} -i ${CX_INDEX} -l ${CX_REFERENCE_LENGTH} -r ${CX_READS_1} -R ${CX_READS_2} 

echo -e "${LCY}[+] Cleanup results directory [PWD is \"${PWD}\"]${NC}"
mv results/detectedSVs/* results/
rmdir results/detectedSVs
mv results/filter_inversions.bed results/inversions.bed

echo -e "${LCY}[+] Run GrassSV.py utils check_sv [PWD is \"${PWD}\"]${NC}"
GrassSV.py utils check_sv -g ../generated_mutations/ -d results/ > grass_sv.results.quality