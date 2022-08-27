#!/bin/bash

RED='\033[0;31m'
LCY='\033[1;36m' # Light Cyan
NC='\033[0m' # No Color


. common_color.sh

Help()
{
   # Display Help
   echo -e "${LCY}This will run all steps for gridss pipeline in a single process."
   echo -e "Prouced output will be stored in ${PWD}/lumpy"
   echo -e 
   echo -e "Syntax: $0 [-g|r|R|l|h]"
   echo -e "options:"
   echo -e "g     Path to genome (fasta format)"
   echo -e "r     Path to reads1 (fastq format)"
   echo -e "R     Path to reads2 (fastq format)"
   echo -e "h     Show this helpe message"
   echo -e ${NC}
}

. common_parse.sh

#[CONFUSING] getopt's ':' means "no argument" for LONG options, and "requires argument" for SHORT option 
LONG_LIST=("allow-overwrite:")
SHORT_LIST=("g:l:r:R:h")

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

#Exit on error
set -o errexit

if [[ -d gridss  ]] && [[ ! -z ${CX_ALLOW_OVERWITE+x} ]]; then
  echo -e "${YELLOW}[Warning] Flag \"allow-overwite\" set to true, recurisvely removing $PWD/gridss directory.${NC}"
  rm -rf gridss
fi 

mkdir gridss || exit 1
cd gridss || exit 1

# TODO: add running BWA

#Alignment using bwa
export CX_SAM=alignments.bwa.sam
export CX_BAM=alignments.bwa.bam
echo -e "${LCY}[+] Run bwa mem${NC}"
urun "bwa mem $CX_REFERENCE $CX_READS_1 $CX_READS_2 > $CX_SAM"
echo -e "${LCY}[+] Run samtools view${NC}"
urun "samtools view -S -b $CX_SAM > $CX_BAM"

#Preprocessing with samtools
export CX_BAM_SORTED=alignments.bwa.sorted
echo -e "${LCY}[+] Run samtools sort${NC}"
urun "samtools sort $CX_BAM -o $CX_BAM_SORTED"

export CX_BAM_SORTED=alignments.bwa.sorted.bam #samtools adds a .bam to the output file
export CX_BAM_BAI=alignments.bwa.sorted.bam.bai
echo -e "${LCY}[+] Run samtools index${NC}"
urun "samtools index $CX_BAM_SORTED $CX_BAM_BAI"

mkdir -p workdir && cd workdir || exit 1

export CX_GRIDSS_OUT=gridss.out.vcf.gz
urun "gridss --reference ${CX_REFERENCE} --output ${CX_GRIDSS_OUT} ../${CX_BAM_SORTED}"
gzip -d ${CX_GRIDSS_OUT}
export CX_GRIDSS_OUT=${CX_GRIDSS_OUT%.gz}


cd ..
mkdir -p results

awk '
/SVTYPE=BND/ {printf  "%s %s %s %s\n", $1, $2, $2+1, $3}
' workdir/${CX_GRIDSS_OUT} > results/breakpoints.bed
