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
   echo -e "${LCY}This will run all steps for manta pipeline in a single process."
   echo -e "Prouced output will be stored in ${PWD}/manta" 
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

#Exit on error
set -o errexit

if [[ -d manta  ]] && [[ ! -z ${CX_ALLOW_OVERWITE+x} ]]; then
  echo -e "${YELLOW}[Warning] Flag \"allow-overwite\" set to true, recurisvely removing $PWD/manta directory.${NC}" 
  rm -rf manta
fi

mkdir manta
d manta


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

#Running manta
export MANTA_CONFIG=manta/bin/configManta.py

echo -e "${LCY}[+] Run manta${NC}"
python $MANTA_CONFIG --bam $CX_BAM_SORTED --referenceFasta $CX_REFERENCE --runDir .
python runWorkflow.py

#Removing trash, unzipping results and sending them to the main folder
rm *orkflow.*
rm -rf workspace
gzip -dv results/variants/candidateSV.vcf.gz
cp results/variants/candidateSV.vcf results.vcf

#AWK - TO DO

echo -e "${LCY}[+] Run awk - replacement${NC}"
awk '
/^#/ {print}
/^.*<RPL>/ {print}
' ${CX_PINDEL_VCF_OUT} > replacement.vcf

echo -e "${LCY}[+] Run awk - duplication${NC}"
awk '
/^#/ {print}
/^.*<DUP:TANDEM>/ {print}
' ${CX_PINDEL_VCF_OUT} > duplication.vcf

echo -e "${LCY}[+] Run awk - deletion${NC}"
awk '
/^#/ {print}
/^.*<DEL>/ {print}
' ${CX_PINDEL_VCF_OUT} > deletion.vcf

echo -e "${LCY}[+] Run awk - inversion${NC}"
awk '
/^#/ {print}
/^.*<INV>/ {print}
' ${CX_PINDEL_VCF_OUT} > inversion.vcf

echo -e "${LCY}[+] Run awk - insertion${NC}"
awk '
/^#/ {print}
/^.*<INS>/ {print}
' ${CX_PINDEL_VCF_OUT} > insertion.vcf

echo -e "${LCY}[+] Step out of ${CX_PINDEL_VCF_DIR} to $(realpath ..) ${NC}"
cd .. || exit

shopt -s extglob


echo -e "${LCY}[+] Converting ${CX_PINDEL_VCF_DIR} to ${CX_PINDEL_VCF_DIR} ${NC}"

for FILE in $(cd ${CX_PINDEL_VCF_DIR}; ls !(out).vcf); do

    NAME=${FILE%.vcf}.bed
    GrassSV.py utils csv2bed -i ${CX_PINDEL_VCF_DIR}/${FILE} -o ${CX_PINDEL_BED_DIR}/${NAME}
done
