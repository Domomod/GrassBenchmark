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

#Exit on error
set -o errexit

#if [[ -d lumpy  ]] && [[ ! -z ${CX_ALLOW_OVERWITE+x} ]]; then
#  echo -e "${YELLOW}[Warning] Flag \"allow-overwite\" set to true, recurisvely removing $PWD/lumpy directory.${NC}"
#  rm -rf lumpy
#fi 

echo -e "${LCY}[+] Creating lumpy directory [PWD is \"${PWD}\"]${NC}"
#mkdir lumpy

echo -e "${LCY}[+] Going into lumpy directory [PWD is \"${PWD}\"]${NC}"
cd lumpy

echo -e "${LCY}[+] Running speedseq align [PWD is \"${PWD}\"]${NC}"
#urun "speedseq align -R '@RG\tID:id\tSM:sample\tLB:lib' $CX_REFERENCE $CX_READS_1 $CX_READS_2"

export CX_READS_BAM=${CX_READS_1##*/}.bam
export CX_DISCORDANTS=${CX_READS_1##*/}.splitters.bam
export CX_SPLITTERS=${CX_READS_1##*/}.discordants.bam
export CX_LUMPY_VCF_OUT=out.vcf
export CX_LUMPY_VCF_DIR=out.vcf
export CX_LUMPY_BED_DIR=out.bed

echo -e "${LCY}[+] Running lumpyexpress [PWD is \"${PWD}\"]${NC}"
mkdir -p ${CX_LUMPY_VCF_DIR} 
#urun "lumpyexpress -B ${CX_READS_BAM} -S ${CX_SPLITTERS} -D ${CX_DISCORDANTS} -o ${CX_LUMPY_VCF_DIR}/${CX_LUMPY_VCF_OUT}"

echo -e "${LCY}[+] Step into ${CX_LUMPY_VCF_DIR}  [PWD is \"${PWD}\"]${NC}"
cd ${CX_LUMPY_VCF_DIR} || exit

echo -e "${LCY}[+] Run awk - bnd [PWD is \"${PWD}\"]${NC}"
awk ' 
/^#/ {print}
/^.*SVTYPE=BND/ {print}
' ${CX_LUMPY_VCF_OUT} > bnds.vcf

#echo -e "${LCY}[+] Run awk - duplication${NC}"
#awk '
#/^#/ {print}
#/^.*SVTYPE=DUP:TANDEM/ {print}
#' ${CX_LUMPY_VCF_OUT} > duplications.vcf

echo -e "${LCY}[+] Run awk - duplication${NC}"
awk '
/^#/ {print}
/^.*SVTYPE=DUP/ {print}
' ${CX_LUMPY_VCF_OUT} > duplications.vcf

echo -e "${LCY}[+] Run awk - deletion${NC}"
awk '
/^#/ {print}
/^.*SVTYPE=DEL/ {print}
' ${CX_LUMPY_VCF_OUT} > deletions.vcf

echo -e "${LCY}[+] Run awk - inversion${NC}"
awk '
/^#/ {print}
/^.*SVTYPE=INV/ {print}
' ${CX_LUMPY_VCF_OUT} > inversions.vcf

echo -e "${LCY}[+] Run awk - insertion${NC}"
awk '
/^#/ {print}
/^.*SVTYPE=INS/ {print}
' ${CX_LUMPY_VCF_OUT} > insertions.vcf

echo -e "${LCY}[+] Run awk - insertion${NC}"
awk '
/^#/ {print}
/^.*SVTYPE=CNV/ {print}
' ${CX_LUMPY_VCF_OUT} > copy_numbers.vcf

echo -e "${LCY}[+] Step out of ${CX_LUMPY_VCF_DIR} to $(realpath ..) ${NC}"
cd .. || exit

shopt -s extglob


echo -e "${LCY}[+] Converting ${CX_LUMPY_VCF_DIR} to ${CX_PINDEL_VCF_DIR} ${NC}"

for FILE in $(cd ${CX_LUMPY_VCF_DIR}; ls !(out).vcf); do
    
    NAME=${FILE%.vcf}.bed    
    GrassSV.py utils csv2bed -i ${CX_LUMPY_VCF_DIR}/${FILE} -o ${CX_LUMPY_BED_DIR}/${NAME}
done 

echo 
echo -e "${B_GRN}[SUCCESS] Running $0 succesfull ${NC}"
