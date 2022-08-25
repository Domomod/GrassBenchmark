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
   echo -e "${LCY}This will run all steps for pindel pipeline in a single process."
   echo -e "Prouced output will be stored in ${PWD}/lumpy"
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

if [[ -d pindel  ]] && [[ ! -z ${CX_ALLOW_OVERWITE+x} ]]; then
  echo -e "${YELLOW}[Warning] Flag \"allow-overwite\" set to true, recurisvely removing $PWD/lumpy directory.${NC}"
  rm -rf pindel
fi 

#mkdir pindel
cd pindel


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

#Running pindel
export CX_PINDEL_DIR=out.pindel
export CX_PINDEL_OUT=out.pindel
export CX_PINDEL_VCF_DIR=out.vcf
export CX_PINDEL_VCF_OUT=out.vcf
export CX_PINDEL_BED_DIR=out.bed

echo -e "${LCY}[+] Run pindel${NC}"
#mkdir -p $CX_PINDEL_DIR
#urun "pindel -f $CX_REFERENCE -i <( echo $CX_BAM_SORTED $CX_READ_LENGTH sample) -o ${CX_PINDEL_DIR}/${CX_PINDEL_OUT}"
#mkdir -p ${CX_PINDEL_VCF_DIR}
echo -e "${LCY}[+] Run pindel2vcf${NC}"
#urun "pindel2vcf -co 2 -P ${CX_PINDEL_DIR}/${CX_PINDEL_OUT} -r $CX_REFERENCE -R yeast -d 20220713 -v ${CX_PINDEL_VCF_DIR}/${CX_PINDEL_VCF_OUT}"

echo -e "${LCY}[+] Step into ${CX_PINDEL_VCF_DIR} ${NC}"
cd ${CX_PINDEL_VCF_DIR} || exit

echo -e "${LCY}[+] Run awk - replacement${NC}"
awk ' 
/^#/ {print}
/SVTYPE=RPL/ {print}
' ${CX_PINDEL_VCF_OUT} > replacements.vcf

echo -e "${LCY}[+] Run awk - duplication${NC}"
awk '
/^#/ {print}
/SVTYPE=DUP/ {print}
' ${CX_PINDEL_VCF_OUT} > duplications.vcf

echo -e "${LCY}[+] Run awk - deletion${NC}"
awk '
/^#/ {print}
/SVTYPE=DEL/ {print}
' ${CX_PINDEL_VCF_OUT} > deletions.vcf

echo -e "${LCY}[+] Run awk - inversion${NC}"
awk '
/^#/ {print}
/SVTYPE=INV/ {print}
' ${CX_PINDEL_VCF_OUT} > inversions.vcf

echo -e "${LCY}[+] Run awk - insertion${NC}"
awk '
/^#/ {print}
/SVTYPE=INS/ {print}
' ${CX_PINDEL_VCF_OUT} > insertions.vcf

echo -e "${LCY}[+] Step out of ${CX_PINDEL_VCF_DIR} to $(realpath ..) ${NC}"
cd .. || exit

shopt -s extglob


echo -e "${LCY}[+] Converting ${CX_PINDEL_VCF_DIR} to ${CX_PINDEL_VCF_DIR} ${NC}"

rm -f ${CX_PINDEL_BED_DIR}/breakpoints.bed
for FILE in $(cd ${CX_PINDEL_VCF_DIR}; ls !(out).vcf); do
    
    NAME=${FILE%.vcf}.bed    
    GrassSV.py utils csv2bed -i ${CX_PINDEL_VCF_DIR}/${FILE} -o ${CX_PINDEL_BED_DIR}/${NAME}
    awk -v record=${FILE%.vcf} '
    /^.*/ {printf  "%s %s %s %s_%s_l\n%s %s %s %s_%s_r\n", $1, $2, $2+1, record, $4, $1, $3, $3+1, record, $4}
    ' ${CX_PINDEL_BED_DIR}/${NAME} >> ${CX_PINDEL_BED_DIR}/breakpoints.bed
done 