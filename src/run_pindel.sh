#!/bin/bash
set -o errexit

RED='\033[0;31m'
LCY='\033[1;36m' # Light Cyan
NC='\033[0m' # No Color

command 2> >(sed $'s,.*,\e[31m&\e[m,'>&2)

Cleanup() {
    echo -e ${NC}
}
trap Cleanup EXIT

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


if [[ $# -eq 0 ]] ; then
    echo -e "${RED}[Warning] No imput parameters provided. Displaying help message:${NC}"
    Help
    exit 0
fi

while getopts 'g:l:i:r:R:h' OPTION; do
  case "$OPTION" in
    
    g)
      echo -e "Genome file: \"$OPTARG\""
      CX_REFERENCE=$OPTARG;;
    r)
      echo -e "Reads 1 file: \"$OPTARG\""
      CX_READS_1=$OPTARG;;
    R)
      echo -e "Reads 2 file: \"$OPTARG\""
      CX_READS_2=$OPTARG;;   
    l)
      echo -e "Read length: \"$OPTARG\""
      CX_READ_LENGTH=$OPTARG;;
    h)
      Help
	  exit 0
      ;;
    ?)
      Help
      exit 1
      ;;
  esac
done

#Sanitize input
if [ -z ${CX_REFERENCE} ] || [ -z ${CX_READS_1} ] || [ -z ${CX_READS_2} ]; then
    echo -e ${RED}
    echo "[Warning] Missing required arguments: "
    [ -z ${CX_REFERENCE} ] &&   echo -e "\t[-g] Path to genome (fasta format)"
    [ -z ${CX_READS_1} ] &&     echo -e "\t[-r] Path to reads1 (fastq format)"
    [ -z ${CX_READS_2} ] &&     echo -e "\t[-R] Path to reads2 (fastq format)"
    [ -z ${CX_READ_LENGTH} ] && echo -e "\t[-l] Value of read length"
    echo -e ${NC}
    CX_BENCHMARK_INPUT_ERROR=true
fi

if [ ! -f ${CX_REFERENCE} ] || [ ! -f ${CX_READS_1} ] || [ ! -f ${CX_READS_2} ]; then
    echo -e ${RED}
    echo "[Warning] Following files do not exist: "
    [ ! -f ${CX_REFERENCE} ] && [ ! -z ${CX_REFERENCE} ] && echo -e "\t[-g] Path to genome (fasta format)"
    [ ! -f ${CX_REFERENCE} ] && [ ! -z ${CX_REFERENCE} ] && echo -e "\t     Path: ${CX_REFERENCE}"
    [ ! -f ${CX_READS_1}   ] && [ ! -z ${CX_READS_1}   ] && echo -e "\t[-r] Path to reads1 (fastq format)"
    [ ! -f ${CX_READS_1}   ] && [ ! -z ${CX_READS_1}   ] && echo -e "\t     Path: ${CX_READS_1}"
    [ ! -f ${CX_READS_2}   ] && [ ! -z ${CX_READS_2}   ] && echo -e "\t[-R] Path to reads2 (fastq format)"
    [ ! -f ${CX_READS_2}   ] && [ ! -z ${CX_READS_2}   ] && echo -e "\t     Path: ${CX_READS_2}"

    echo -e ${NC}
fi

if [ ! -z ${CX_BENCHMARK_INPUT_ERROR} ] && [ ${CX_BENCHMARK_INPUT_ERROR} = true  ]; then
    Help
    exit 1
fi

CX_REFERENCE=$(realpath ${CX_REFERENCE})
CX_READS_1=$(realpath ${CX_READS_1})
CX_READS_2=$(realpath ${CX_READS_2})
echo -e "Genome file (realpath): \"$CX_REFERENCE\""
echo -e "Read1 file (realpath): \"$CX_READS_1\""
echo -e "Read2 file (realpath): \"$CX_READS_2\""

mkdir pindel
cd pindel

echo -e ${NC}

#Alignment using bwa
export CX_SAM=alignments.bwa.sam
export CX_BAM=alignments.bwa.bam
echo -e "[+] Run bwa mem"
urun "bwa mem $CX_REFERENCE $CX_READS_1 $CX_READS_2 > $CX_SAM"
echo -e "[+] Run samtools view"
urun "samtools view -S -b $CX_SAM > $CX_BAM"

#Preprocessing with samtools
export CX_BAM_SORTED=alignments.bwa.sorted
echo -e "[+] Run samtools sort"
urun "samtools sort $CX_BAM -o $CX_BAM_SORTED"

export CX_BAM_SORTED=alignments.bwa.sorted.bam #samtools adds a .bam to the output file
export CX_BAM_BAI=alignments.bwa.sorted.bam.bai
echo -e "[+] Run samtools index"
urun "samtools index $CX_BAM_SORTED $CX_BAM_BAI"

#Running pindel
export CX_PINDEL_DIR=out.pindel
export CX_PINDEL_OUT=out.pindel
export CX_PINDEL_VCF_DIR=out.vcf
export CX_PINDEL_VCF_OUT=out.vcf
echo -e "[+] Run pindel"
mkdir -p $CX_PINDEL_DIR
urun "pindel -f $CX_REFERENCE -i <( echo $CX_BAM_SORTED $CX_READ_LENGTH sample) -o ${CX_PINDEL_DIR}/${CX_PINDEL_OUT}"

mkdir -p ${CX_PINDEL_VCF_DIR}
echo -e "[+] Run pindel2vcf"
urun "pindel2vcf -co 2 -P ${CX_PINDEL_DIR}/${CX_PINDEL_OUT} -r $CX_REFERENCE -R yeast -d 20220713 -v ${CX_PINDEL_VCF_DIR}/${CX_PINDEL_VCF_OUT}"

cd ${CX_PINDEL_VCF_DIR} || exit

echo -e "[+] Run awk - replacement"
awk ' 
/^##/ {print}
/^.*<RPL>/ {print}
' ${CX_PINDEL_VCF_OUT} > replacement.vcf

echo -e "[+] Run awk - duplication"
awk '
/^##/ {print}
/^.*<DUP:TANDEM>/ {print}
' ${CX_PINDEL_VCF_OUT} > duplication.vcf

echo -e "[+] Run awk - deletion"
awk '
/^##/ {print}
/^.*<DEL>/ {print}
' ${CX_PINDEL_VCF_OUT} > deletion.vcf

echo -e "[+] Run awk - inversion"
awk '
/^##/ {print}
/^.*<INV>/ {print}
' ${CX_PINDEL_VCF_OUT} > inversion.vcf

echo -e "[+] Run awk - insertion"
awk '
/^##/ {print}
/^.*<INS>/ {print}
' ${CX_PINDEL_VCF_OUT} > insertion.vcf

