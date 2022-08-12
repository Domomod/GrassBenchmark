#!/bin/bash

RED='\033[0;31m'
LCY='\033[1;36m' # Light Cyan
NC='\033[0m' # No Color

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


if [[ $# -eq 0 ]] ; then
    echo -e "${RED}[Warning] No imput parameters provided. Displaying help message:${NC}"
    Help
    exit 0
fi

while getopts 'g:l:i:r:R:h' OPTION; do
  case "$OPTION" in
    
    g)
      echo "Genome file: \"$OPTARG\""
      CX_REFERENCE=$OPTARG;;
    r)
      echo "Reads 1 file: \"$OPTARG\""
      CX_READS_1=$OPTARG;;
    R)
      echo "Reads 2 file: \"$OPTARG\""
      CX_READS_2=$OPTARG;;   
    l)
      echo "Read length: \"$OPTARG\""
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

if [ -f ${CX_REFERENCE} ] || [ -f ${CX_READS_1} ] || [ -f ${CX_READS_2} ]; then
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

if [ ${CX_BENCHMARK_INPUT_ERROR} = true  ]; then
    Help
    exit 1
fi

CX_REFERENCE=$(realpath ${CX_REFERENCE})
CX_READS_1=$(realpath ${CX_READS_1})
CX_READS_2=$(realpath ${CX_READS_2})

#Exit on error
set -o errexit

mkdir lumpy
cd lumpy || exit 1

urun "speedseq align -R '@RG\tID:id\tSM:sample\tLB:lib' $CX_REFERENCE $CX_READS_1 $CX_READS_2"

export CX_READS_BAM=${CX_READS_1##*/}.fq.bam
export CX_DISCORDANTS=${CX_READS_1##*/}.splitters.bam
export CX_SPLITTERS=${CX_READS_1##*/}.discordants.bam
export CX_LUMPY_VCF=lumpy.vcf

urun "lumpyexpress -B ${CX_READS_BAM} -S ${CX_SPLITTERS} -D ${CX_DISCORDANTS} -o ${CX_LUMPY_VCF}"

cd lumpy