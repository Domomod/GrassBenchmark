#!/bin/bash

RED='\033[0;31m'
LCY='\033[1;36m' # Light Cyan
NC='\033[0m' # No Color

Help()
{
   # Display Help
   echo -e "${LCY}This will run all steps for gridss pipeline in a single process."
   echo -e "Produced output will be stored in ${PWD}/gridss"
   echo -e 
   echo -e "Syntax: $0 [-g|a|r|R|s|o|h] - provide either (-r and -R) or -s"
   echo -e "options:"
   echo -e "g     Path to genome (fasta format)"
   echo -e "a     Assembly (bam format)"
   echo -e "r     Path to reads1 (fastq format)"
   echo -e "R     Path to reads2 (fastq format)"
   echo -e "s     Alignment sorted (bam format)"
   echo -e "o     Output file name (vcf.gz format)"
   echo -e "h     Show this helpe message"
   echo -e ${NC}
}


if [[ $# -eq 0 ]] ; then
    echo -e "${RED}[Warning] No input parameters provided. Displaying help message:${NC}"
    Help
    exit 0
fi

while getopts 'g:a:r:R:s:o:h' OPTION; do
  case "$OPTION" in
    g)
      echo "Genome file: \"$OPTARG\""
      GENOME_FILE_PATH=$OPTARG;;
    a)
      echo "Assembly file: \"$OPTARG\""
      ASSEMBLY_FILE_PATH=$OPTARG;;
    r)
      echo "Reads1 file: \"$OPTARG\""
      READS1_FILE_PATH=$OPTARG;;
    R)
      echo "Reads2 file: \"$OPTARG\""
      READS2_FILE_PATH=$OPTARG;;
    s)
      echo "Alignment sorted file: \"$OPTARG\""
      ALIGNMENT_FILE_PATH=$OPTARG;;
    o)
      echo "Output file name: \"$OPTARG\""
      OUTPUT_NAME=$OPTARG;;
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
if [ -z ${GENOME_FILE_PATH} ] || [ -z ${ALIGNMENT_FILE_PATH} ] || [ -z ${READS1_FILE_PATH} ] || [ -z ${READS2_FILE_PATH} ] || [ -z ${ASSEMBLY_FILE_PATH} ] || [ -z ${OUTPUT_NAME} ]; then
    echo -e ${RED}
    echo "[Warning] Missing required arguments: "
    [ -z ${GENOME_FILE_PATH} ]           && echo -e "\t[-g] Path to genome (fasta format)"          && CX_BENCHMARK_INPUT_ERROR=true
    [ -z ${ASSEMBLY_FILE_PATH} ]         && echo -e "\t[-a] Assembly (bam format)"                  && CX_BENCHMARK_INPUT_ERROR=true
    [ -z ${READS1_FILE_PATH} ]           && echo -e "\t[-r] Path to reads1 (fasta format)"          && CX_BENCHMARK_INPUT_WARN=true
    [ -z ${READS2_FILE_PATH} ]           && echo -e "\t[-R] Path to reads2 (fasta format)"          && CX_BENCHMARK_INPUT_WARN=true
    [ -z ${ALIGNMENT_FILE_PATH} ]        && echo -e "\t[-s] Path to alignment sorted (bam format)"  && CX_BENCHMARK_INPUT_WARN=true
    [ -z ${OUTPUT_NAME} ]                && echo -e "\t[-o] Output file name (vcf.gz format)"       && CX_BENCHMARK_INPUT_ERROR=true
    echo -e ${NC}
fi

if [ -f ${GENOME_FILE_PATH} ] || [ -f ${ASSEMBLY_FILE_PATH} ] || [ -f ${ALIGNMENT_FILE_PATH} ] || [ -f ${READS1_FILE_PATH} ] || [ -f ${READS2_FILE_PATH} ]; then
    echo -e ${RED}
    echo "[Warning] Following files do not exist: "
    [ ! -f ${GENOME_FILE_PATH} ] && [ ! -z ${GENOME_FILE_PATH} ] \
        && echo -e "\t[-g] Path to genome (fasta format)" \
        && echo -e "\t\t\tPath: ${GENOME_FILE_PATH}"

    [ ! -f ${ALIGNMENT_FILE_PATH} ] && [ ! -z ${ALIGNMENT_FILE_PATH}   ] \
        && echo -e "\t[-af] Path to alignment (bam format)" \
        && echo -e "\t\t\tPath: ${ALIGNMENT_FILE_PATH}"

    [ ! -f ${READS1_FILE_PATH} ] && [ ! -z ${READS1_FILE_PATH} ] \
        && echo -e "\t[-r1] Path to reads1 (fasta format)" \
        && echo -e "\t\t\tPath: ${READS1_FILE_PATH}"

    [ ! -f ${READS2_FILE_PATH} ] && [ ! -z ${READS2_FILE_PATH} ] \
        && echo -e "\t[-r2] Path to reads2 (fasta format)" \
        && echo -e "\t\t\tPath: ${READS2_FILE_PATH}"

    [ ! -f ${ASSEMBLY_FILE_PATH}   ] && [ ! -z ${ASSEMBLY_FILE_PATH} ] \
        && echo -e "\t[-as] Path to assembly (bam format)" \
        && echo -e "\t\t\tPath: ${ASSEMBLY_FILE_PATH}"
    
    echo -e ${NC} 
fi

if [ ${CX_BENCHMARK_INPUT_ERROR} ]; then
    Help
    echo -e "Stopping"
    exit 1
fi

if [ ${CX_BENCHMARK_INPUT_WARN} ]; then
    Help
    echo -e "Continuing..."
fi


GENOME_FILE_PATH=$(realpath ${GENOME_FILE_PATH})
ALIGNMENT_FILE_PATH=$(realpath ${ALIGNMENT_FILE_PATH})
READS1_FILE_PATH=$(realpath ${READS1_FILE_PATH})
READS2_FILE_PATH=$(realpath ${READS2_FILE_PATH})
ASSEMBLY_FILE_PATH=$(realpath ${ASSEMBLY_FILE_PATH})
OUTPUT_FILE_PATH=$(realpath ${OUTPUT_NAME})

echo -e "Genome file (realpath): \"$GENOME_FILE_PATH\""
echo -e "Alignment file (realpath): \"$ALIGNMENT_FILE_PATH\""
echo -e "Read1 file (realpath): \"$READS1_FILE_PATH\""
echo -e "Read2 file (realpath): \"$READS2_FILE_PATH\""
echo -e "Assembly file (realpath): \"$ASSEMBLY_FILE_PATH\""
echo -e "Output file (realpath): \"$OUTPUT_FILE_PATH\""

#Exit on error
set -o errexit

mkdir gridss && cd gridss || exit 1

# TODO: add running BWA

./urun "gridss --reference ${GENOME_FILE_PATH} --output ${OUTPUT_FILE_PATH} --assembly ${ALIGNMENT_FILE_PATH}"
