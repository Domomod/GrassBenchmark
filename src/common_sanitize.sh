#!/bin/bash

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
    CX_BENCHMARK_INPUT_ERROR=true
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