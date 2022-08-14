#!/bin/bash

command 2> >(sed $'s,.*,\e[31m&\e[m,'>&2)

Cleanup() {
    echo -e ${NC}
}
trap Cleanup EXIT

err_report() {
  echo -e "${RED}Encountered an error, exiting with 1." >&2
}

trap err_report ERR

if [[ $# -eq 0 ]] ; then
    echo -e "${RED}[Warning] No imput parameters provided. Displaying help message:${NC}"
    Help
    exit 0
fi

ParseGetOpts( ) {
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
    allow-overwrite)
      CX_ALLOW_OVERWITE=true
      ;;
    ?)
      Help
      exit 1
      ;;
    *)
      PARSE_BREAK=true
      ;;
  esac
}

export -f ParseGetOpts
