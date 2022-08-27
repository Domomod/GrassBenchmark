#!/bin/bash


. common_color.sh

CX_GRIDSS_DIR=gridss/results
CX_GRASS_DIR=grass_sv/results
CX_PINDEL_DIR=pindel/out.bed
CX_LUMPY_DIR=lumpy/out.bed


for DATA_DIR in $(ls -d fastadna_*); do
    (
    cd $DATA_DIR || exit 1
    echo -e $B_GRN "======================================================================================" $NC
    echo -e $B_GRN "[PROCESSING] ${DATA_DIR}" $NC
    echo -e $B_GRN "======================================================================================\n" $NC

    cd generated_mutations
    . export_breakpoints.sh
    cd ..

    echo -e $LCY   \
    '##########\n' \
    '# GRIDSS #\n' \
    '##########\n' $NC

    GrassSV.py utils check_sv -g generated_mutations -d ${CX_GRIDSS_DIR} --bp_mode

    echo -e $LCY   \
    '##########\n' \
    '# LUMPY  #\n' \
    '##########\n' $NC

    GrassSV.py utils check_sv -g generated_mutations -d ${CX_LUMPY_DIR}
    GrassSV.py utils check_sv -g generated_mutations -d ${CX_LUMPY_DIR} --bp_mode


    echo -e $LCY   \
    '##########\n' \
    '# PINDEL #\n' \
    '##########\n' $NC

    GrassSV.py utils check_sv -g generated_mutations -d ${CX_PINDEL_DIR}
    GrassSV.py utils check_sv -g generated_mutations -d ${CX_PINDEL_DIR} --bp_mode


    echo -e $LCY   \
    '##########\n' \
    '# GRASS  #\n' \
    '##########\n' $NC

    GrassSV.py utils check_sv -g generated_mutations -d ${CX_GRASS_DIR}
    GrassSV.py utils check_sv -g generated_mutations -d ${CX_GRASS_DIR} --bp_mode

    ) | tee out/$DATA_DIR.check_all.log.temp
    more out/$DATA_DIR.check_all.log.temp | sed 's/\x1b\[[0-9;]*[mGKH]//g'   > out/$DATA_DIR.check_all.log
done

cd out
rm *.temp