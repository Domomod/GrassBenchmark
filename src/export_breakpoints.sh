mkdir -p breakpoints
mkdir -p breakpoints/margin1
mkdir -p breakpoints/margin5
mkdir -p breakpoints/margin10
mkdir -p breakpoints/margin70

awk -v m=1 '
/TRANS:TO/ {print}
' translocations.bed > translocations.to.bed

awk -v m=1 '
/TRANS:FROM/ {print}
' translocations.bed > translocations.from.bed

#
# Translocation split
#

awk -v m=1 '
/TRANS:FROM/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
/TRANS:TO/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' translocations.bed > breakpoints/margin1/translocation.breakpoints.bed

awk -v m=1 '
/INS/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' insertions.bed > breakpoints/margin1/insertions.breakpoints.bed

awk -v m=1 '
/INV/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' inversions.bed > breakpoints/margin1/inversions.breakpoints.bed

awk -v m=1 '
/DEL/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' deletions.bed > breakpoints/margin1/deletions.breakpoints.bed

awk -v m=1 '
/DUP/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' duplications.bed > breakpoints/margin1/duplications.breakpoints.bed

#
# 5
#

awk -v m=5 '
/TRANS:FROM/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
/TRANS:TO/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' translocations.bed > breakpoints/margin5/translocation.breakpoints.bed

awk -v m=5 '
/INS/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' insertions.bed > breakpoints/margin5/insertions.breakpoints.bed

awk -v m=5 '
/INV/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' inversions.bed > breakpoints/margin5/inversions.breakpoints.bed

awk -v m=5 '
/DEL/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' deletions.bed > breakpoints/margin5/deletions.breakpoints.bed

awk -v m=5 '
/DUP/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' duplications.bed > breakpoints/margin5/duplications.breakpoints.bed

#
# 10
#


awk -v m=10 '
/TRANS:FROM/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
/TRANS:TO/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' translocations.bed > breakpoints/margin10/translocation.breakpoints.bed

awk -v m=10 '
/INS/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' insertions.bed > breakpoints/margin10/insertions.breakpoints.bed

awk -v m=10 '
/INV/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' inversions.bed > breakpoints/margin10/inversions.breakpoints.bed

awk -v m=10 '
/DEL/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' deletions.bed > breakpoints/margin10/deletions.breakpoints.bed

awk -v m=10 '
/DUP/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' duplications.bed > breakpoints/margin10/duplications.breakpoints.bed

#
# 70
#


awk -v m=70 '
/TRANS:FROM/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
/TRANS:TO/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' translocations.bed > breakpoints/margin70/translocation.breakpoints.bed

awk -v m=70 '
/INS/ {printf  "%s %s %s %s\n", $1, $2-m, $3+m, $4}
' insertions.bed > breakpoints/margin70/insertions.breakpoints.bed

awk -v m=70 '
/INV/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' inversions.bed > breakpoints/margin70/inversions.breakpoints.bed

awk -v m=70 '
/DEL/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' deletions.bed > breakpoints/margin70/deletions.breakpoints.bed

awk -v m=70 '
/DUP/ {printf  "%s %s %s %s_l\n%s %s %s %s_r\n", $1, $2-m, $2+m, $4, $1, $3-m, $3+m, $4}
' duplications.bed > breakpoints/margin70/duplications.breakpoints.bed