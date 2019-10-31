library(rehh)
data("haplohh_cgu_bta12")
class(haplohh_cgu_bta12@haplo)
slotNames(haplohh_cgu_bta12)


dim(haplohh_cgu_bta12@haplo)

## these are also kept in nhap and nsnp

haplohh_cgu_bta12@haplo[1:10,1:5]
## kept as 0s and 1s

data("haplohh_cgu_bta12")
#layout(matrix(1:2, 2, 1))

bifurcation.diagram(haplohh_cgu_bta12, mrk_foc=456, all_foc=1, nmrk_l=20, nmrk_r=20)
source("bifurcation_diagram.R")
bifurcation.diagram2(haplohh_cgu_bta12, mrk_foc=456, all_foc=1, nmrk_l=20, nmrk_r=20)
