data(wgscan.cgu)
data("wgscan.eut")
res.rsb <- ies2rsb(wgscan.cgu, wgscan.eut, "CGU", "EUT")
rsbplot(res.rsb)

res.xpehh <- ies2xpehh(wgscan.cgu, wgscan.eut, "CGU", "EUT")
xpehhplot(res.xpehh)

## Bifurcation diagrams

data("haplohh_cgu_bta12")
layout(matrix(1:2, 2, 1))
bifurcation.diagram(haplohh_cgu_bta12, mrk_foc=456, all_foc=1, nmrk_l=20, nmrk_r=20)
bifurcation.diagram(haplohh_cgu_bta12,  mrk_foc=456, all_foc=2, nmrk_l=20, nmrk_r=20) 
