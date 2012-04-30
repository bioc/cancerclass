setClass("nvalidation", representation(
ngenes="numeric", 
method="character", 
dist="character", 
ntrain="numeric", 
nrep="numeric", 
hparam="numeric", 
misclass="list", 
nselected="matrix", 
samples="matrix", 
classifier="character", 
fdata="data.frame" ) )
