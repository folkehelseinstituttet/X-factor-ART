# do some bootstrapping to check the consistency of significant CpGs 
# with 30% of the ART samples and the control samples
# Author: Haakon Nustad

library(Rfast)
library(stringr)
library(bacon)


# set wd and load data:
pheno = readRDS("data/covar_mat_used.rds")
data = readRDS("data/DNAm_matrix_used.rds")

dim(data)
data[1:10,1:4]
dim(pheno)

# load also phenotype, family relationships

plate = key_child[rownames(pheno),"AMP_Plate"]
plate = setNames(as.integer(as.factor(plate)), rownames(pheno))
identical(names(plate), rownames(pheno))

# --- Model setup ----
# function specifications:
lme_fast <- function(x, V, id, seq = 2){nx <- complete.cases(x,V);return(do.call(cbind,rint.reg(x[nx], V[nx,], as.factor(id[nx]))[2:3])[seq,])}
get_var_name <- function(ff)str_replace_all(unlist(strsplit(as.character(ff)[3], "[+]")), " ", "")


# main model and results:
table(pheno[,"isMale"])
phenoDF = data.frame(pheno)
table(phenoDF[phenoDF$ART == 1,]$isMale)
#males
boys = rownames(phenoDF[phenoDF$isMale == 1,])
art_boys = phenoDF[phenoDF$isMale == 1,]$ART
mMat_boys = log2(data[,boys]/(1 - data[,boys]))
covar.mat.use_boys = pheno[boys,]
plate_boys = plate[rownames(covar.mat.use_boys)]

ggplot(data = data.frame(betas = as.vector(data[rownames(trueBoys),boys]), x = rep(rownames(trueBoys), times = dim(data[,boys])[2]), art = factor(rep(art_boys, each = 2))), aes(y = betas, x = x, col = art)) + geom_boxplot()

identical(colnames(mMat_boys), rownames(covar.mat.use_boys))
identical(colnames(mMat_boys), names(plate_boys))

ff.1c_boys <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
table(na.omit(covar.mat.use_boys[,get_var_name(ff.1c_boys)])[,"ART"])
res.1c_boys <- t(apply(mMat_boys,1, lme_fast, V = covar.mat.use_boys[,get_var_name(ff.1c_boys)], id = plate_boys)) 

z_stat = res.1c_boys[,1]/res.1c_boys[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1c_boys = cbind(res.1c_boys, z_stat, ps, ps_adj_BH)
sum(res.1c_boys[,5] < 0.05)
#trueBoys = rownames(res.1c_boys[res.1c_boys[,5] < 0.05,])

bc = bacon(res.1c_boys[,3])
estimates(bc)
inflation(bc)
bias(bc)
res.1c_boys = cbind(res.1c_boys, bacon_pval = pval(bc))
sum(p.adjust(res.1c_boys[,6], method = "BH") < 0.05)
rownames(res.1c_boys[p.adjust(res.1c_boys[,6], method = "BH") < 0.05,])

res.1c_boys = cbind(res.1c_boys, p.adjust(res.1c_boys[,6], method = "BH"))
trueBoys = res.1c_boys[which(res.1c_boys[,7] < 0.01),]
#female:
girls = rownames(phenoDF[phenoDF$isMale == 0,])
art_girls = phenoDF[phenoDF$isMale == 0,]$ART
mMat_girls = log2(data[,girls]/(1 - data[,girls]))
covar.mat.use_girls = pheno[girls,]
plate_girls = plate[rownames(covar.mat.use_girls)]

ff.1c_girls <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
table(na.omit(covar.mat.use_girls[,get_var_name(ff.1c_girls)])[,"ART"])
res.1c_girls <- t(apply(mMat_girls,1, lme_fast, V = covar.mat.use_girls[,get_var_name(ff.1c_girls)], id = plate_girls)) 

z_stat = res.1c_girls[,1]/res.1c_girls[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1c_girls = cbind(res.1c_girls, z_stat, ps, ps_adj_BH)
sum(res.1c_girls[,5] < 0.01)

bc = bacon(res.1c_girls[,3])
estimates(bc)
inflation(bc)
bias(bc)
res.1c_girls = cbind(res.1c_girls, bacon_pval = pval(bc))
sum(p.adjust(res.1c_girls[,6], method = "BH") < 0.01)
res.1c_girls = cbind(res.1c_girls, p.adjust(res.1c_girls[,6], method = "BH"))
trueGirls = res.1c_girls[which(res.1c_girls[,7] < 0.01),]
fit(bc, n = 100)

ggplot(data = data.frame(betas = as.vector(data[rownames(trueGirls),girls]), x = rep(rownames(trueGirls), times = dim(data[,girls])[2]), art = factor(rep(art_girls, each = 3))), aes(y = betas, x = x, col = art)) + geom_boxplot()
ggplot(data = data.frame(betas = as.vector(data[rownames(trueGirls)[3],girls]), x = rep(rownames(trueGirls)[3], times = dim(data[,girls])[2]), art = factor(rep(art_girls, each = 1))), aes(y = betas, x = x, col = art)) + geom_boxplot()


# bootstrapping::
# ----------------------------------
# males:
n = 1000 #bootstraps with replacement
prop = 1 # how large prop of the original size
sART = round(sum(pheno[, "isMale"] == 1 & pheno[, "ART"] == 1)*prop)
s = round(sum(pheno[, "isMale"] == 1 & pheno[, "ART"] == 0)*prop)
CpGListBoys = list()
t = Sys.time()
for(i in 1:n){
  bootC = sample(rownames(phenoDF[phenoDF$isMale == 1 & phenoDF$ART == 0,]), size = s, replace = T)
  bootA = sample(rownames(phenoDF[phenoDF$isMale == 1 & phenoDF$ART == 1,]), size = sART, replace = T)
  bootS = c(bootC, bootA)
  
  mMat_boys = log2(data[,bootS]/(1 - data[,bootS]))
  covar.mat.use_boys = pheno[bootS,]
  plate_boys = plate[rownames(covar.mat.use_boys)]
  
  identical(colnames(mMat_boys), rownames(covar.mat.use_boys))
  identical(colnames(mMat_boys), names(plate_boys))
  
  ff.1c_boys <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
  table(na.omit(covar.mat.use_boys[,get_var_name(ff.1c_boys)])[,"ART"])
  res.1c_boys <- t(apply(mMat_boys,1, lme_fast, V = covar.mat.use_boys[,get_var_name(ff.1c_boys)], id = plate_boys)) 
  
  z_stat = res.1c_boys[,1]/res.1c_boys[,2]
  ps = pnorm(-abs(z_stat))*2
  ps_adj_BH = p.adjust(ps, method = "BH")
  res.1c_boys = cbind(res.1c_boys, z_stat, ps, ps_adj_BH)
  sum(res.1c_boys[,5] < 0.01)
  
  bc = bacon(res.1c_boys[,3])
  estimates(bc)
  inflation(bc)
  bias(bc)
  res.1c_boys = cbind(res.1c_boys, bacon_pval = pval(bc))
  sum(p.adjust(res.1c_boys[,6], method = "BH") < 0.01)
  CpGListBoys[[i]] = names(which(p.adjust(res.1c_boys[,6], method = "BH") < 0.01))
  
}
Sys.time() - t
sort(table(Reduce(c, CpGListBoys)), decreasing = T)[1:20]/n
intersect(trueBoys, names(sort(table(Reduce(c, CpGListBoys)), decreasing = T)[1:10]))


# females:
table(phenoDF[phenoDF$ART == 1,]$isMale)
# females:
n = 1000 # bootstraps with replacement
prop = 1 # how large prop of the original size
sART = round(sum(pheno[, "isMale"] == 0 & pheno[, "ART"] == 1)*prop)
s = round(sum(pheno[, "isMale"] == 0 & pheno[, "ART"] == 0)*prop)

CpGListGirls = list()
t = Sys.time()
for(i in 1:n){
  bootC = sample(rownames(phenoDF[phenoDF$isMale == 0 & phenoDF$ART == 0,]), size = s, replace = T)
  bootA = sample(rownames(phenoDF[phenoDF$isMale == 0 & phenoDF$ART == 1,]), size = sART, replace = T)
  bootS = c(bootC, bootA)
  
  mMat_girls = log2(data[,bootS]/(1 - data[,bootS]))
  covar.mat.use_girls = pheno[bootS,]
  plate_girls = plate[rownames(covar.mat.use_girls)]
  
  ff.1c_girls <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
  table(na.omit(covar.mat.use_girls[,get_var_name(ff.1c_girls)])[,"ART"])
  res.1c_girls <- t(apply(mMat_girls,1, lme_fast, V = covar.mat.use_girls[,get_var_name(ff.1c_girls)], id = plate_girls)) 
  
  z_stat = res.1c_girls[,1]/res.1c_girls[,2]
  ps = pnorm(-abs(z_stat))*2
  ps_adj_BH = p.adjust(ps, method = "BH")
  res.1c_girls = cbind(res.1c_girls, z_stat, ps, ps_adj_BH)
  sum(res.1c_girls[,5] < 0.01)
  
  bc = bacon(res.1c_girls[,3])
  estimates(bc)
  inflation(bc)
  bias(bc)
  res.1c_girls = cbind(res.1c_girls, bacon_pval = pval(bc))
  sum(p.adjust(res.1c_girls[,6], method = "BH") < 0.01)
  CpGListGirls[[i]] = names(which(p.adjust(res.1c_girls[,6], method = "BH") < 0.01))
  
}
Sys.time() - t
sort(table(Reduce(c, CpGListGirls)), decreasing = T)[1:40]/n


intersect(rownames(trueGirls), names(sort(table(Reduce(c, CpGListGirls)), decreasing = T)[1:40]))


# check with previous results:
#load("data/XchromosomeResults.RData")

which(Boys.1c[[2]][,5] < 0.01)
trueBoys
plot(Boys.1c[[2]][,3], res.1c_boys[,3])

trueGirls
which(p.adjust(Girls.1c[[2]][,6], method = "BH") < 0.01)
plot(Girls.1c[[2]][,3], res.1c_girls[,3])

save(trueBoys, trueGirls, CpGListBoys, CpGListGirls, file = "data/bootstrapp_res100percent.RData")


# load data and check:
load("data/bootstrapp_res100percent.RData")
girls_res = sort(table(Reduce(c, CpGListGirls)), decreasing = T)/1000
boys_res = sort(table(Reduce(c, CpGListBoys)), decreasing = T)/1000
girls_res2 = data.frame(CpG = names(girls_res), proportion = as.numeric(girls_res))
boys_res2 = data.frame(CpG = names(boys_res), proportion = as.numeric(boys_res))

write.csv(girls_res2, file = "data/bootstrap_results_girls_final.csv", row.names = F)
write.csv(boys_res2, file = "data/bootstrap_results_boys_final.csv", row.names = F)



