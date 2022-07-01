# Analysis on X chromosome
# Author: Haakon Nustad
# Code adapted by Julia Romanowska

#---- SETUP ----
library(here)
library(stringr)
library(dplyr)
library(plyr)
library(Rfast)
library(venn)
library(bacon)
library(ggplot2)

# ---- Load data ----
# This part is specific for MoBa directories at TSD, so won't be made public
# base_dir - start directory
# dnam_dir - directory with QCed DNAm beta values

# loading DNAm data - bmat.child.chrX.BMIQ.Random.noTed, bmat.child.chrX.BMIQ.START.noTed (and for parents)
# loading probe information - annot.epic.b5

# loading questionnaire data and phenotype data - MoBa_Q1_Q4_MBRN

# --- Prepare data ----
chrX_probes = annot.epic.b5[annot.epic.b5$CHR == "X",]
rm(annot.epic.b5)

xCpGsSTART = rownames(bmat.child.chrX.BMIQ.Random.noTed)
chrX_probes = chrX_probes[intersect(rownames(chrX_probes), xCpGsSTART),]
table(chrX_probes$Infinium_Design_Type)

# fix meth data samples:
child_rand = colnames(bmat.child.chrX.BMIQ.Random.noTed)
child_start = colnames(bmat.child.chrX.BMIQ.START.noTed)
father_rand = colnames(bmat.father.chrX.BMIQ.Random.noTed)
father_start = colnames(bmat.father.chrX.BMIQ.START.noTed)
mother_rand = colnames(bmat.mother.chrX.BMIQ.Random.noTed)
mother_start = colnames(bmat.mother.chrX.BMIQ.START.noTed)

art_child = child_start
art_child = c(art_child, child_rand[as.logical(start_art[child_rand,])])
control_child = child_rand[!as.logical(start_art[child_rand,])]

# bmat gathering:
bmat = cbind(bmat.child.chrX.BMIQ.Random.noTed, bmat.child.chrX.BMIQ.START.noTed)
bmat = bmat[,c(control_child, art_child)]


## Remove unrealistic BMI measurments
MoBa_Q1_Q4_MBRN$AA87[which(MoBa_Q1_Q4_MBRN$AA87 == 16)] <- NA 
MoBa_Q1_Q4_MBRN$AA87[which(MoBa_Q1_Q4_MBRN$AA87 == 1)] <- NA
MoBa_Q1_Q4_MBRN$mat_bmi_pre <- MoBa_Q1_Q4_MBRN$AA85/((MoBa_Q1_Q4_MBRN$AA87/100)^2)

covariates <- data.frame(
  ART		= MoBa_Q1_Q4_MBRN$ANY_ART,
  TTP_planned	= MoBa_Q1_Q4_MBRN$ttp_planned,
  mat_age	= MoBa_Q1_Q4_MBRN$MORS_ALDER,
  pat_age	= MoBa_Q1_Q4_MBRN$FARS_ALDER, 
  isPrimiparous  = ifelse(MoBa_Q1_Q4_MBRN$PARITET_5 == "0 (primiparous)", 0,1),
  isMale 	= ifelse(MoBa_Q1_Q4_MBRN$kjonn_epic == "Male", 1, 0), # "Epigenetic" sex when discordant
  mat_edu    	= factor(MoBa_Q1_Q4_MBRN$mat_edu, levels = c(3,1,2,4)),
  mat_bmi    	= MoBa_Q1_Q4_MBRN$mat_bmi_pre,
  mat_smk    	= factor(MoBa_Q1_Q4_MBRN$daily_before_during_preg),
  GA         	= MoBa_Q1_Q4_MBRN$SVLEN_START,
  BW		= MoBa_Q1_Q4_MBRN$VEKT, 
  TTP_8     	= ifelse(MoBa_Q1_Q4_MBRN$ttp_planned > 8, 1,0), 
  TTP_12     	= ifelse(MoBa_Q1_Q4_MBRN$ttp_planned > 12, 1, 0), 
  isFresh	= ifelse(MoBa_Q1_Q4_MBRN$art_iiff %in% c("fresh icsi", "fresh ivf"), 1, 0),
  isFrozen  	= ifelse(MoBa_Q1_Q4_MBRN$art_iiff %in% c("frozen icsi", "frozen ivf"), 1, 0),
  plate     	= key_child[rownames(MoBa_Q1_Q4_MBRN),"AMP_Plate"],
  row.names 	= rownames(MoBa_Q1_Q4_MBRN)
)

use <- Reduce(intersect, list(rownames(covariates), colnames(bmat)))

covar.use   <- covariates[use,]
mMat        <- log2(bmat[,use]/(1-bmat[,use]))
rm(bmat, covariates)

plate	    <- setNames(as.integer(as.factor(covar.use$plate)), use)

# change types and colnames of covariates:
covar.mat.use <- as.matrix(
  cbind(
    covar.use[,c("ART", "mat_age","pat_age","mat_bmi")], 
    ifelse(covar.use$mat_smk == 1, 1, 0),
    ifelse(covar.use$mat_smk == 2, 1, 0),
    ifelse(covar.use$mat_smk == 3, 1, 0),
    ifelse(covar.use$mat_edu == 1, 1, 0),
    ifelse(covar.use$mat_edu == 2, 1, 0),
    ifelse(covar.use$mat_edu == 4, 1, 0),
    covar.use[,c("isPrimiparous", "isMale","GA","BW","TTP_8", "TTP_12", "isFresh", "isFrozen")]
  )
)

colnames(covar.mat.use) <-  c(
  "ART", "mat_age", "pat_age", "mat_bmi", "smk1", "smk2", "smk3",
  "edu1", "edu2", "edu4", "isPrimiparous", "isMale", "GA","BW","TTP_8","TTP_12",
  "isFresh", "isFrozen"
)

# --- Model setup ----
# function specifications:
lme_fast <- function(x, V, id, seq = 2){
  nx <- complete.cases(x,V)
  return(
    do.call(
      cbind, 
      rint.reg(x[nx], V[nx,], as.factor(id[nx]))[2:3])[seq,]
    )
}
get_var_name <- function(ff){
  str_replace_all(unlist(strsplit(as.character(ff)[3], "[+]")), " ", "")
}

# BW and GA model ART paper: 
ff.1d.bw <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + isMale + GA + BW)
table(na.omit(covar.mat.use[,get_var_name(ff.1d.bw)])[, "ART"])
res.1d.bw <- t(apply(
  mMat, 1, lme_fast, V = covar.mat.use[,get_var_name(ff.1d.bw)], id = plate
  )
) 

z_stat = res.1d.bw[,1]/res.1d.bw[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1d.bw = cbind(res.1d.bw, z_stat, ps, ps_adj_BH)

sum(res.1d.bw[, 5] < 0.01)
comb.1d.bw = rownames(res.1d.bw[res.1d.bw[, 5] < 0.01, ])

# main model ART paper:
ff.1c <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + isMale)
table(na.omit(covar.mat.use[,get_var_name(ff.1c)])[,"ART"])
res.1c <- t(apply(mMat,1, lme_fast, V = covar.mat.use[,get_var_name(ff.1c)], id = plate)) 

z_stat = res.1c[,1]/res.1c[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1c = cbind(res.1c, z_stat, ps, ps_adj_BH)

sum(res.1c[, 5] < 0.01)
comb.1c = rownames(res.1c[res.1c[, 5] < 0.01, ])

# overlap:
res.1c[rownames(res.1c[res.1c[, 5] < 0.01, ]),]
res.1d.bw[rownames(res.1d.bw[res.1d.bw[, 5] < 0.01, ]),]



# --- divide on sex, and check results:----
  # males:
mMat_boys = mMat[,covar.mat.use[,"isMale"] == 1]
covar.mat.use_boys = covar.mat.use[covar.mat.use[,"isMale"] == 1,]
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
boys.1c = rownames(res.1c_boys[res.1c_boys[,5] < 0.01,])

ff.1d.bw_boys <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + GA + BW)
table(na.omit(covar.mat.use_boys[,get_var_name(ff.1d.bw_boys)])[,"ART"])
res.1d.bw_boys <- t(apply(mMat_boys,1, lme_fast, V = covar.mat.use_boys[,get_var_name(ff.1d.bw_boys)], id = plate_boys)) 

z_stat = res.1d.bw_boys[,1]/res.1d.bw_boys[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1d.bw_boys = cbind(res.1d.bw_boys, z_stat, ps, ps_adj_BH)
sum(res.1d.bw_boys[,5] < 0.01)
boys.1d.bw = rownames(res.1d.bw_boys[res.1d.bw_boys[,5] < 0.01,])


  # Females:
mMat_girls = mMat[,covar.mat.use[,"isMale"] == 0]
covar.mat.use_girls = covar.mat.use[covar.mat.use[,"isMale"] == 0,]
plate_girls = plate[rownames(covar.mat.use_girls)]

identical(colnames(mMat_girls), rownames(covar.mat.use_girls))
identical(colnames(mMat_girls), names(plate_girls))

ff.1c_girls <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
table(na.omit(covar.mat.use_girls[,get_var_name(ff.1c_girls)])[,"ART"])
res.1c_girls <- t(apply(mMat_girls,1, lme_fast, V = covar.mat.use_girls[,get_var_name(ff.1c_girls)], id = plate_girls)) 

z_stat = res.1c_girls[,1]/res.1c_girls[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1c_girls = cbind(res.1c_girls, z_stat, ps, ps_adj_BH)
sum(res.1c_girls[,5] < 0.01)
girls.1c = rownames(res.1c_girls[res.1c_girls[,5] < 0.01,])

ff.1d.bw_girls <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + GA + BW)
table(na.omit(covar.mat.use_girls[,get_var_name(ff.1d.bw_girls)])[,"ART"])
res.1d.bw_girls <- t(apply(mMat_girls,1, lme_fast, V = covar.mat.use_girls[,get_var_name(ff.1d.bw_girls)], id = plate_girls)) 

z_stat = res.1d.bw_girls[,1]/res.1d.bw_girls[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1d.bw_girls = cbind(res.1d.bw_girls, z_stat, ps, ps_adj_BH)
sum(res.1d.bw_girls[,5] < 0.01)
girls.1d.bw = rownames(res.1d.bw_girls[res.1d.bw_girls[,5] < 0.01,])

# ------------------------------------------------------------------------------------------


# controling for parent methylation: ----
bmat.father <- cbind(bmat.father.chrX.BMIQ.Random.noTed, bmat.father.chrX.BMIQ.START.noTed)
rm(bmat.father.chrX.BMIQ.Random.noTed, bmat.father.chrX.BMIQ.START.noTed)

bmat.mother <- cbind(bmat.mother.chrX.BMIQ.Random.noTed, bmat.mother.chrX.BMIQ.START.noTed)
rm(bmat.mother.chrX.BMIQ.Random.noTed, bmat.mother.chrX.BMIQ.START.noTed)

length(intersect(intersect(colnames(bmat.mother), use), colnames(bmat.father)))

mMat.mother <- log2(bmat.mother/(1-bmat.mother))
mMat.father <- log2(bmat.father/(1-bmat.father))
use_mp = intersect(intersect(use, colnames(mMat.father)), colnames(mMat.mother))

mMat.mother <- mMat.mother[, intersect(use_mp, colnames(mMat.mother))]
mMat.father <- mMat.father[, intersect(use_mp, colnames(mMat.father))]

mMat_kids = mMat[,intersect(use_mp, colnames(mMat))]

cg  <- Reduce(intersect, list(rownames(mMat_kids),rownames(mMat.mother), rownames(mMat.father)))
ind <- Reduce(intersect, list(colnames(mMat_kids),colnames(mMat.mother), colnames(mMat.father)))


ff.1fc  <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + isMale)
res.1fc <- t(sapply(cg, function(x,mMat,mMat_2, mMat_3,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,], mMat_3[x,]), id = id),
                    mMat = mMat_kids[,ind], mMat_2 = mMat.mother[,ind], mMat_3 = mMat.father[,ind],
                    cov = covar.mat.use[ind,get_var_name(ff.1fc)], id = plate[ind]))

z_stat = res.1fc[,1]/res.1fc[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc = cbind(res.1fc, z_stat, ps, ps_adj_BH)
sum(res.1fc[,5] < 0.01)

comb.1c.parent = rownames(res.1fc[res.1fc[,5] < 0.01,])

ff.1fc_bw  <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + isMale + GA + BW)
res.1fc_bw <- t(sapply(cg, function(x,mMat,mMat_2, mMat_3,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,], mMat_3[x,]), id = id),
                    mMat = mMat_kids[,ind], mMat_2 = mMat.mother[,ind], mMat_3 = mMat.father[,ind],
                    cov = covar.mat.use[ind,get_var_name(ff.1fc_bw)], id = plate[ind]))

z_stat = res.1fc_bw[,1]/res.1fc_bw[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc_bw = cbind(res.1fc_bw, z_stat, ps, ps_adj_BH)
sum(res.1fc_bw[,5] < 0.01)
comb.1d.bw.parent = rownames(res.1fc_bw[res.1fc_bw[,5] < 0.01,])


  # Boys:
ind_boys <- Reduce(intersect, list(colnames(mMat_boys),colnames(mMat.mother), colnames(mMat.father)))

ff.1fc_boys  <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
res.1fc_boys <- t(sapply(cg, function(x,mMat,mMat_2,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,]), id = id),
                    mMat = mMat_boys[,ind_boys], mMat_2 = mMat.mother[,ind_boys],
                    cov = covar.mat.use[ind_boys,get_var_name(ff.1fc_boys)], id = plate[ind_boys]))

table(na.omit(covar.mat.use[ind_boys,get_var_name(ff.1fc_boys)])[,"ART"])

z_stat = res.1fc_boys[,1]/res.1fc_boys[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc_boys = cbind(res.1fc_boys, z_stat, ps, ps_adj_BH)
sum(res.1fc_boys[,5] < 0.01)

boys.1c.parent = rownames(res.1fc_boys[res.1fc_boys[,5] < 0.01,])

ff.1fc_bw_boys  <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + GA + BW)
res.1fc_bw_boys <- t(sapply(cg, function(x,mMat,mMat_2,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,]), id = id),
                       mMat = mMat_boys[,ind_boys], mMat_2 = mMat.mother[,ind_boys],
                       cov = covar.mat.use[ind_boys,get_var_name(ff.1fc_bw_boys)], id = plate[ind_boys]))

table(na.omit(covar.mat.use[ind_boys,get_var_name(ff.1fc_bw_boys)])[,"ART"])

z_stat = res.1fc_bw_boys[,1]/res.1fc_bw_boys[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc_bw_boys = cbind(res.1fc_bw_boys, z_stat, ps, ps_adj_BH)
sum(res.1fc_bw_boys[,5] < 0.01)
boys.1d.bw.parent = rownames(res.1fc_bw_boys[res.1fc_bw_boys[,5] < 0.01,])


  # Girls:
ind_girls <- Reduce(intersect, list(colnames(mMat_girls),colnames(mMat.mother), colnames(mMat.father)))

ff.1fc_girls  <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous)
res.1fc_girls <- t(sapply(cg, function(x,mMat,mMat_2, mMat_3,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,], mMat_3[x,]), id = id),
                         mMat = mMat_girls[,ind_girls], mMat_2 = mMat.mother[,ind_girls], mMat_3 = mMat.father[,ind_girls],
                         cov = covar.mat.use[ind_girls,get_var_name(ff.1fc_girls)], id = plate[ind_girls]))

table(na.omit(covar.mat.use[ind_girls,get_var_name(ff.1fc_girls)])[,"ART"])

z_stat = res.1fc_girls[,1]/res.1fc_girls[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc_girls = cbind(res.1fc_girls, z_stat, ps, ps_adj_BH)
sum(res.1fc_girls[,5] < 0.01)
girls.1c.parent = rownames(res.1fc_girls[res.1fc_girls[,5] < 0.01,])


ff.1fc_bw_girls  <- as.formula(cg ~ ART + mat_age + smk1 + smk2 + smk3 + mat_bmi + isPrimiparous + GA + BW)
res.1fc_bw_girls <- t(sapply(cg, function(x,mMat,mMat_2, mMat_3,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,], mMat_3[x,]), id = id),
                            mMat = mMat_girls[,ind_girls], mMat_2 = mMat.mother[,ind_girls], mMat_3 = mMat.father[,ind_girls],
                            cov = covar.mat.use[ind_girls,get_var_name(ff.1fc_bw_girls)], id = plate[ind_girls]))

table(na.omit(covar.mat.use[ind_girls,get_var_name(ff.1fc_bw_girls)])[,"ART"])

z_stat = res.1fc_bw_girls[,1]/res.1fc_bw_girls[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc_bw_girls = cbind(res.1fc_bw_girls, z_stat, ps, ps_adj_BH)
sum(res.1fc_bw_girls[,5] < 0.01)
girls.1d.bw.parent = rownames(res.1fc_bw_girls[res.1fc_bw_girls[,5] < 0.01,])

# Lets do bacon: ----
# girls:
bc = bacon(res.1c_girls[,3])
estimates(bc)
inflation(bc)
bias(bc)
res.1c_girls = cbind(res.1c_girls, pval(bc))
sum(p.adjust(res.1c_girls[,6], method = "BH") < 0.01)
which(p.adjust(res.1c_girls[,6], method = "BH") < 0.01)
bc2 = bacon(res.1d.bw_girls[,3])
res.1d.bw_girls = cbind(res.1d.bw_girls, pval(bc2))
sum(p.adjust(res.1d.bw_girls[,6], method = "BH") < 0.01)
which(p.adjust(res.1d.bw_girls[,6], method = "BH") < 0.01)

bc = bacon(res.1fc_girls[,3])
estimates(bc)
inflation(bc)
bias(bc)
res.1fc_girls = cbind(res.1fc_girls, pval(bc))
sum(p.adjust(res.1fc_girls[,6], method = "BH") < 0.01)
which(p.adjust(res.1fc_girls[,6], method = "BH") < 0.01)
bc2 = bacon(res.1fc_bw_girls[,3])
res.1fc_bw_girls = cbind(res.1fc_bw_girls, pval(bc2))
sum(p.adjust(res.1fc_bw_girls[,6], method = "BH") < 0.01)
which(p.adjust(res.1fc_bw_girls[,6], method = "BH") < 0.01)


# boys, should not alter anything:
bc = bacon(res.1c_boys[,3])
estimates(bc)
inflation(bc)
bias(bc)
res.1c_boys = cbind(res.1c_boys, pval(bc))
sum(p.adjust(res.1c_boys[,6], method = "BH") < 0.01)
which(p.adjust(res.1c_boys[,6], method = "BH") < 0.01)
bc2 = bacon(res.1d.bw_boys[,3])
res.1d.bw_boys = cbind(res.1d.bw_boys, pval(bc2))
sum(p.adjust(res.1d.bw_boys[,6], method = "BH") < 0.01)
which(p.adjust(res.1d.bw_boys[,6], method = "BH") < 0.01)

bc = bacon(res.1fc_boys[,3])
estimates(bc)
inflation(bc)
bias(bc)
res.1fc_boys = cbind(res.1fc_boys, pval(bc))
sum(p.adjust(res.1fc_boys[,6], method = "BH") < 0.01)
which(p.adjust(res.1fc_boys[,6], method = "BH") < 0.01)
bc2 = bacon(res.1fc_bw_boys[,3])
res.1fc_bw_boys = cbind(res.1fc_bw_boys, pval(bc2))
sum(p.adjust(res.1fc_bw_boys[,6], method = "BH") < 0.01)
which(p.adjust(res.1fc_bw_boys[,6], method = "BH") < 0.01)


#lets just try something:
plateAndSentrixPos = setNames(as.integer(as.factor(paste0(key_child$AMP_Plate, key_child$SentrixPosition_A))), rownames(key_child))
res.1fc_bw_girls_adjusted <- t(sapply(cg, function(x,mMat,mMat_2, mMat_3,cov, id)lme_fast(mMat[x,],V = cbind(cov,mMat_2[x,], mMat_3[x,]), id = id),
                             mMat = mMat_girls[,ind_girls], mMat_2 = mMat.mother[,ind_girls], mMat_3 = mMat.father[,ind_girls],
                             cov = covar.mat.use[ind_girls,get_var_name(ff.1fc_bw_girls)], id = plateAndSentrixPos[ind_girls]))
z_stat = res.1fc_bw_girls_adjusted[,1]/res.1fc_bw_girls_adjusted[,2]
ps = pnorm(-abs(z_stat))*2
ps_adj_BH = p.adjust(ps, method = "BH")
res.1fc_bw_girls_adjusted = cbind(res.1fc_bw_girls_adjusted, z_stat, ps, ps_adj_BH)


# Save data ----
# save the data in RData object:
Both.1c = list(ff.1c, res.1c)
Both.1d.bw = list(ff.1d.bw, res.1d.bw)
Boys.1c = list(ff.1c_boys, res.1c_boys)
Boys.1d.bw = list(ff.1d.bw_boys, res.1d.bw_boys)
Girls.1c = list(ff.1c_girls, res.1c_girls)
Girls.1d.bw = list(ff.1d.bw_girls, res.1d.bw_girls)
Both.parents.1c = list(ff.1fc, res.1fc)
Both.parents.1d.bw = list(ff.1fc_bw, res.1fc_bw)
Boys.parents.1c = list(ff.1fc_boys, res.1fc_boys)
Boys.parents.1d.bw = list(ff.1fc_bw_boys, res.1fc_bw_boys)
Girls.parents.1c = list(ff.1fc_girls, res.1fc_girls)
Girls.parents.1d.bw = list(ff.1fc_bw_girls, res.1fc_bw_girls)
save(Both.1c, Both.1d.bw, Boys.1c, Boys.1d.bw, Girls.1c, Girls.1d.bw, 
     Both.parents.1c, Both.parents.1d.bw, Boys.parents.1c, Boys.parents.1d.bw, 
     Girls.parents.1c, Girls.parents.1d.bw,
     file = file.path(here::here("data"), "XchromosomeResults.RData"))

# check intersecting CpGs:
# all intersecting for comb and boys:
head(chrX_probes)
chrX_probes[comb.1c,c(11, 12, 49, 15, 19, 43, 44)]
chrX_probes[boys.1c,c(11, 12, 49, 15, 19, 43, 44)]

girlOverlap = Reduce(intersect, list(girls.1c, girls.1c.parent, girls.1d.bw, girls.1d.bw.parent))
chrX_probes[girlOverlap,c(11, 12, 49, 15)]
chrX_probes[girlOverlap,c(19, 43, 44)]

# Plotting: ----
# Is venndiagram enough?
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#tiff(filename ="Fig_S2bw.tiff", width = 18, height = 15, units = "cm", res = 300, compression = "lzw")
vennModel1a4a = venn(x = list("1" = comb.1c,
                              "2" = comb.1c.parent,
                              "3" = comb.1d.bw, 
                              "4" = comb.1d.bw.parent), zcolor = cbPalette[1:4], ilcs = 1.3, sncs = 1.2, box = FALSE)
#dev.off()

vennModel1a4a = venn(x = list("1" = boys.1c,
                              "2" = boys.1c.parent,
                              "3" = boys.1d.bw, 
                              "4" = boys.1d.bw.parent), zcolor = cbPalette[1:4], ilcs = 1.3, sncs = 1.2, box = FALSE)


vennModel1a4a = venn(x = list("1" = girls.1c,
                              "2" = girls.1c.parent,
                              "3" = girls.1d.bw, 
                              "4" = girls.1d.bw.parent), zcolor = cbPalette[1:4], ilcs = 1.3, sncs = 1.2, box = FALSE)

ggplot(data = data.frame(y = -log10(c(sort(res.1c[,4], decreasing = T), 
                         sort(res.1fc[,4], decreasing = T), 
                         sort(res.1d.bw[,4], decreasing = T),
                         sort(res.1fc_bw[,4], decreasing = T))), model = factor(rep(c("1", "2", "3", "4"), each = dim(res.1c)[1])), 
                         x = -log10(rep(sort(1:dim(res.1c)[1]/(dim(res.1c)[1]), decreasing = T), times = 4))), aes(y = y, col = model, x = x)) + geom_point(alpha = 0.7) + geom_abline(intercept = 0, slope = 1)

#boys
ggplot(data = data.frame(y = -log10(c(sort(res.1c_boys[,6], decreasing = T), 
                                      sort(res.1fc_boys[,6], decreasing = T), 
                                      sort(res.1d.bw_boys[,6], decreasing = T),
                                      sort(res.1fc_bw_boys[,6], decreasing = T))), model = factor(rep(c("1", "2", "3", "4"), each = dim(res.1c_boys)[1])), 
                         x = -log10(rep(sort(1:dim(res.1c_boys)[1]/(dim(res.1c_boys)[1]), decreasing = T), times = 4))), aes(y = y, col = model, x = x)) + geom_point(alpha = 0.7) + geom_abline(intercept = 0, slope = 1)
#Girls
ggplot(data = data.frame(y = -log10(c(sort(res.1c_girls[,4], decreasing = T), 
                                      sort(res.1fc_girls[,4], decreasing = T), 
                                      sort(res.1d.bw_girls[,4], decreasing = T),
                                      sort(res.1fc_bw_girls[,4], decreasing = T))), model = factor(rep(c("1", "2", "3", "4"), each = dim(res.1c_girls)[1])), 
                         x = -log10(rep(sort((1:dim(res.1c_girls)[1]/(dim(res.1c_girls)[1])), decreasing = T), times = 4))), aes(y = y, col = model, x = x)) + geom_point(alpha = 0.7) + geom_abline(intercept = 0, slope = 1)



plot(x=-log10(rep(sort(1:dim(res.1c)[1]/(dim(res.1c)[1]), decreasing = T), times = 4)),y =-log10(rep(sort(pnorm(rnorm(n = dim(res.1c)[1])), decreasing = T), times = 4)))


ggplot(data = data.frame(y = -log10(c(sort(res.1c_girls[,4], decreasing = T), 
                                      sort(res.1fc_girls[,4], decreasing = T), 
                                      sort(res.1d.bw_girls[,4], decreasing = T),
                                      sort(res.1fc_bw_girls[,4], decreasing = T), 
                                      sort(res.1fc_bw_girls_adjusted[,4], decreasing = T))), model = factor(rep(c("1", "2", "3", "4", "5"), each = dim(res.1c)[1])), 
                         x = -log10(rep(sort((1:dim(res.1c)[1]/(dim(res.1c)[1])), decreasing = T), times = 5))), aes(y = y, col = model, x = x)) + geom_point(alpha = 0.7) + geom_abline(intercept = 0, slope = 1)


ggplot(data = data.frame(x = res.1c_boys[,1]), aes(x = x)) + geom_histogram(bins = 50) + geom_vline(xintercept = mean(res.1c_boys[,1]))
ggplot(data = data.frame(x = res.1c_girls[,1]), aes(x = x)) + geom_histogram(bins = 50)+ geom_vline(xintercept = mean(res.1c_girls[,1]))
ggplot(data = data.frame(x = res.1c[,1]), aes(x = x)) + geom_histogram(bins = 50)+ geom_vline(xintercept = mean(res.1c[,1]))

ggplot(data = data.frame(boys = res.1c_boys[,3], girls = res.1c_girls[,3]), aes(x = boys, y = girls)) + geom_point(alpha = 0.5) #+ ylim(c(-0.6, 0.6)) + xlim(c(-0.6, 0.6))


gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 19, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

gg_qqplot(res.1fc_bw_girls[,4])


