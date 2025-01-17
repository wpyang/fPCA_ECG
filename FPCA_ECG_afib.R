
require(tidyverse)
require(fda)
require(StepReg)
require(DescTools)
require(gtsummary)
require(survival)
require(mediation)
library(ggpubr)



#####################################################################
### Step 1: FPCA
#####################################################################
npcs <- 8

ecg.fpca.summary2 <- list()
ecg.fpca.scores2 <- list()

### fd parameters
time <- 0:400/400
ntime <- length(time)
norder <- 6
nbasis <- ntime + norder - 2
rng <- range(time)
wbasis <- create.bspline.basis(rangeval=rng, nbasis=nbasis, norder=norder, breaks=time)
fd0 <- fd(matrix(0,nbasis,1), wbasis)
Wnbasis <- 4
Wbasis <- create.bspline.basis(rng, Wnbasis)
Wfd0 <- fd(matrix(0,Wnbasis,1),Wbasis)
fdPar <- fdPar(fd0, Lfdobj=2, lambda=1e-5)

for (i in 1:length(ecgvar)){

    lead <- ecgvar[i]
    ecg.pca.tmp <- pca.fd(ecgfd_cs_train[[lead]]$ecgfd,nharm=npcs,fdPar)
    ecg.varmx.pca.tmp <- varmx.pca.fd(ecg.pca.tmp)
  
    op <- par(mfrow=c(2,4))
    plot.pca.fd(ecg.varmx.pca.tmp)
    par(op)

    ecg.fpca.summary2[[lead]] <- ecg.varmx.pca.tmp    
    
    ## calculate fPC score using original data
    fpca.w <- eval.fd(time,ecg.varmx.pca.tmp$harmonics) 
    orig.dat <- ecgfd_cs_train[[lead]]$ecg.ori.sqrt
    orig.dat.demean <- orig.dat -rowMeans(orig.dat)
    pc.scores <- t(orig.dat.demean) %*% fpca.w/length(time)
    
    ecg.fpca.scores2[[lead]] <- ecg.varmx.pca.tmp$scores
    
    
}

ecg.fpca.scores2.df <- as.data.frame(do.call(cbind,ecg.fpca.scores2))
colnames(ecg.fpca.scores2.df) <- apply(expand.grid(1:npcs,ecgvar), 1, function(x) paste(x[2], x[1], sep="."))
ecg.fpca.scores2.df$SEQ <- as.numeric(ecgfd_cs_train$I$ecgfd$fdnames$reps)

ecg.fpca.scores2.afib <- merge(ecg.fpca.scores2.df,list_ecgs_crosssectional_training,by="SEQ")


#####################################################################
### Step 2: fPCA feature screening
#####################################################################

ttest.crosssectional <- function(data){
    
    datatmp <- data[with(data,order(PID)),]
    
    ## t-test
    ttest.summary1 <- NULL
    for (i in 1:length(ecgvar)){
        for (j in 1:npcs){

            ## test between groups
            tmp.result <- t.test(as.formula(paste0(ecgvar[i],".",j, "~afib")),data=datatmp)
            tmp.result <- as.data.frame(t(c(tmp.result$estimate,tmp.result$statistic,tmp.result$p.value)))
            names(tmp.result) <- c("mean0","mean1","tstat","pvalue")
                      
            tmp.result$lead <- ecgvar[i]
            tmp.result$pc <- j
            
            if (i ==1 & j==1 ){
                ttest.summary1 <-  tmp.result
            }
            else {
                ttest.summary1 <- rbind(ttest.summary1,tmp.result)
            }
        }
    }
    return(ttest.summary1)
}

ecg.fpca.varlist <- NULL
for (i in 1:length(ecgvar)){
  ecg.fpca.varlist <- c(ecg.fpca.varlist,paste0(ecgvar[i],".",c(1:npcs)))
}

ttest.summary2.afib <- ttest.crosssectional(ecg.fpca.scores2.afib)
ttest.summary2.afib$ecgvar <- paste0(ttest.summary2.afib$lead, ".", ttest.summary2.afib$pc)

ttest.summary2.afib.sig <- ttest.summary2.afib[ttest.summary2.afib$pvalue<=5e-4,]
ttest.summary2.afib.sig$ecgvar <- paste0(ttest.summary2.afib.sig$lead, ".", ttest.summary2.afib.sig$pc)

ecgvar.siglist <- paste0(ttest.summary2.afib.sig$lead, ".", ttest.summary2.afib.sig$pc)



##### stepwise selection  #####
model2 <- stepwise(as.formula(paste("afib", paste(ecgvar.siglist,collapse="+"), sep = "~")),
                        data=ecg.fpca.scores2.afib, 
                   type="logit",
                   strategy = "bidirection", metric="SL", sle=0.001)

model2$`Summary of coefficients for the selected model with afib under bidirection and SL `$Variable


################################################################
## Step 3: check PC across all 12 leads
################################################################

gen.w.func <- function(ecgvar, ecg.fpca.summary){
  lead.tmp <- gsub("\\..*","",ecgvar)
  pc.tmp <- as.numeric( gsub(".*\\.","",ecgvar) )
  
  ecg.fpca <- ecg.fpca.summary[[lead.tmp]]
  pc.w <- ecg.fpca$harmonics[pc.tmp]
  pc.w
}

gen.score.func <- function(lead, pc.w, ecgfd){
  t( inprod( pc.w, center.fd( ecgfd[[lead]]$ecgfd ) ) )
}



gen.plot.func <- function(lead.tmp, pc.w, ecgfd, data){
  ecgfd.tmp <- ecgfd[[lead.tmp]]$ecgfd
  pc.score.tmp <- data[,lead.tmp]
  
  ## fpca data
  tmp0 <- eval.fd(mean.fd(ecgfd.tmp),0:100/100)
  tmp1 <- eval.fd(pc.w,0:100/100)
  
  tmp.result <- t.test(as.formula(paste0(lead.tmp, "~afib")),data=data)
  
  ## plot data
  mean0 <- mean(pc.score.tmp)
  mean1 <- mean(pc.score.tmp[data$afib==0])
  mean2 <- mean(pc.score.tmp[data$afib==1])
  
  plotdat0 <- sqfunc(tmp1*mean0 + tmp0)
  plotdat1 <- sqfunc(tmp1*mean1 + tmp0)
  plotdat2 <- sqfunc(tmp1*mean2 + tmp0)
  
  plot(0:100/100, plotdat1,col="blue",type='l',lty=1, ylab="",xlab="", main=lead.tmp, ylim=c(min(c(plotdat0,plotdat1,plotdat2)),
                                                                                      max(c(plotdat0,plotdat1,plotdat2))) )
  lines(0:100/100, plotdat2,col="red",lty=2)    
  
}




for (j in 1:length(ecgvar.selected1)){
  pc.w.tmp <- gen.w.func(ecgvar.selected1[j], ecg.fpca.summary2)
  
  tmp <- sapply(ecgvar, gen.score.func, pc.w=pc.w.tmp, ecgfd=ecgfd_cs_train)
  df.tmp <- data.frame(SEQ=as.numeric(ecgfd_cs_train$I$ecgfd$fdnames$reps),
                       tmp )
  
  ## healthy subjects
  df.tmp.afib <- merge(df.tmp, list_ecgs_crosssectional_training, by="SEQ")

  

  par(mfrow=c(4,3))
  for ( i in 1:length(ecgvar) ){
    gen.plot.func(ecgvar[i], pc.w.tmp, ecgfd_cs_train, df.tmp.afib)
  }
}

#####################################################################
### Step 4: analysis of change in training data
#####################################################################
for (j in 1:length(ecgvar.selected1)){
  pc.w.tmp <- gen.w.func(ecgvar.selected1[j], ecg.fpca.summary2)
  
  lead.tmp <- gsub("\\..*","",ecgvar.selected1[j])
  tmp <- t( inprod( pc.w.tmp, center.fd( ecgfd_change_train[[lead.tmp]]$ecgfd ) ) )
  
  if (j==1){
    df.pair <- data.frame(SEQ=as.numeric(ecgfd_change_train$I$ecgfd$fdnames$reps), tmp)
  } else {
    df.pair <- cbind(df.pair,tmp)
  }
}
colnames(df.pair) <- c("SEQ", ecgvar.selected1)

alldat.change.train <- merge(df.pair,list_ecgs_change_training,by="SEQ")
alldat.change.train1 <- merge(alldat.change.train, ecg_muse, by="SEQ", all.x = TRUE)

## calculate the average of the two ECG measures
alldat.change.train.first <- alldat.change.train1 %>%
  group_by(PID) %>%
  arrange(ecg.date) %>%
  filter(row_number()==1)


tmp.pred <- alldat.change.train.first[,c(ecgvar.selected1,ecg.muse.var.cont)]
tmp.pred.scale <- scale(tmp.pred)*(-1)
ecgvar.scale <- paste0(c(ecgvar.selected1,ecg.muse.var.cont), ".scale")
colnames(tmp.pred.scale) <- ecgvar.scale
alldat.change.train.first.1 <- cbind(alldat.change.train.first, tmp.pred.scale)

## calculate the difference between the two ECG measures
alldat.change.train.dif <- alldat.change.train1 %>%
  group_by(PID) %>%
  arrange(ecg.date) %>%
  dplyr::summarise(across(all_of(c(ecgvar.selected1,ecg.muse.var)),diff)) 

changelist <- paste0(c(ecgvar.selected1,ecg.muse.var.cont), ".change")
changelist.cat <- paste0(ecg.muse.var.cat,".change")
colnames(alldat.change.train.dif)  <- c("PID", changelist, changelist.cat)


tmp.pred <- alldat.change.train.dif[,changelist]
tmp.pred.scale <- scale(tmp.pred)*(-1)

changelist.scale <- paste0(changelist, ".scale")
colnames(tmp.pred.scale) <- changelist.scale
alldat.change.train.dif.1 <- cbind(alldat.change.train.dif, tmp.pred.scale)

alldat.change.train.both <- merge(alldat.change.train.first.1, alldat.change.train.dif.1, by="PID")

## baseline covariates
cov.std.list <- c("AGE", "race2", "race3", "race4", "female", "SMOKENOW",
                  "BMI", "SYSTOLIC", "DIASTOLIC", "DIABETES", "EGFR_CRIC", "logpcr")
cov.data <- cricdata.phs1 %>%
  select(PID, all_of(cov.std.list), TIME_AFIB)  

alldat.change.train.both.1 <- merge(alldat.change.train.both, cov.data, by="PID")


res.cox64 <- coxph(as.formula(paste("Surv(TIME_AFIB,SA_AFIB)", paste(c(ecgvar.scale,changelist.scale[1]),collapse="+"), sep = "~") ) , 
                   data = alldat.change.train.both.1)
summary(res.cox64)