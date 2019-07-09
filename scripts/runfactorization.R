##### runfactoization runs all the considered multi-omics factorization
### the required inputs are:
### "folder" corresponding to the path to the folder where the input files are contained, in the idea that all omics matrices are organized inside a unique folder
### "file.names" corresponding to a vector containing the names of all the omics files
### "num.factors" containing the number of factors in which we desire to decompose the matrices
### "sep=" "" corresponding to the separator used in the omics files required to properly read them
### "single.cell" indicating if the data are single cell data. In this case the filtering of the data will be more intense in respect to other data types
### "filtering"

### the input files need to be log2 transformed before running the analysis


library("RGCCA")
library("r.jive")
library("IntNMF")
library("omicade4")
library("MSFA")
library("GPArotation")
library("MOFAtools")
library("tensorBSS")
source("tICA.R")
library("iCluster")

runfactorization <- function(folder,file.names,num.factors,sep=" ",filtering="none"){
  factorizations<-list()
  t<-1
  method<-numeric(0)
  
  num.factors<-as.numeric(num.factors)

  
  ##creating list of omics
  omics <- list()
  for(i in 1:length(file.names)){
    omics[[i]]<-as.matrix(read.table(paste(folder,file.names[i],sep="/"),sep=sep,row.names=1,header=T))
  }
  
  ####
  #omics<-lapply(omics, function(x) t(x))
  ######
  
  ##restricting to common samples and filtering
  samples<-colnames(omics[[1]])
  for(j in 1:length(omics)){
    samples<-intersect(samples,colnames(omics[[j]]))
  }
  for(j in 1:length(omics)){
    omics[[j]]<-omics[[j]][,samples]
    if(filtering!="none"){
      x<-apply( omics[[j]],1,sd)
      x<-as.matrix(sort(x, decreasing = T))
      w<-which(x>0)
      if(filtering=="stringent"){
        selected<-rownames(x)[1:min(w[length(w)],5000)]
      }else{
        selected<-rownames(x)[1:min(w[length(w)],6000)]
      }
      m<-match(rownames(omics[[j]]),selected)
      w<-which(!is.na(m))
      omics[[j]]<-omics[[j]][w,]
    }else{
      omics[[j]]<-omics[[j]][which(apply(omics[[j]],2,sd)>0),]
    }
  }  
  
 
  ##RGCCA 
  factorizations_RGCCA<-rgcca(lapply(omics, function(x) t(x)), ncomp = rep(num.factors, length(omics)), scheme = "centroid", scale = TRUE, init = "svd",bias = TRUE, tol = 1e-08, verbose = F)
  factors_rgcca<-as.matrix(factorizations_RGCCA$Y[[1]])
  metagenes_rgcca <- list()
  for(j in 1:length(omics)){
    metagenes_rgcca[[j]]<-as.matrix(factorizations_RGCCA$a[[j]])
    rownames(metagenes_rgcca[[j]])<-rownames(omics[[j]])
    colnames(metagenes_rgcca[[j]])<-1:num.factors
  }
  factorizations[[t]]<-list(factors_rgcca,metagenes_rgcca)
  t<-t+1
  method<-c(method,"RGCCA")
  
  ###MCIA
  omics_pos<-list()
  for(j in 1:length(omics)){
    if(min(omics[[j]])<0){
      omics_pos[[j]]<-omics[[j]]+abs(min(omics[[j]]))
    }else{
      omics_pos[[j]]<-omics[[j]]
    }
    omics_pos[[j]]<-omics_pos[[j]]/max(omics_pos[[j]])
  }
  factorizations_mcia<-mcia(omics_pos, cia.nf = num.factors)
  factors_mcia<-as.matrix(factorizations_mcia$mcoa$SynVar)
  metagenes_mcia<-list()
  for(j in 1:length(omics)){
    metagenes_mcia[[j]]<-as.matrix(factorizations_mcia$mcoa$axis[1:dim(omics[[j]])[1],])
    rownames(metagenes_mcia[[j]])<-rownames(omics[[j]])
    colnames(metagenes_mcia[[j]])<-1:num.factors
  }
  factorizations[[t]]<-list(factors_mcia,metagenes_mcia)
  t<-t+1
  method<-c(method,"MCIA")
  
  ###MOFA
  MOFAobject <- createMOFAobject(omics)
  DataOptions <- getDefaultDataOptions()
  ModelOptions <- getDefaultModelOptions(MOFAobject)
  ModelOptions$numFactors <- num.factors
  TrainOptions <- getDefaultTrainOptions()
  
  MOFAobject <- prepareMOFA(
    MOFAobject,
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )
  MOFAobject <- runMOFA(MOFAobject)
  metagenes_mofa <- getWeights(
    MOFAobject
  )
  
  factors_mofa <- getFactors(
    MOFAobject
  )
  factorizations[[t]]<-list(factors_mofa,metagenes_mofa)
  t<-t+1
  method<-c(method,"MOFA")
  
  # # ###MSFA
  # if(filtering=="none" ||  filtering=="stringent"){
  # start_value <- start_msfa(X_s = omics_pos, k = num.factors, j_s = c(num.factors,num.factors,num.factors))
  # mle <-  ecm_msfa(omics_pos, start_value, trace = FALSE)
  # factors_msfa<-as.matrix(mle$Phi)
  # metagenes_msfa<-list()
  # for(j in 1:length(omics)){
  #    metagenes_msfa[[j]]<-t(omics[[j]]) %*% solve(cor(t(omics[[j]]))) %*% factors_msfa
  #    rownames(metagenes_msfa[[j]])<-rownames(omics[[j]])
  #    colnames(metagenes_msfa[[j]])<-1:num.factors
  #  }
  #  factorizations[[t]]<-list(factors_msfa,metagenes_msfa)
  #  t<-t+1
  # method<-c(method,"MSFA")
  # }
  
  ###iCluster
  factorizations_icluster<-iCluster2(lapply(omics, function(x) t(x)), k=num.factors+1)
  factors_icluster<-as.matrix(t(factorizations_icluster$expZ))
  metagenes_icluster<-list()
  for(j in 1:length(omics)){
    metagenes_icluster[[j]]<-as.matrix(factorizations_icluster$W[1:dim(omics[[j]])[1],])
    rownames(metagenes_icluster[[j]])<-rownames(omics[[j]])
    colnames(metagenes_icluster[[j]])<-1:num.factors
  }
  factorizations[[t]]<-list(factors_icluster,metagenes_icluster)
  t<-t+1
  method<-c(method,"iCluster")
  
  ###intNMF
  factorizations_intnmf<-nmf.mnnals(dat=lapply(omics_pos, function(x) t(x)), k=num.factors)
  factors_intNMF<-as.matrix(factorizations_intnmf$W)
  colnames(factors_intNMF)<-1:num.factors
  metagenes_intNMF<-list()
  for(j in 1:length(omics)){
    metagenes_intNMF[[j]]<-t(factorizations_intnmf$H[[j]]) 
    rownames(metagenes_intNMF[[j]])<-rownames(omics[[j]])
    colnames(metagenes_intNMF[[j]])<-1:num.factors
  }
  factorizations[[t]]<-list(factors_intNMF,metagenes_intNMF)
  t<-t+1
  method<-c(method,"intNMF")
  
  
  ### JIVE
  factorizations_jive<-jive(omics, rankJ=num.factors, rankA = rep(num.factors, length(omics)), method = "given", conv = "default", maxiter = 100, showProgress=FALSE)
  rankJV <- factorizations_jive$rankJ;
  rankIV.v <- factorizations_jive$rankA;
  J<-numeric(0)
  ng<-0
  metagenes_jive <- list();
  for(j in 1:length(omics)){
    J <- rbind(J,factorizations_jive$joint[[j]]);
    ng<-c(ng,dim(factorizations_jive$joint[[j]])[1])
  }
  svd.o <- svd(J);
  jV <- svd.o$v %*% diag(svd.o$d);
  for(j in 1:length(omics)){
    metagenes_jive[[j]] <- svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]; ###error in dimension
    rownames(metagenes_jive[[j]])<-rownames(omics[[j]])
    colnames(metagenes_jive[[j]])<-1:num.factors
  }
  factors_jive=jV[,1:rankJV]
  rownames(factors_jive)<-colnames(omics[[1]])
  colnames(factors_jive)<-1:num.factors
  factorizations[[t]]<-list(factors_jive,metagenes_jive)
  t<-t+1
  method<-c(method,"JIVE")
  

  

  


  

  
  ###tICA
    omics_tensor<-list()
    for(j in 1:length(omics)){
      omics_tensor[[j]]<-cor(omics[[j]], method = "spearman")
    }
    
    S<-vector(length = dim(omics[[1]])[2]*dim(omics[[1]])[2]*length(omics))
    dim(S) <- c(length(omics), dim(omics[[1]])[2], dim(omics[[1]])[2])
    for(j in 1:length(omics)){
      S[j,,]<-t(omics_tensor[[j]])
    }
    tICA<-DoTICA(S,num.factors,method="FOBI")
    factors_tica<-as.matrix(tICA$signals)
    rownames(factors_tica)<-colnames(omics[[1]])
    metagenes_tica<-list()
    for(j in 1:length(omics)){
      metagenes_tica[[j]]<-omics[[j]] %*% ginv(t(tICA$signals))
      rownames(metagenes_tica[[j]])<-rownames(omics[[j]])
      colnames(metagenes_tica[[j]])<-1:num.factors
    }
    factorizations[[t]]<-list(factors_tica,metagenes_tica)
    t<-t+1
    method<-c(method,"tICA")


  


  
  # ###scikit-fusion
  system("mkdir ../data/scikit/")
  temp.folder<-'../data/scikit/'
  files<-""
  for(j in 1:length(omics)){
    write.table(omics[[j]],paste(temp.folder,"/omics",j,".txt",sep=""),sep=" ",col.names=T, row.names=T)
    files<-paste(files,paste("omics",j,".txt",sep=""),sep=" ")
  }
  system(paste("python ./scikit_fusion.py", "'../data/scikit/'","'../data/scikit/'",num.factors,"' '",files,sep=" "))
  factors_scikit<-as.matrix(read.table(paste(temp.folder,"signals.txt",sep=""),sep="\t",header=F))
  colnames(factors_scikit)<-1:num.factors
  rownames(factors_scikit)<-colnames(omics[[1]])
  metagenes_scikit<-list()
  for(j in 1:length(omics)){
    metagenes_scikit[[j]]<-as.matrix(read.table(paste(temp.folder,"projomics",j,".txt",sep=""),sep="\t",header=F))
    rownames(metagenes_scikit[[j]])<-rownames(omics[[j]])
    colnames(metagenes_scikit[[j]])<-1:num.factors
  }
  factorizations[[t]]<-list(factors_scikit,metagenes_scikit)
  t<-t+1
  method<-c(method,"scikit-fusion")
  system("rm -r ../data/scikit/")
  
  
  out<-list(factorizations=factorizations,method=method,icluster.clusters=as.matrix(factorizations_icluster$clusters),intNMF.clusters=as.matrix(factorizations_intnmf$clusters))
  
  return(out)
}