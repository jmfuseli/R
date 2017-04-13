# This script adjusts the imaging phenotypes by fitting a linear model
# to each phenotype as a function of the covariates and returns the
# residuals as the adjusted phenotype

#-------------------------------------------------------------------------

#get the phenotype file names
pheno.file.names = dir(path, pattern = "_pheno.txt")

#take off the extension to make the new file name
fnames = strsplit(pheno.file.names, ".txt")

#get the covariate file names
cov.file.names = dir(path, pattern = "_COV.txt")

#read in the fam file
fam = read.table(file.choose(), header = FALSE)

adjust_phenos = function( n ){
  #read in files and set missing info to NA
  phenos = read.table(pheno.file.names[n], header = TRUE, na.strings = -9)
  cov = read.table(cov.file.names[n], header = TRUE, na.strings = -9)
  
  #check for NAs
  cov_NA_ind = which(is.na(cov), arr.ind=TRUE)
  pheno_NA_ind = which(is.na(phenos), arr.ind=TRUE)
  
  if ((length(cov_NA_ind) == 0) && (length(pheno_NA_ind) == 0)) {
    #take the intersection of the IIDs in the phenotype and covariate files
    same_ids = intersect(phenos$IID, cov$IID)
    
    #take the intersection with the fam file IID; these are the IDs plink analyzes
    ids = intersect(same_ids, fam$V2)
    
    #subset cov and pheno files on these ids and then run linear regression
    final_phenos = phenos[as.character(phenos$IID) %in% ids,]
    final_cov = cov[as.character(cov$IID) %in% ids,]
    
    new_phenos = cbind(final_phenos$FID, as.character(final_phenos$IID))
    #run linear regression: f(phenotype) = covariates 
    for(i in 3:dim(final_phenos)[2]){
      pheno = final_phenos[,i]
      pheno.lm = lm(pheno ~ final_cov$Age + final_cov$ICV + final_cov$Gender, 
                    data = final_cov)
      new_pheno = pheno.lm$residuals
      new_phenos = cbind(new_phenos, new_pheno)
    }
  }
  #make column names
  colnames(new_phenos) = colnames(final_phenos)
  out_fname = paste(fnames[[n]], "mv_plink.txt", sep="_")
  
  #save residuals as the adjusted phenotype for input in MV-Plink
  write.table(new_phenos, file = out_fname, row.names=FALSE, quote=FALSE)
}

for (i in 1:length(pheno.file.names)){
  adjust_phenos(i) 
}


