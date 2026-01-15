Impgeno <- read.csv(paste0(pathIn, 'GenoDosageData.csv'), header=T)
SitesPassing <- read.table(paste0(pathIn, 'Allsitesfreq.frq'), header=T)
SitesPassing$pq <- SitesPassing$A1*SitesPassing$A2 
ImpGeno <- as.matrix(Impgeno[,-1])
Z  <- scale(ImpGeno, scale=FALSE)
ZZ <- Z%*%t(Z)
fq <- sum(SitesPassing$pq)
GRM <- ZZ/(2*fq)
print(mean(diag(GRM)))
qr(GRM)$rank - ncol(GRM)         # this should be 0 if fGRM is full rank,but its not so, it will be -1 if it has rank 1 less than geno
diag(GRM) <- diag(GRM) + 0.01   # try this "bending" and see if it fixes the fit.
qr(GRM)$rank - ncol(GRM)        # this should be 0 if mGRM is full rank after the 'bending'.
print(mean(diag(GRM)))
#####
colnames(GRM) <- rownames(GRM) <- Impgeno[,1]
GCAgrm <- GRM/2
print(mean(diag(GCAgrm)))

writeGrm<-function(G){
  nLines=length(G[1,])
  print(nLines)
  for(i in c(1:nLines)){
    for(j in c(1:i)){
      forout=c(i,j,G[i,j])
      if(i==1&&j==1){
        write(t(forout),file=paste0(pathIn, "femGRM.grm", sep=""), ncol=3)
      } else {
        write(t(forout),file=paste0(pathIn, "femGRM.grm", sep=""), ncol=3,append=TRUE)
      }
    }
  }
}
writeGrm(GCAgrm)