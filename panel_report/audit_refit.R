# Usage: Rscript audit_refit.R TRACE_DIRECTORY
folder <- commandArgs(trailingOnly=TRUE)[1]
read_tsv <- function(stem) read.delim(file.path(folder,paste0(stem,".tsv")),stringsAsFactors=FALSE)
d <- read_tsv("sample_counts")
families <- read_tsv("family_presence")
fit_once <- function(rows,scenario) {
  x <- aggregate(cbind(indel_reads,reads)~donor+condition,rows,sum)
  x$condition <- factor(x$condition,levels=c("NoEP","HiFi","WT"))
  stopifnot(nrow(x)==9,all(x$reads>0),all(x$indel_reads>=0 & x$indel_reads<=x$reads))
  m <- glm(cbind(indel_reads,reads-indel_reads)~0+condition+donor,data=x,family=binomial(),
           method=brglm2::brglmFit,type="AS_mean",control=brglm2::brglmControl(maxit=300,epsilon=1e-9))
  stopifnot(m$converged)
  x$rate <- fitted(m)*100
  ds <- c(sort(unique(x$donor)),"equal_donor_mean")
  do.call(rbind,lapply(ds,function(donor){
    z <- if(donor=="equal_donor_mean") x else x[x$donor==donor,]
    p <- tapply(z$rate,z$condition,mean)
    data.frame(target=rows$target[1],scenario=scenario,donor=donor,noep_pct=unname(p["NoEP"]),
      wt_pct=unname(p["WT"]),wt_delta_pp=unname(p["WT"]-p["NoEP"]))
  }))
}
result <- list()
for(t in unique(d$target)) {
  for(rule in unique(d$rule)) result[[length(result)+1]] <- fit_once(d[d$target==t & d$rule==rule,],rule)
  a <- d[d$target==t & d$rule=="published_pBC",]
  f <- families[families$target==t,]
  for(fam in unique(f$family)) {
    rr <- f[f$family==fam,];id <- function(x) paste(x$donor,x$replicate,x$condition)
    ix <- match(id(a),id(rr));stopifnot(!anyNA(ix))
    tmp <- a;tmp$reads <- tmp$reads-rr$aligned_reads[ix];tmp$indel_reads <- tmp$indel_reads-rr$indel_reads[ix]
    result[[length(result)+1]] <- fit_once(tmp,paste0("omit_",fam))
  }
}
write.table(do.call(rbind,result),file.path(folder,"sensitivity.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
