# Run on a matched three-arm panel: Rscript report.R COUNTS.tsv NEW_OUTPUT_DIR [SUPPLEMENT_FOLDER]
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) %in% 2:3) stop("Usage: Rscript report.R COUNTS.tsv NEW_OUTPUT_DIR [SUPPLEMENT_FOLDER]")
input <- normalizePath(args[1])
output <- args[2]
if (file.exists(output)) stop("Choose a new output directory")
if (!requireNamespace("brglm2",quietly=TRUE)) stop("brglm2 is required")
d <- read.delim(input,stringsAsFactors=FALSE,check.names=FALSE)
arms <- c("NoEP","HiFi","WT")
stages <- c("Dedud","Unfiltered","Dedup","Dedud+Dedup")
samplekey <- c("donor","replicate","target","strategy")
fields <- unlist(lapply(arms,function(c) paste0(c,c("_indel","_denominator"))))
stopifnot(all(c(samplekey,fields) %in% names(d)), !anyNA(d[samplekey]), !anyDuplicated(d[samplekey]))
targets <- sort(unique(d$target))
pairs <- unique(d[c("donor","replicate")])
stopifnot(length(targets)==79, setequal(unique(d$strategy),stages),length(unique(d$donor))>=2,
          nrow(d)==nrow(pairs)*length(targets)*length(stages))
for(c in arms) {
  n <- d[[paste0(c,"_denominator")]]; k <- d[[paste0(c,"_indel")]]
  stopifnot(all(is.finite(n)&is.finite(k)&n>=k&k>=0&k==floor(k)&n==floor(n)))
}
don <- aggregate(d[fields],by=d[c("target","strategy","donor")],FUN=sum)
stratum <- function(t,s) paste(t,s,sep="|")
g <- split(d,stratum(d$target,d$strategy)); gd <- split(don,stratum(don$target,don$strategy))
rows <- list(); evidence <- list(); diagnostics <- list()
for (s in stages) for (t in targets) {
  donors <- gd[[stratum(t,s)]]; pair <- g[[stratum(t,s)]]
  cells <- do.call(rbind,lapply(arms,function(c) data.frame(donor=donors$donor,arm=c,
         k=donors[[paste0(c,"_indel")]],n=donors[[paste0(c,"_denominator")]])))
  cells$arm <- factor(cells$arm,levels=arms)
  if(any(cells$n==0)) stop(paste("No donor/condition coverage",t,s))
  notes <- character()
  m <- withCallingHandlers(glm(cbind(k,n-k)~0+arm+donor,data=cells,family=binomial(),method=brglm2::brglmFit,type="AS_mean",control=brglm2::brglmControl(maxit=300,epsilon=1e-9)),
    warning=function(w){notes <<- c(notes,conditionMessage(w));invokeRestart("muffleWarning")})
  diagnostics[[stratum(t,s)]] <- data.frame(target=t,strategy=s,converged=m$converged, iterations=m$iter,warning=paste(unique(notes),collapse=" | "))
  stopifnot(m$converged)
  X <- model.matrix(m);pred <- fitted(m);V <- vcov(m);co <- coef(m)
  mu <- tapply(pred,cells$arm,mean)
  grad <- function(c) {
     j <- cells$arm==c
     colMeans(X[j,,drop=FALSE]*as.numeric(pred[j]*(1-pred[j])))
  }
  src <- if(s %in% c("Dedud","Dedud+Dedup")) c("Dedud","Dedud+Dedup") else c("Unfiltered","Dedup")
  reads <- gd[[stratum(t,src[1])]]; fam <- gd[[stratum(t,src[2])]]
  for (c in arms[-1]) {
    diff <- 100*unname(mu[c]-mu["NoEP"])
    jac <- 100*(grad(c)-grad("NoEP"))
    se <- sqrt(drop(t(jac)%*%V%*%jac))
    contrast <- setNames(numeric(length(co)),names(co))
    contrast[paste0("arm",c)] <- 1; contrast["armNoEP"] <- -1
    z <- drop(sum(contrast*co)/sqrt(drop(t(contrast)%*%V%*%contrast)))
    # Cross-products preserve ties between count fractions even at tiny rates.
    positive <- function(a) {
      n <- as.double(a[[paste0(c,"_denominator")]]); n0 <- as.double(a$NoEP_denominator)
      k <- as.double(a[[paste0(c,"_indel")]]); k0 <- as.double(a$NoEP_indel)
      valid <- n>0 & n0>0
      stopifnot(all(k*n0<2^53),all(k0*n<2^53))
      ifelse(valid,sign(k*n0-k0*n),NA_real_)
    }
    sg <- positive(pair)
    id <- paste(t,s,c,sep="|")
    rows[[id]] <- data.frame(target=t,strategy=s,editor=c,control_rate_pct=100*unname(mu["NoEP"]),
      edited_rate_pct=100*unname(mu[c]),delta_pp=diff,ci95_low_pp=diff-qnorm(.975)*se,
      ci95_high_pp=diff+qnorm(.975)*se,wald_z=z,p_one_sided=pnorm(z,lower.tail=FALSE),
      pairs_valid=sum(!is.na(sg)),pairs_positive=sum(sg>0,na.rm=TRUE),pairs_zero=sum(sg==0,na.rm=TRUE),
      pairs_negative=sum(sg<0,na.rm=TRUE),donors_positive=sum(positive(donors)>0),
      edited_indel_reads=sum(reads[[paste0(c,"_indel")]]),edited_total_reads=sum(reads[[paste0(c,"_denominator")]]),
      edited_indel_families=sum(fam[[paste0(c,"_indel")]]),edited_total_families=sum(fam[[paste0(c,"_denominator")]]),
      noep_indel_reads=sum(reads$NoEP_indel),noep_total_reads=sum(reads$NoEP_denominator),
      noep_indel_families=sum(fam$NoEP_indel),noep_total_families=sum(fam$NoEP_denominator))
    evidence[[id]] <- do.call(rbind,lapply(donors$donor,function(dd){
      i <- reads$donor==dd;j <- fam$donor==dd;sgd <- sg[pair$donor==dd]
      data.frame(target=t,strategy=s,editor=c,donor=dd,pairs_valid=sum(!is.na(sgd)),pairs_positive=sum(sgd>0,na.rm=TRUE),
        edited_indel_reads=reads[[paste0(c,"_indel")]][i],edited_total_reads=reads[[paste0(c,"_denominator")]][i],
        edited_indel_families=fam[[paste0(c,"_indel")]][j],edited_total_families=fam[[paste0(c,"_denominator")]][j],
        noep_indel_reads=reads$NoEP_indel[i],noep_total_reads=reads$NoEP_denominator[i],
        noep_indel_families=fam$NoEP_indel[j],noep_total_families=fam$NoEP_denominator[j])
    }))
  }
}
panel <- do.call(rbind,rows);rownames(panel)<-NULL
adjust <- function(p) ave(p,panel$strategy,panel$editor,FUN=function(x) p.adjust(x,method="BH"))
panel$q_bh <- adjust(panel$p_one_sided)
band <- function(delta,q) ifelse(q>=.05,"q>=0.05",ifelse(delta>.1,"q<0.05; >0.1 pp",
                    ifelse(delta>.01,"q<0.05; >0.01 to 0.1 pp","q<0.05; <=0.01 pp")))
panel$effect_band <- band(panel$delta_pp,panel$q_bh)
if(length(args)==3) {
  tabs <- setNames(2:5,stages)
  reference <- do.call(rbind,lapply(stages,function(s){
    z <- read.csv(file.path(args[3],paste0("SupplementaryTable",tabs[[s]],".csv")),check.names=FALSE)
    stopifnot(!anyDuplicated(z$barcode),setequal(z$barcode,targets))
    do.call(rbind,lapply(arms[-1],function(c){
      tag <- if(c=="WT") "WTCas9" else c
      data.frame(target=z$barcode,strategy=s,editor=c,supp_table=tabs[[s]],
        manuscript_delta_pp=100*(z[[paste0("predicted-",tag)]]-z[["predicted-NoEP"]]),
        manuscript_p=z[[paste0("asymPval-treat",tag)]])
    }))
  }))
  id <- function(x) paste(x$target,x$strategy,x$editor)
  ref <- reference[match(id(panel),id(reference)),]
  stopifnot(identical(id(ref),id(panel)),all(is.finite(ref$manuscript_p)),all(ref$manuscript_p>=0 & ref$manuscript_p<=1))
  panel <- cbind(panel,ref[c("supp_table","manuscript_delta_pp","manuscript_p")])
  panel$manuscript_q_bh <- adjust(panel$manuscript_p)
  panel$change_from_manuscript_pp <- panel$delta_pp-panel$manuscript_delta_pp
  panel$manuscript_effect_band <- band(panel$manuscript_delta_pp,panel$manuscript_q_bh)
}
parent <- dirname(output)
if(!dir.exists(parent)) dir.create(parent,recursive=TRUE)
stage <- tempfile("panel-",tmpdir=parent);dir.create(stage)
export <- function(data,name) write.table(data,file.path(stage,name),sep="\t",quote=FALSE,row.names=FALSE,na="")
export(panel,"panel_evidence.tsv")
export(do.call(rbind,diagnostics),"model_diagnostics.tsv")
bydonor <- do.call(rbind,evidence)
report_id <- function(x) paste(x$target,x$strategy,x$editor)
matched <- match(report_id(bydonor),report_id(panel))
stopifnot(!anyNA(matched))
bydonor <- cbind(bydonor,panel[matched,c("delta_pp","ci95_low_pp","ci95_high_pp","q_bh","effect_band")])
export(bydonor,"donor_evidence.tsv")
writeLines(c(paste("counts",input),paste("brglm2",packageVersion("brglm2")),
  paste("complete_pairs",nrow(pairs)),paste("donors",paste(sort(unique(d$donor)),collapse=",")),
  "95% intervals: pointwise response-scale delta-method Wald. One-sided coefficient Wald tests for increases.",
  "BH family: 79 targets per editor and strategy. Family support is summed over libraries.",
  "Coverage-conditional binomial inference; pBC family counts and replicate direction accompany each estimate."),file.path(stage,"model_info.txt"))
stopifnot(!file.exists(output),file.rename(stage,output))
