# Manuscript step 7 for its additive treatment/donor model. Fit using sufficient
# binomial cell counts, then restore the original pBC-row frequencies for ggeffect.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L, requireNamespace("brglm2",quietly=TRUE))
stage <- args[1]; refs <- args[2]
d <- read.delim(file.path(stage,"library_counts.tsv"),stringsAsFactors=FALSE)
stage_table <- c(Dedud=2L,Unfiltered=3L,Dedup=4L,"Dedud+Dedup"=5L)
arms <- c("NoEP","HiFi","WTCas9")
d$treat <- factor(ifelse(d$arm=="WT","WTCas9",d$arm),levels=arms)
d$donor <- factor(d$donor)
key <- c("target","strategy","donor","replicate","arm")
stopifnot(!anyDuplicated(d[key]), !anyNA(d), nrow(d)==26*79*4, all(d$indel>=0 & d$denominator>=d$indel),
          all(d$denominator>=d$family_rows & d$family_rows>=0))
tbls <- lapply(stage_table,function(i) read.csv(file.path(refs,paste0("SupplementaryTable",i,".csv")),check.names=FALSE))
tabulate_rates <- function(X, beta, V, tags) {
    eta <- drop(X %*% beta)
    se <- sqrt(rowSums((X %*% V)*X))
    values <- list(predicted=plogis(eta),stdError=se,confHigh=plogis(eta+qnorm(.975)*se),confLow=plogis(eta-qnorm(.975)*se))
    unlist(lapply(names(values),function(k) setNames(values[[k]],paste(k,tags,sep="-"))))
}
targets <- sort(unique(d$target))
stopifnot(length(targets)==79)
for (tab in tbls) stopifnot(nrow(tab)==79,!anyDuplicated(tab$barcode),setequal(tab$barcode,targets))
rows <- list(); details <- list(); weights <- list()
for (s in names(stage_table)) for (target in targets) {
    sub <- d[d$strategy==s & d$target==target,]
    x <- aggregate(cbind(indel,denominator,family_rows)~donor+treat,sub,sum)
    stopifnot(nrow(x)==9,all(x$denominator>0))
    warns <- character()
    model <- withCallingHandlers(glm(cbind(indel,denominator-indel)~treat+donor,data=x,family=binomial(),
                    method=brglm2::brglmFit,type="AS_mean",control=brglm2::brglmControl(maxit=300,epsilon=1e-9)),
                    warning=function(w){warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
    stopifnot(model$converged,all(is.finite(coef(model))))
    X <- model.matrix(model); beta <- coef(model); V <- vcov(model)
    names_treat <- paste0("treat",arms[-1])
    z <- beta[names_treat]/sqrt(diag(V)[names_treat])
    # Match the submitted one-sided Wald convention, including zeros from subtraction.
    p <- 1-pnorm(z)
    Xmar <- X[match(arms,x$treat),,drop=FALSE]
    donor_cols <- startsWith(colnames(X),"donor")
    for (j in which(donor_cols)) Xmar[,j] <- weighted.mean(X[,j],x$family_rows)
    donor_weights <- tapply(x$family_rows,x$donor,sum)/sum(x$family_rows)
    id <- paste(target,s)
    weights[[id]] <- data.frame(target=target,strategy=s,donor=names(donor_weights),
                family_row_weight=as.numeric(donor_weights),iterations=model$iter,warning=paste(unique(warns),collapse=" | "))
    calculated <- c(setNames(beta,paste0("estim-",names(beta))),setNames(z,paste0("zScore-",names(z))),
        setNames(p,paste0("asymPval-",names(p))),tabulate_rates(X,beta,V,paste(x$donor,x$treat,sep="-")),
        tabulate_rates(Xmar,beta,V,arms))
    tab <- tbls[[s]]
    reference <- unlist(tab[match(target,tab$barcode),setdiff(names(tab),"barcode"),drop=FALSE],use.names=TRUE)
    stopifnot(setequal(names(reference),names(calculated)))
    reference <- reference[names(calculated)]
    details[[id]] <- data.frame(target=target,strategy=s,table=stage_table[[s]],field=names(calculated),
        agent=unname(calculated),submitted=unname(reference),difference=unname(calculated-reference))
    for (a in arms[-1]) {
        pick <- function(vec,arm) unname(vec[paste0("predicted-",arm)])*100
        ctrl <- pick(calculated,"NoEP"); treated <- pick(calculated,a)
        ref_ctrl <- pick(reference,"NoEP"); ref_treated <- pick(reference,a)
        term <- paste0("asymPval-treat",a)
        rows[[paste(id,a)]] <- data.frame(target=target,strategy=s,editor=ifelse(a=="WTCas9","WT",a),
            supp_table=stage_table[[s]],control_rate_pct=ctrl,edited_rate_pct=treated,delta_pp=treated-ctrl,
            manuscript_control_pct=ref_ctrl,manuscript_edited_pct=ref_treated,manuscript_delta_pp=ref_treated-ref_ctrl,
            difference_pp=(treated-ctrl)-(ref_treated-ref_ctrl),p_one_sided=unname(calculated[term]),manuscript_p=unname(reference[term]))
    }
}
panel <- do.call(rbind,rows)
adjust <- function(p) ave(p,panel$strategy,panel$editor,FUN=function(x) p.adjust(x,"BH"))
panel$q_bh <- adjust(panel$p_one_sided); panel$manuscript_q_bh <- adjust(panel$manuscript_p)
hit <- function(delta,q) delta>.1 & q<.05
panel$agent_hit <- hit(panel$delta_pp,panel$q_bh)
panel$submitted_hit <- hit(panel$manuscript_delta_pp,panel$manuscript_q_bh)
panel$hit_agrees <- panel$agent_hit==panel$submitted_hit
summary <- do.call(rbind,lapply(split(panel,list(panel$strategy,panel$editor),drop=TRUE),function(z) {
    abs_diff <- abs(z$difference_pp)
    data.frame(strategy=z$strategy[1],editor=z$editor[1],targets=nrow(z),agent_hits=sum(z$agent_hit),
       submitted_hits=sum(z$submitted_hit),hit_agreement=sum(z$hit_agrees),agent_only=sum(z$agent_hit & !z$submitted_hit),
       submitted_only=sum(!z$agent_hit & z$submitted_hit),median_abs_delta_diff_pp=median(abs_diff),
       max_abs_delta_diff_pp=max(abs_diff),largest_difference_target=z$target[which.max(abs_diff)])
}))
export <- function(x,name) write.table(x,file.path(stage,paste0(name,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
export(panel,"panel_comparison"); export(summary,"comparison_summary")
export(do.call(rbind,details),"model_field_comparison"); export(do.call(rbind,weights),"donor_weights")
writeLines(c(paste("brglm2",packageVersion("brglm2")), "Model: Firth logistic regression, treatment + donor, NoEP reference.",
    "Treatment margin: inverse logit after averaging donor design columns over original sample/target/pBC rows across treatments.",
    "Tests: one-sided coefficient Wald. BH within each 79-target editor/strategy combination.",
    "The field-level comparison covers submitted model coefficients, z, p, donor/treatment and treatment estimates, SE and 95% limits.",
    "The separate simulation-based power columns are outside the indel-rate comparison."),file.path(stage,"model_info.txt"))
