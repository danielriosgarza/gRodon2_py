# Regenerate portable numeric parity fixtures with R itself.
# Usage: Rscript scripts/generate_ioc_fitting_reference.R <gRodon2> <output-dir>
# Requires only base/recommended R packages (MASS); uses synthetic numeric traits.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Expected the R package directory and output directory")
source_root <- normalizePath(args[1], winslash="/")
output <- args[2]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(MASS))
source(file.path(source_root, "R/index_of_copiotrophy.R"))
loaded <- new.env(parent=baseenv())
load(file.path(source_root, "R/sysdata.rda"), envir=loaded)
IoC_models <- loaded$IoC_models
options(digits=17)

traits <- data.frame(dCUB=c(-.1,-.2,-.5,-.3,-.04,-.7),
                     nCAZy=c(10,50,80,45,3,130),
                     nGenes=c(1500,2200,4000,3000,800,5100))
heldout <- data.frame(dCUB=c(-.2,-.5,-.01), nCAZy=c(30,90,0), nGenes=c(2000,4000,500))
fit <- fitIoC(traits)
fit_complete <- fitIoC(traits, completeness=c(90,100,50,75,100,80))
write.csv(traits, file.path(output, "ioc_traits.csv"), row.names=FALSE)
write.csv(heldout, file.path(output, "ioc_heldout.csv"), row.names=FALSE)
write.csv(data.frame(IoC=fit$IoC, complete_IoC=fit_complete$IoC),
          file.path(output, "ioc_fit.csv"), row.names=FALSE)
predictions <- data.frame(custom=predictIoC(heldout,model=fit))
for (name in c("permafrost","pacific","soil","human","all_habitat")) {
  predictions[[name]] <- predictIoC(heldout,model=name)
}
write.csv(predictions, file.path(output,"ioc_predictions.csv"), row.names=FALSE)
write.csv(data.frame(variance=fit$model$sdev^2, center=fit$model$center, scale=fit$model$scale),
          file.path(output,"ioc_pca.csv"), row.names=FALSE)

i <- seq_len(60)
training <- data.frame(
  CUBHE=.25+.1*sin(i*.7), ConsistencyHE=.1+.04*cos(i*.37), CPB=.02+.01*sin(i*1.3),
  OGT=15+(i %% 30), CUB=.2+.03*cos(i*.25), GCdiv=.05+.01*sin(i*.13),
  MILC=.8+.15*sin(i*.21), ENCprime=40+2*cos(i*.5), B=.2+.1*sin(i*.41),
  SCUO=.25+.1*cos(i*.37), MCB=.3+.05*cos(i*.71),
  d=exp(1.2+.4*sin(i*.7)-.2*cos(i*.37)+.1*sin(i*1.3)+.13*sin(i*2.2))
)
extremo <- training
extremo$d <- training$d*exp(.005*(extremo$OGT-30))
write.csv(training,file.path(output,"training.csv"),row.names=FALSE)
write.csv(extremo,file.path(output,"training_extremo.csv"),row.names=FALSE)
training$dCUB <- (training$CUB-training$CUBHE)/training$CUB
extremo$dCUB <- (extremo$CUB-extremo$CUBHE)/extremo$CUB
pdf(file.path(output,"boxcox-diagnostic.pdf"))
best_lambda <- function(terms) {
  profile <- MASS::boxcox(reformulate(terms, response="d"), data=training)
  profile$x[which.max(profile$y)]
}
lambda_full <- best_lambda(c("CUBHE","ConsistencyHE","CPB"))
lambda_euk <- best_lambda("CUBHE")
lambda_gc <- best_lambda(c("dCUB","GCdiv"))
lambda_trim <- best_lambda(c("MILC","ENCprime","B","SCUO","MCB"))
dev.off()
unlink(file.path(output,"boxcox-diagnostic.pdf"))
models <- list(
  gRodon_model_base=list(c("CUBHE","ConsistencyHE","CPB"),lambda_full,FALSE),
  gRodon_model_temp=list(c("CUBHE","ConsistencyHE","CPB","OGT"),lambda_full,TRUE),
  gRodon_model_partial=list(c("CUBHE","ConsistencyHE"),lambda_full,FALSE),
  gRodon_model_partial_temp=list(c("CUBHE","ConsistencyHE","OGT"),lambda_full,TRUE),
  gRodon_model_meta=list(c("CUBHE"),lambda_full,FALSE),
  gRodon_model_meta_temp=list(c("CUBHE","OGT"),lambda_full,TRUE),
  gRodon_model_euk=list(c("CUBHE"),lambda_euk,FALSE),
  gRodon_model_euk_temp=list(c("CUBHE","OGT"),lambda_euk,TRUE),
  meta_model_base=list(c("dCUB","GCdiv"),lambda_gc,FALSE),
  meta_model_temp=list(c("dCUB","GCdiv","OGT"),lambda_gc,TRUE),
  meta_nogc_model_base=list(c("dCUB"),lambda_gc,FALSE),
  meta_nogc_model_temp=list(c("dCUB","OGT"),lambda_gc,TRUE),
  meta_old_model_base=list(c("CUBHE"),lambda_euk,FALSE),
  meta_old_model_temp=list(c("CUBHE","OGT"),lambda_euk,TRUE),
  gRodon_model_trim=list(c("MILC","ENCprime","B","SCUO","MCB"),lambda_trim,FALSE),
  gRodon_model_trim_temp=list(c("MILC","ENCprime","B","SCUO","MCB","OGT"),lambda_trim,TRUE)
)
parameters <- list()
predictions <- list()
for (name in names(models)) {
  spec <- models[[name]]
  dat <- if (spec[[3]]) extremo else training
  exponent <- spec[[2]]
  dat$transformed <- if (exponent == 0) log(dat$d) else (dat$d^exponent-1)/exponent
  fitted <- lm(reformulate(spec[[1]],response="transformed"),data=dat)
  covariance <- vcov(fitted)
  values <- c(lambda=exponent,sigma=summary(fitted)$sigma,df=fitted$df.residual,coef(fitted))
  parameters[[name]] <- rbind(
    data.frame(model=name,kind="parameter",row=names(values),column="",value=unname(values)),
    data.frame(model=name,kind="covariance",row=rep(rownames(covariance),ncol(covariance)),
               column=rep(colnames(covariance),each=nrow(covariance)),value=as.vector(covariance))
  )
  prediction <- predict(fitted,newdata=training[c(3,17,39),],interval="confidence",level=.95)
  predictions[[name]] <- data.frame(model=name,observation=c(3,17,39),prediction)
}
write.csv(do.call(rbind,parameters),file.path(output,"fitted_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,predictions),file.path(output,"fitted_predictions.csv"),row.names=FALSE)
writeLines(c("Synthetic numeric parity fixtures; no genome sequences.",
             R.version.string, paste("MASS",as.character(packageVersion("MASS"))),
             "PCA uses upstream gRodon fitIoC and predictIoC verbatim.",
             "Regression uses MASS::boxcox and stats::lm/predict.lm matching gRodon formulas.",
             "Regenerate using scripts/generate_ioc_fitting_reference.R."),
           file.path(output,"provenance.txt"))
