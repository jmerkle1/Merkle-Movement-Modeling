### adapted from Mathieu's code for the new name of the lme4 models
### updated october 2013

kfoldRSF <- function(mod, k = 5, nrepet = 100, nbins = 10, jitter = FALSE,
                     random = TRUE, reproducible = TRUE)
{
  if (!inherits(mod, c("glm", "glmerMod")))
    stop("Model of class '(g)lm' or '(g)lmerMod' expected")
  if (inherits(mod, c("glmerMod")))
    require(lme4)
  dt <- model.frame(mod)
  kfold <- rd <- numeric(length = nrepet)
  resp <- as.character(attr(terms(mod), "variables"))[attr(terms(mod),
                                                           "response") + 1]
  for (i in 1:nrepet) {
    dt$sets <- "train"
    if(reproducible)
      set.seed(i)
    dt$sets[sample(which(dt[, resp] == 1), sum(dt[, resp] == 1)/k)] <- "test"
    reg <- update(mod, data = subset(dt, sets == "train"))
    if (inherits(mod, "glm"))
      predall <- exp(as.numeric(model.matrix(terms(reg),
                                             dt) %*% coef(reg)))
    else if (inherits(mod, "glmerMod"))
      predall <- exp(as.numeric(model.matrix(terms(reg),
                                             dt) %*% lme4::fixef(reg)))
    if (jitter) {
      if(reproducible)
        set.seed(i)
      predall <- jitter(predall)
    }
    quant <- quantile(predall[dt[, resp] == 0], probs = seq(from = 0,
                                                            to = 1, length.out = nbins + 1))
    quant[1] <- 0
    quant[length(quant)] <- Inf
    int <- factor(findInterval(predall[dt$sets == "test"],
                               quant), levels = 1:nbins)
    kfold[i] <- cor(1:nbins, table(int), method = "spearman")
    if (random) {
      if (reproducible)
        set.seed(i)
      dt$sets[sample(which(dt[, resp] == 0), sum(dt[, resp] ==
                                                   1)/k)] <- "rd"
      int <- factor(findInterval(predall[dt$sets == "rd"],
                                 quant), levels = 1:nbins)
      rd[i] <- cor(1:nbins, table(int), method = "spearman")
    }
  }
  if (random)
    return(data.frame(kfold = c(kfold, rd), type = rep(c("obs",
                                                         "rand"), each = nrepet)))
  else return(kfold)
}
