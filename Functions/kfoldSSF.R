kfoldSSF <- function(mod, #fitted model
                     k = 5,   # how many folds to split the data into (e.g., 5 is analogous to 80% training, and 20% test)
                     nrepet = 100, # how many times to repeat the procedure
                     jitter = FALSE,  #should the predicted values be jittered? Helps deal with ties
                     reproducible = TRUE,   # sets set.seed so results are reproducible
                     details = FALSE,     # should details be spit out with the results?
                     sampleByCluster=FALSE)  #is the sampling based on strata (FALSE) or by cluster (TRUE)? 
{
  if (!inherits(mod, c("coxph", "coxme")))
    stop("Model of class 'coxph' or 'coxme' expected.")
  if (inherits(mod, "coxph"))
    require(survival)
  if (inherits(mod, "coxme"))
    stop("Not yet implemented for 'coxme' models.")
  dt <- try(model.frame(mod), silent = TRUE)
  if (class(dt) == "try-error")
    stop("'model.frame' was unable to retrieve the data.",
         "Use 'model = TRUE' in the 'coxph' call.")
  if(any(grepl("cluster",names(dt))) == FALSE & sampleByCluster==TRUE)
    stop("You must have a cluster term in your model to specify sampleByCluster=TRUE.")
  
  names(dt)[1] <- "srv"
  nstr <- attr(terms(mod), "specials")$strata
  names(dt)[nstr] <- namestr <- sub("strata\\((.*)\\)", "\\1", names(dt)[nstr])
  
  # get the name of the cluster column
  nameclu <- grepRaw("cluster", as.character(mod$userCall)[2])
  nameclu <- substr(as.character(mod$userCall)[2], nameclu, nchar(as.character(mod$userCall)[2]))
  nameclu <- sub("cluster\\((.*)\\)", "\\1", nameclu)
  names(dt)[names(dt)=="(cluster)"] <- nameclu
  
  # if (!is.null(attr(terms(mod), "specials")$cluster)) {
  #   nclu <- attr(terms(mod), "specials")$cluster
  #   names(dt)[nclu] <- nameclu <-  sub("cluster\\((.*)\\)", "\\1",
  #                          names(dt)[nclu])
  # }
  if(sampleByCluster==TRUE){
    if(k > length(unique(dt[, nameclu])))
      stop("You don't have enough clusters for your specification of k")
  }
  ncase <- table(tapply(dt$srv[, 2], dt[, nstr], function(x) sum(x == 1)))
  if (any(names(ncase) == "0"))
    stop(paste("Some stratas had no case.",
               "It is likely that NAs were present in the variables for some cases."))
  kfold <- rd <- warn <- numeric(length = nrepet)
  if (details)
    dbg <- list()
  for (i in 1:nrepet) {
    dt$sets <- "train"
    if (reproducible)
      set.seed(i)
    
    if(sampleByCluster==TRUE){
      dt$sets[dt[, nameclu] %in% sample(unique(dt[, nameclu]),
                                        length(unique(dt[, nameclu]))/k)] <- "test"
    }else{
      dt$sets[dt[, namestr] %in% sample(unique(dt[, namestr]),
                                        length(unique(dt[, namestr]))/k)] <- "test"
    }
    
    reg <- update(mod, srv ~ ., data = subset(dt, sets == "train"), model = TRUE)
    dtest <- droplevels(subset(dt, sets == "test"))
    if (inherits(mod, "coxph"))
      dtest$predall <- exp(predict(reg, type = "lp", newdata = dtest,
                                   reference = "sample"))
    if (jitter) {
      if (reproducible)
        set.seed(i)
      dtest$predall <- jitter(dtest$predall)
    }
    samplepred <- function(df) {
      nrand <- sum(df$srv[, 2] == 0)
      obs <- rank(df$predall)[df$srv[, 2] == 1]
      if (reproducible)
        set.seed(i)
      rand <- sample(rank(df$predall[df$srv[, 2] == 0]),1)
      return(data.frame(obs = obs, rand = rand, nrand = nrand))
    }
    ranks <- do.call(rbind, by(dtest, dtest[, namestr], samplepred))
    nrand <- unique(ranks$nrand)
    if (length(nrand) != 1) {
      nrand <- max(nrand)
      warn[i] <- 1
    }
    kfold[i] <- cor(1:(nrand+1), table(factor(ranks$obs,
                                              levels = 1:(nrand+1))), method = "spearman")
    rd[i] <- cor(1:(nrand), table(factor(ranks$rand, levels = 1:(nrand))), method = "spearman")
    if (details)
      dbg[[i]] <- ranks
  }
  res <- data.frame(kfold = c(kfold, rd), type = rep(c("obs","rand"), each = nrepet))
  if (details)
    attr(res, "details") <- dbg
  if (sum(warn) > 0)
    warning(paste("The number of controls was not equal among stratas for",
                  sum(warn), "repetitions. Correlations might be biased.",
                  "Use 'details = TRUE' to get more details."))
  return(res)
}