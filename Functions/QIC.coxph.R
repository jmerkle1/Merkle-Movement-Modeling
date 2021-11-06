####mathieu's QIC function....

QIC.coxph <- function(mod, details = FALSE)
{
  if (!exists("naive.var", mod))
    stop("QIC can be computed only if robust variances are estimated.")
  trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
  quasi <- mod$loglik[2]
  if (details)
    return(data.frame(QICR = -2 * quasi + 2 * trace, QICI = -2 *
                        quasi + 2 * length(mod$coefficients), QuasiLL = quasi,
                      n = mod$n, nevent = mod$nevent, K = length(mod$coefficients),
                      Trace = trace))
  else return(-2 * quasi + 2 * trace)
}