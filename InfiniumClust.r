InfiniumClust
function (tumor.data, purity, K = 3, maxiter = 100, tol = 0.001) 
{
    tumor.data = na.omit(as.matrix(tumor.data))
    vars = rowVars(tumor.data)
    names(vars) = rownames(tumor.data)
    dmp = names(sort(vars, decreasing = TRUE)[1:100])
    sample = intersect(colnames(tumor.data), names(purity))
    tumor.data = tumor.data[dmp, sample]
    purity = as.vector(purity[sample])
    par.intial = .initializeParameter(tumor.data, K)
    p = par.intial$p
    mu = par.intial$mu
    mu.n = par.intial$mu.n
    sd = par.intial$sd
    sd.n = par.intial$sd.n
    diff = 1
    iter = 0
    N_cpg = nrow(tumor.data)
    N_sample = ncol(tumor.data)
    Z = matrix(0, nrow = N_sample, ncol = K)
    while (diff > tol & iter < maxiter) {
        tol.ll2 = 0
        for (j in 1:nrow(Z)) {
            tmp1 <- rep(0, K)
            for (s in 1:K) {
                tmp1[s] <- log(p[s]) + sum(dnorm(tumor.data[, 
                  j], mean = purity[j] * mu[, s] + (1 - purity[j]) * 
                  mu.n, sd = sqrt(purity[j]^2 * sd[, s]^2 + (1 - 
                  purity[j])^2 * sd.n^2), log = T))
            }
            tmp1.sum.log <- .Rsumlog(tmp1)
            tol.ll2 <- tol.ll2 + tmp1.sum.log
            Z[j, ] = exp(tmp1 - tmp1.sum.log)
        }
        p.new = colSums(Z)/sum(Z)
        mu.n.new <- unlist(lapply(1:N_cpg, function(i) {
            tmp1 = 0
            tmp2 = 0
            for (j in 1:length(purity)) {
                for (k in 1:K) {
                  tmp1 = tmp1 + Z[j, k] * (1 - purity[j]) * (tumor.data[i, 
                    j] - purity[j] * mu[i, k])/(purity[j]^2 * 
                    sd[i, k]^2 + (1 - purity[j])^2 * sd.n[i]^2)
                  tmp2 = tmp2 + Z[j, k] * (1 - purity[j])^2/(purity[j]^2 * 
                    sd[i, k]^2 + (1 - purity[j])^2 * sd.n[i]^2)
                }
            }
            return(tmp1/tmp2)
        }))
        mu.new <- t(sapply(1:N_cpg, function(i) {
            mu.new.row = c()
            for (k in 1:ncol(Z)) {
                tmp1 <- sum(Z[, k] * purity * tumor.data[i, ]/(purity^2 * 
                  sd[i, k]^2 + (1 - purity)^2 * sd.n[i]^2)) - 
                  sum(Z[, k] * purity * (1 - purity)/(purity^2 * 
                    sd[i, k]^2 + (1 - purity)^2 * sd.n[i]^2)) * 
                    mu.n.new[i]
                tmp2 <- sum(Z[, k] * purity^2/(purity^2 * sd[i, 
                  k]^2 + (1 - purity)^2 * sd.n[i]^2))
                mu.new.row = append(mu.new.row, tmp1/tmp2)
            }
            return(mu.new.row)
        }))
        sd.new <- t(sapply(1:N_cpg, function(i) {
            sd.new.row = c()
            for (k in 1:K) {
                func_var <- function(x) {
                  sum(Z[, k] * (log(p[k]) + dnorm(tumor.data[i, 
                    ], mean = purity * mu[i, k] + (1 - purity) * 
                    mu.n[i], sd = sqrt(purity^2 * x + (1 - purity)^2 * 
                    sd.n[i]^2), log = T)))
                }
                temp = optimize(func_var, interval = c(0.001, 
                  100), maximum = T, tol = tol)
                sd.new.row = append(sd.new.row, sqrt(temp$maximum))
            }
            return(sd.new.row)
        }))
        sd.n.new <- unlist(lapply(1:N_cpg, function(i) {
            func_sd.n <- function(x) {
                obj <- 0
                for (j in 1:N_sample) {
                  for (k in 1:ncol(Z)) {
                    obj <- obj + Z[j, k] * (log(p[k]) + dnorm(tumor.data[i, 
                      j], mean = purity[j] * mu[i, k] + (1 - 
                      purity[j]) * mu.n[i], sd = sqrt(purity[j]^2 * 
                      sd[i, k]^2 + (1 - purity[j])^2 * x), log = T))
                  }
                }
                return(obj)
            }
            temp = optimize(func_sd.n, interval = c(0.001, 100), 
                maximum = T, tol = tol)
            return(sqrt(temp$maximum))
        }))
        diff <- sqrt(mean((mu.new - mu)^2 + (sd.new - sd)^2 + 
            (mu.n.new - mu.n)^2 + (sd.n.new - sd.n)^2))
        if (is.nan(diff)) 
            browser()
        p = p.new
        mu = mu.new
        mu.n = mu.n.new
        sd = sd.new
        sd.n = sd.n.new
        iter = iter + 1
        cat("Iter", iter, ", diff=", diff, "\n")
    }
    cat("# of clusters chosen:", K, "\n")
    cat("Probability of each cluster = ", p, "\n")
    cat("Total log-likelihood = ", tol.ll2, "\n")
    res = list(tol.ll = tol.ll2, Z = Z)
    return(res)
}
<environment: namespace:InfiniumPurify>




########################################################
myasin <- function(x) asin(2*x-1)

InfiniumPurify<-function (tumor.data, normal.data, purity) 
{
    if (missing(tumor.data) | missing(normal.data) | missing(purity)) {
        stop("'tumor.data', 'normal.data' and 'purity' are required.")
    }
    probes = intersect(rownames(tumor.data), rownames(normal.data))
    tumor.sample = intersect(colnames(tumor.data), names(purity))
    normal.sample = colnames(normal.data)
    purity = purity[tumor.sample]
    if (length(normal.sample) < 20 | length(tumor.sample) < 20) {
        stop("tumor and normal samples should be more than 20!")
    }
    .get_corrBeta <- function(input) {
        x = as.numeric(input[tumor.sample])
        y = as.numeric(input[normal.sample])
        type = c(rep("Tumor", length(x)), rep("Normal", length(y)))
        data = data.frame(beta = c(x, y), type = type)
        Y = .myasin(data$beta)
        X = c(1 - purity, rep(0, length(y)))
        fit = lm(Y ~ X)
        tmp = resid(fit) + coef(fit)[1]
        beta.pred = (sin(tmp) + 1)/2
        beta.pred
    }
    all.data = cbind(tumor.data[probes, tumor.sample], normal.data[probes, 
        ])
    probes.rmna = probes[rowSums(is.na(all.data)) == 0]
    all.data.corr = t(apply(all.data[probes.rmna, ], 1, .get_corrBeta))
    tumor.data.corr = all.data.corr[, 1:length(tumor.sample)]
    colnames(tumor.data.corr) = tumor.sample
    tumor.data.corr
}

<environment: namespace:InfiniumPurify>
