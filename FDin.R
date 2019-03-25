library(FD)

traitrix <- simulind(5, 5, 5, tr.method = "norm")
plot(traitrix)
write.csv(traitrix, "traitrix.csv")
distrix <- as.dist(traitrix)
FDin(distrix)

#compare to original fdisp function's results
fdisp(distrix)

simulind <- function (t = 3, r = 10, p = 100, tr.method = c("unif", "norm", "lnorm"))
{
	tr.method <- match.arg(tr.method)
	traitrix <- matrix(NA, p, t)
	if (tr.method == "unif") 
        traitrix <- apply(traitrix, 2, function(p) runif(p))
    if (tr.method == "norm") 
        traitrix <- apply(traitrix, 2, function(p) rnorm(p))
    if (tr.method == "lnorm") 
        traitrix <- apply(traitrix, 2, function(p) rlnorm(p))
	tr.numbs <- 1:t
	tr <- c("tr")
	tr.names <- paste(tr, tr.numbs, sep = "")
	ind.numbs <- 1:p
	ind <- c("ind")
	ind.names <- paste(ind, ind.numbs, sep = "")
	dimnames(traitrix) <- list(ind.names, tr.names)
	return(traitrix)
	
	#no.com <- 1:(r * l.s)
    #names.com <- paste("com", no.com, sep = "")
    #dimnames(abun) <- list(names.com, names.sp)
    values <- dbFD(traitrix, calc.FDiv = T, w.abun = F, 
        messages = F, calc.CWM = F)
    results <- cbind(nb.sp, values$FDin, values$FRic, values$FEve, 
        values$FDiv, values$RaoQ)
    names.var <- c("nb.sp", "FDin", "FRic", "FEve", "FDiv", "RaoQ")
    dimnames(results) <- list(names.com, names.var)
    #abun.gamma <- apply(abun, 2, mean)
    dbFD.gamma <- dbFD(traits, abun.gamma, calc.FDiv = F, w.abun = w.abun, 
        messages = F)
    FDis.gamma <- dbFD.gamma$FDis
    FDis.mean <- mean(values$FDis)
    nbsp.FDis <- cor.test(results[, 1], results[, 2])
    cor1 <- bquote(italic(r) == .(round(nbsp.FDis$estimate, 3)))
    FRic.FDis <- cor.test(results[, 2], results[, 3])
    cor2 <- bquote(italic(r) == .(round(FRic.FDis$estimate, 3)))
    FEve.FDis <- cor.test(results[, 2], results[, 4])
    cor3 <- bquote(italic(r) == .(round(FEve.FDis$estimate, 3)))
    FDiv.FDis <- cor.test(results[, 2], results[, 5])
    cor4 <- bquote(italic(r) == .(round(FDiv.FDis$estimate, 3)))
    RaoQ.FDis <- cor.test(results[, 6], results[, 2])
    cor6 <- bquote(italic(r) == .(round(RaoQ.FDis$estimate, 3)))
    nbsp.RaoQ <- cor.test(results[, 1], results[, 6])
    cor5 <- bquote(italic(r) == .(round(nbsp.RaoQ$estimate, 3)))
    par(mar = c(5, 5, 4, 2), mfrow = c(3, 2), las = 1)
    plot(results[, 3], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    title(xlab = "FRic", ylab = "FDis", cex.lab = 1.5)
    mtext("a", line = 1, adj = 0, cex = 1.5, font = 2)
    mtext(cor2, side = 3, line = 1, cex = 1.2)
    plot(results[, 5], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    mtext("b", line = 1, adj = 0, cex = 1.5, font = 2)
    title(xlab = "FDiv", ylab = "FDis", cex.lab = 1.5)
    mtext(cor4, side = 3, line = 1, cex = 1.2)
    plot(results[, 4], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    title(xlab = "FEve", ylab = "FDis", cex.lab = 1.5)
    mtext("c", line = 1, adj = 0, cex = 1.5, font = 2)
    mtext(cor3, side = 3, line = 1, cex = 1.2)
    plot(results[, 6], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    title(xlab = "Rao's Q", ylab = "FDis", cex.lab = 1.5)
    mtext("d", line = 1, adj = 0, cex = 1.5, font = 2)
    mtext(cor6, side = 3, line = 1, cex = 1.2)
    boxplot(results[, 2] ~ results[, 1], cex.axis = 1.2)
    mtext("e", line = 1, adj = 0, cex = 1.5, font = 2)
    title(xlab = "Species richness", ylab = "FDis", cex.lab = 1.5)
    mtext(cor1, side = 3, line = 1, cex = 1.2)
    boxplot(results[, 6] ~ results[, 1], cex.axis = 1.2)
    mtext("f", line = 1, adj = 0, cex = 1.5, font = 2)
    title(xlab = "Species richness", ylab = "Rao's Q", cex.lab = 1.5)
    mtext(cor5, side = 3, line = 1, cex = 1.2)
    res <- list()
    res$results <- results
    res$traits <- traits
    res$abun <- abun
    res$abun.gamma <- abun.gamma
    res$FDis.gamma <- FDis.gamma
    res$FDis.mean <- FDis.mean
    invisible(res)

}

# calculate FDin
# this is adapted from the fdisp function in the FD package

FDin <- function (d, tol = 1e-07) 
{
#    if (!inherits(d, "dist")) 
#        stop("'d' must be a 'dist' object.")
    n <- attr(d, "Size")
#    if (is.null(attr(d, "Labels"))) 
#        stop("'d' must have labels.", "\n")
#    else sn.d <- attr(d, "Labels")
#    if (any(is.na(d))) 
#        stop("NA's in the distance matrix.", "\n")
    d <- as.matrix(d)
    com <- nrow(d)
    A <- matrix(0, ncol = n, nrow = n)
    A <- -0.5 * d^2
    A <- A + t(A)
    G <- bicenter.wt(A)
    e <- eigen(G, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    w0 <- eig[n]/eig[1]
    if (w0 > -tol) 
        r <- sum(eig > (eig[1] * tol))
    else r <- length(eig)
    vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
        r)
    pos <- eig > 0
    avg.dist.cent <- rep(NA, nrow(d))
    names(avg.dist.cent) <- row.names(d)
        for (i in 1:com) {
        pres <- which(d[i, ] > 0)
        nb.in <- nrow((unique(vec <- vectors[pres, , drop = F])))
        if (nb.in >= 2) {
            w <- d[i, pres]
            centroid <- apply(vec, 2, weighted.mean, w = w)
            dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
            dist.pos <- rowSums(dist.pos^2)
            if (any(!pos)) {
                dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
                dist.neg <- rowSums(dist.neg^2)
            }
            else dist.neg <- 0
            zij <- sqrt(abs(dist.pos - dist.neg))
            avg.dist.cent[i] <- weighted.mean(zij, w)
        }
        else avg.dist.cent[i] <- 0
        communityFDin <- mean(avg.dist.cent)
    }
    return(list(FDin = avg.dist.cent, eig = eig, vectors = vectors, communityFDin = communityFDin))
}

#almost intact fdisp function, with some modifications to allow unlabeled matrices
FDis <- function (d, a, tol = 1e-07) 
{
    if (!inherits(d, "dist")) 
        stop("'d' must be a 'dist' object.")
    n <- attr(d, "Size")
    #if (is.null(attr(d, "Labels"))) 
    #    stop("'d' must have labels.", "\n")
    #else sn.d <- attr(d, "Labels")
    #if (missing(a)) {
    #    ab.names <- list("Community1", sn.d)
    #    a <- matrix(1, 1, n, dimnames = ab.names)
    #}
    com <- nrow(a)
    if (!is.matrix(a)) 
        stop("'a' must be a matrix.")
    if (ncol(a) != n) 
        stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
    if (is.null(colnames(a))) 
        stop("'a' must have column names", "\n")
    else sn.a <- colnames(a)
    if (any(sn.d != sn.a)) 
        stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
            "\n")
    a[which(is.na(a))] <- 0
    abun.sum <- apply(a, 1, sum)
    if (any(abun.sum == 0)) 
        stop("At least one community has zero-sum abundances (no species).", 
            "\n")
    abun.sum2 <- apply(a, 2, sum)
    if (any(abun.sum2 == 0)) 
        stop("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
    if (any(is.na(d))) 
        stop("NA's in the distance matrix.", "\n")
    A <- matrix(0, ncol = n, nrow = n)
    A[row(A) > col(A)] <- -0.5 * d^2
    A <- A + t(A)
    G <- bicenter.wt(A)
    e <- eigen(G, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    w0 <- eig[n]/eig[1]
    if (w0 > -tol) 
        r <- sum(eig > (eig[1] * tol))
    else r <- length(eig)
    vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
        r)
    dimnames(vectors) <- list(colnames(a), NULL)
    pos <- eig > 0
    avg.dist.cent <- rep(NA, nrow(a))
    names(avg.dist.cent) <- row.names(a)
    for (i in 1:com) {
        pres <- which(a[i, ] > 0)
        nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
        if (nb.sp >= 2) {
            w <- a[i, pres]
            centroid <- apply(vec, 2, weighted.mean, w = w)
            dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
            dist.pos <- rowSums(dist.pos^2)
            if (any(!pos)) {
                dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
                dist.neg <- rowSums(dist.neg^2)
            }
            else dist.neg <- 0
            zij <- sqrt(abs(dist.pos - dist.neg))
            avg.dist.cent[i] <- weighted.mean(zij, w)
        }
        else avg.dist.cent[i] <- 0
    }
    return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}


# simulate a dataset of communities with trait values for individuals
# this is adapted from simul.dbFD from the FD package

simulind <- function (i = c(5, 10, 15, 20, 25, 30, 35, 40), t = 3, r = 10, 
	p = 100, tr.method = c("unif", "norm", "lnorm")) 
{
    if (p < max(i)) 
        stop("'p' must be greater than the maximum number of individuals in 'i'.")
    l.i <- length(i)
    r.rep <- rep(r, l.i)
    nb.in <- rep(i, r.rep)
    tr.method <- match.arg(tr.method)
    traits <- matrix(NA, p, t)
    if (tr.method == "unif") 
        traits <- apply(traits, 2, function(p) runif(p))
    if (tr.method == "norm") 
        traits <- apply(traits, 2, function(p) rnorm(p))
    if (tr.method == "lnorm") 
        traits <- apply(traits, 2, function(p) rlnorm(p))
    fill.abun <- function(x, y, z) {
        set <- sample(1:length(x), size = y)
        if (z == "lnorm") 
            x[set] <- rlnorm(length(set))
        if (z == "norm") 
            x[set] <- rnorm(length(set))
        if (z == "unif") 
            x[set] <- runif(length(set))
        return(x)
    }
    abun <- mapply(fill.abun, abun, nb.sp, MoreArgs = list(z = abun.method))
    abun <- t(abun)
    no.tr <- 1:t
    tr <- c("tr")
    names.tr <- paste(tr, no.tr, sep = "")
    no.sp <- 1:p
    sp <- c("sp")
    names.sp <- paste(sp, no.sp, sep = "")
    dimnames(traits) <- list(names.sp, names.tr)
    no.com <- 1:(r * l.s)
    names.com <- paste("com", no.com, sep = "")
    dimnames(abun) <- list(names.com, names.sp)
    values <- dbFD(traits, abun, calc.FDiv = T, w.abun = w.abun, 
        messages = F, calc.CWM = F)
    results <- cbind(nb.sp, values$FDis, values$FRic, values$FEve, 
        values$FDiv, values$RaoQ)
    names.var <- c("nb.sp", "FDis", "FRic", "FEve", "FDiv", "RaoQ")
    dimnames(results) <- list(names.com, names.var)
    abun.gamma <- apply(abun, 2, mean)
    dbFD.gamma <- dbFD(traits, abun.gamma, calc.FDiv = F, w.abun = w.abun, 
        messages = F)
    FDis.gamma <- dbFD.gamma$FDis
    FDis.mean <- mean(values$FDis)
    nbsp.FDis <- cor.test(results[, 1], results[, 2])
    cor1 <- bquote(italic(r) == .(round(nbsp.FDis$estimate, 3)))
    FRic.FDis <- cor.test(results[, 2], results[, 3])
    cor2 <- bquote(italic(r) == .(round(FRic.FDis$estimate, 3)))
    FEve.FDis <- cor.test(results[, 2], results[, 4])
    cor3 <- bquote(italic(r) == .(round(FEve.FDis$estimate, 3)))
    FDiv.FDis <- cor.test(results[, 2], results[, 5])
    cor4 <- bquote(italic(r) == .(round(FDiv.FDis$estimate, 3)))
    RaoQ.FDis <- cor.test(results[, 6], results[, 2])
    cor6 <- bquote(italic(r) == .(round(RaoQ.FDis$estimate, 3)))
    nbsp.RaoQ <- cor.test(results[, 1], results[, 6])
    cor5 <- bquote(italic(r) == .(round(nbsp.RaoQ$estimate, 3)))
    par(mar = c(5, 5, 4, 2), mfrow = c(3, 2), las = 1)
    plot(results[, 3], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    title(xlab = "FRic", ylab = "FDis", cex.lab = 1.5)
    mtext("a", line = 1, adj = 0, cex = 1.5, font = 2)
    mtext(cor2, side = 3, line = 1, cex = 1.2)
    plot(results[, 5], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    mtext("b", line = 1, adj = 0, cex = 1.5, font = 2)
    title(xlab = "FDiv", ylab = "FDis", cex.lab = 1.5)
    mtext(cor4, side = 3, line = 1, cex = 1.2)
    plot(results[, 4], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    title(xlab = "FEve", ylab = "FDis", cex.lab = 1.5)
    mtext("c", line = 1, adj = 0, cex = 1.5, font = 2)
    mtext(cor3, side = 3, line = 1, cex = 1.2)
    plot(results[, 6], results[, 2], xlab = "", ylab = "", pch = ".", 
        cex.axis = 1.2)
    title(xlab = "Rao's Q", ylab = "FDis", cex.lab = 1.5)
    mtext("d", line = 1, adj = 0, cex = 1.5, font = 2)
    mtext(cor6, side = 3, line = 1, cex = 1.2)
    boxplot(results[, 2] ~ results[, 1], cex.axis = 1.2)
    mtext("e", line = 1, adj = 0, cex = 1.5, font = 2)
    title(xlab = "Species richness", ylab = "FDis", cex.lab = 1.5)
    mtext(cor1, side = 3, line = 1, cex = 1.2)
    boxplot(results[, 6] ~ results[, 1], cex.axis = 1.2)
    mtext("f", line = 1, adj = 0, cex = 1.5, font = 2)
    title(xlab = "Species richness", ylab = "Rao's Q", cex.lab = 1.5)
    mtext(cor5, side = 3, line = 1, cex = 1.2)
    res <- list()
    res$results <- results
    res$traits <- traits
    res$abun <- abun
    res$abun.gamma <- abun.gamma
    res$FDis.gamma <- FDis.gamma
    res$FDis.mean <- FDis.mean
    invisible(res)
}