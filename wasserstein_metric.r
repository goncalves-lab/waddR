
#R Code:
wasserstein1d <- function (a, b, p = 1, wa = NULL, wb = NULL) {
    m <- length(a)
    n <- length(b)
    stopifnot(m > 0 && n > 0)
    if (m == n && is.null(wa) && is.null(wb)) {
        return( mean( abs( sort(b) - sort(a) )^p )^(1/p))
    }
    if (is.null(wa)) {
        wa <- rep(1, m)
    }
    if (is.null(wb)) {
        wb <- rep(1, n)
    }
    stopifnot(length(wa) == m && length(wb) == n)
    ua <- (wa/sum(wa))[-m]
    ub <- (wb/sum(wb))[-n]
    cua <- c(cumsum(ua))
    cub <- c(cumsum(ub))
    temp <- cut(cub, breaks = c(-Inf, cua, Inf))
    arep <- table(temp) + 1
    temp <- cut(cua, breaks = c(-Inf, cub, Inf))
    brep <- table(temp) + 1
    aa <- rep(sort(a), times = arep)
    bb <- rep(sort(b), times = brep)
    uu <- sort(c(cua, cub))
    uu0 <- c(0, uu)
    uu1 <- c(uu, 1)
    areap <- sum((uu1 - uu0) * abs(bb - aa)^p)^(1/p)
    return(areap)
}