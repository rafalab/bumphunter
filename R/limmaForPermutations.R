.getEstimate <- function(mat, design, coef, B=NULL, full=FALSE) {
    v <- design[,coef]
    A <- design[,-coef, drop=FALSE]
    qa <- qr(A)
    S <- diag(nrow(A)) - tcrossprod(qr.Q(qa))

    vv <- if (is.null(B)) matrix(v, ncol=1) else replicate(B, sample(v))
    sv <- S %*% vv
    vsv <- diag(crossprod(vv,sv))
    
    b <- (mat %*% crossprod(S, vv)) / vsv
    if(!is.matrix(b))
        b <- matrix(b, ncol = 1)
    if (full) {
        sy <- mat %*% S
        df.residual <- ncol(mat) - qa$rank - 1
        if (is.null(B)) {
            sigma <- matrix(sqrt(rowSums((sy - tcrossprod(b, sv))^2) / df.residual), ncol=1)
        } else {
            sigma <- b
            tmp <- sy
            for (j in 1:B) {
                tmp <- tcrossprod(b[,j], sv[,j])
                sigma[,j] <- rowSums((sy-tmp)^2)
            }
            sigma <- sqrt(sigma/df.residual)
        }
        out <- list(coef=b,
                    sigma=sigma,
                    stdev.unscaled=sqrt(1/vsv),
                    df.residual=df.residual)
        if (is.null(B)) out$stdev <- as.numeric(out$stdev)
    } else {
        out <- b
    }
    return(out)
}

.getModT <- function(obj) {
    s2 <- apply(obj$sigma^2, 2, limma::squeezeVar, obj$df.residual)
    out <- obj$coef
    for (j in 1:ncol(out)) {
        out[,j] <- out[,j] / obj$stdev.unscaled[j] / s2[[j]]$var.post
    }
    df.total <- obj$df.residual + sapply(s2,"[[","df.prior")
    return(list(t=out, df.total=df.total))
}
