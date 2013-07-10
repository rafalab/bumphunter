greaterOrEqual <- function(x,y) {
    precision <- sqrt(.Machine$double.eps)
    (x >= y) | (abs(x-y) <= precision)
}

closeSockets <- function() {
    allCon <- showConnections()
    socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
    sapply(socketCon, function(ii) close.connection(getConnection(ii)) )
}

foreachCleanup <- function() {
    if (exists(".revoDoParCluster", where=doParallel:::.options) && !is.null(doParallel:::.revoDoParCluster)) {
        stopCluster(doParallel:::.revoDoParCluster)
        remove(".revoDoParCluster", where=doParallel:::.options)
    }
}
