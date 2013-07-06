greaterOrEqual <- function(x,y) {
    precision <- sqrt(.Machine$double.eps)
    (x >= y) | (abs(x-y) <= precision)
}

closeSockets <- function() {
    allCon <- showConnections()
    socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
    sapply(socketCon, function(ii) close.connection(getConnection(ii)) )
}
