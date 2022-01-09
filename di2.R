
# creating a function that pumps out maximum possible sample coverage without exceeding a double in sample size.
# We'll use this value (SC2) to determine if we should drop any data and to set a sample coverage point we extrapolate/rarefy to.

di2 = function (x, datatype = "abundance") {
  TYPE <- c("abundance", "incidence", "incidence_freq", 
            "incidence_raw")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  Fun.abun <- function(x) {
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x == k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x > 0)
    f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, 
                     (n - 1)/n * f1^2/2/f2)
    A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
    Chat <- round(1 - f1/n * A, 4)
    ME_CHANGE <- 1 - (1 - Chat) * exp(-2*f2/f1)
    ME_CHANGE <- ifelse(is.nan(ME_CHANGE), 1.0000000, ME_CHANGE)
    c(n, Sobs, Chat, ME_CHANGE, fk)
  }
  Fun.ince <- function(x) {
    nT <- x[1]
    x <- x[-1]
    U <- sum(x)
    Qk <- sapply(1:10, function(k) sum(x == k))
    Q1 <- Qk[1]
    Q2 <- Qk[2]
    Sobs <- sum(x > 0)
    Q0.hat <- ifelse(Q2 == 0, (nT - 1)/nT * Q1 * (Q1 - 1)/2, 
                     (nT - 1)/nT * Q1^2/2/Q2)
    A <- ifelse(Q1 > 0, nT * Q0.hat/(nT * Q0.hat + Q1), 1)
    Chat <- round(1 - Q1/U * A, 4)
    out <- c(nT, U, Sobs, Chat, Qk)
  }
  if (datatype == "abundance") {
    if (class(x) == "numeric" | class(x) == "integer") {
      out <- matrix(Fun.abun(x), nrow = 1)
    }
    else if (class(x) == "list") {
      out <- do.call("rbind", lapply(x, Fun.abun))
    }
    else if (class(x) == "matrix" | class(x) == "data.frame") {
      out <- t(apply(as.matrix(x), 2, Fun.abun))
    }
    if (nrow(out) > 1) {
      out <- data.frame(site = rownames(out), out)
      colnames(out) <- c("site", "n", "S.obs", 
                         "SC", "SC2", paste("f", 1:10, sep = ""))
      rownames(out) <- NULL
    }
    else {
      out <- data.frame(site = "site.1", out)
      colnames(out) <- c("site", "n", "S.obs", 
                         "SC", "SC2", paste("f", 1:10, sep = ""))
    }
    as.data.frame(out)
  }
  else if (datatype == "incidence") {
    if (class(x) == "numeric" | class(x) == "integer") {
      out <- matrix(Fun.ince(x), nrow = 1)
    }
    else if (class(x) == "list") {
      out <- do.call("rbind", lapply(x, Fun.ince))
    }
    else if (class(x) == "matrix" | class(x) == "data.frame") {
      out <- t(apply(as.matrix(x), 2, Fun.ince))
    }
    if (nrow(out) > 1) {
      out <- data.frame(site = rownames(out), out)
      colnames(out) <- c("site", "T", "U", 
                         "S.obs", "SC", paste("Q", 1:10, 
                                              sep = ""))
      rownames(out) <- NULL
    }
    else {
      out <- data.frame(site = "site.1", out)
      colnames(out) <- c("site", "T", "U", 
                         "S.obs", "SC", paste("Q", 1:10, 
                                              sep = ""))
    }
    as.data.frame(out)
  }
}
