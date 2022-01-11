
#This code 
#A. sets working directory
#B. creates di2 function
#C. creates gginext2 function
#notes available for each function 


#### A. setting working directory ####
setwd("C:\\Users\\nerc-user\\OneDrive\\Documents\\GitHub\\HumDiv")

# setting seed so that I get the same result whenever I start with this seed
set.seed(12345)

# creating an "out" function cause it's useful
'%!in%' <- function(x,y)!('%in%'(x,y))

# formatting to html so that complex tables can be used
options(knitr.table.format = "html") 

#### B. di2 function ####
# creating a function that pumps out maximum possible sample coverage 
# without exceeding a double in sample size.
# We'll use this value (SC2) to determine if we should drop any data, 
# and to set a sample coverage point we extrapolate/rarefy to.

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

#### C. gginext function ####
# creating a ggiNEXT alternative that drops the legend and CI intervals 
# It's easier to see what's going on in overcrowded figures this way

gginext2 = function (x, type = 1, se = FALSE, facet.var = "none", color.var = "site", grey = FALSE) {
  TYPE <- c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if (is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1) 
    stop("invalid plot type")
  if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == 
      -1) 
    stop("invalid facet variable")
  if (is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == 
      -1) 
    stop("invalid color variable")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if (facet.var == "order") 
    color.var <- "site"
  if (facet.var == "site") 
    color.var <- "order"
  options(warn = -1)
  z <- fortify(x, type = type)
  options(warn = 0)
  if (ncol(z) == 7) {
    se <- FALSE
  }
  datatype <- unique(z$datatype)
  if (color.var == "none") {
    if (levels(factor(z$order)) > 1 & "site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep = "-")
    }
    else if ("site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }
    else if (levels(factor(z$order)) > 1) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }
    else {
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }
  else if (color.var == "order") {
    z$col <- z$shape <- factor(z$order)
  }
  else if (color.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }
  else if (color.var == "both") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep = "-")
  }
  zz = z
  z$method[z$method == "observed"] = "interpolated"
  z$lty <- z$lty <- factor(z$method, levels = unique(c("interpolated", 
                                                       "extrapolated"), c("interpolation", "interpolation", 
                                                                          "extrapolation")))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$method == "observed"), ]
  g <- ggplot(z, aes_string(x = "x", y = "y", colour = "col")) + 
    geom_point(aes_string(shape = "shape"), size = 5, 
               data = data.sub)
  g <- g + geom_line(aes_string(linetype = "lty"), lwd = 1.5) + 
    guides(linetype = guide_legend(title = "Method"), 
           colour = guide_legend(title = "Guides"), fill = guide_legend(title = "Guides"), 
           shape = guide_legend(title = "Guides")) + theme(legend.position = "bottom", 
                                                           legend.title = element_blank(), text = element_text(size = 18), 
                                                           legend.key.width = unit(1.2, "cm"))
  if (type == 2L) {
    g <- g + labs(x = "Number of sampling units", y = "Sample coverage")
    if (datatype == "abundance") 
      g <- g + labs(x = "Number of individuals", 
                    y = "Sample coverage")
  }
  else if (type == 3L) {
    g <- g + labs(x = "Sample coverage", y = "Species diversity")
  }
  else {
    g <- g + labs(x = "Number of sampling units", y = "Species diversity")
    if (datatype == "abundance") 
      g <- g + labs(x = "Number of individuals", 
                    y = "Species diversity")
  }
  if (facet.var == "order") {
    if (length(levels(factor(z$order))) == 1 & type != 2) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }
    else {
      g <- g + facet_wrap(~order, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides", 
                                              ncol = length(levels(factor(z$order))), byrow = TRUE), 
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }
    else {
      g <- g + facet_wrap(~site, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides", 
                                              nrow = length(levels(factor(z$order)))), fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "both") {
    if (length(levels(factor(z$order))) == 1 | !"site" %in% 
        names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }
    else {
      g <- g + facet_wrap(site ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides", 
                                              nrow = length(levels(factor(z$site))), byrow = TRUE), 
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (grey) {
    g <- g + theme_bw(base_size = 18) + scale_fill_grey(start = 0, 
                                                        end = 0.4) + scale_colour_grey(start = 0.2, end = 0.2) + 
      guides(linetype = guide_legend(title = "Method"), 
             colour = guide_legend(title = "Guides"), 
             fill = guide_legend(title = "Guides"), 
             shape = guide_legend(title = "Guides")) + 
      theme(legend.position = "bottom", legend.title = element_blank())
  }
  g <- g + theme(legend.box = "vertical") + theme(legend.position = "none")
  return(g)
}
