



# function to minimize for MIRF or MTRF method
obj <- function(par, T, ni, ab, aj1T, bj1T, cj1T, met, itmp, wt, D = D, base) 
{
  A <- rep(1, T)
  A[-base] <- par[1:(T - 1)]
  B <- rep(0, T)
  B[-base] <- par[T:(2 * T - 2)]
  if (itmp == 1)
    A[1:T] <- 1
  bj1Ts <- bj1T
  for (t in 1:T) bj1Ts[, t] <- (bj1Ts[, t] * A[t] + B[t])
  aj1Ts <- aj1T
  for (t in 1:T) aj1Ts[, t] <- aj1Ts[, t] / A[t]
  bj <- rowMeans(bj1Ts, na.rm = TRUE)
  aj <- rowMeans(aj1Ts, na.rm = TRUE)

  f <- 0
  for (t in 1:T) 
  {
    bjts <- (bj - B[t]) / A[t]
    ajts <- aj * A[t]
    ajt <- aj1T[, t]
    bjt <- bj1T[, t]
    cjt <- cj1T[, t]
    Pt <- matrix(NA, length(ab), ni)
    for (i in 1:ni) Pt[, i] <- equateIRT:::irtp1(ab, diff = bjt[i], discr = ajt[i], guess = cjt[i], D = D)
    Pts <- matrix(NA, length(ab), ni)
    for (i in 1:ni) Pts[, i] <- equateIRT:::irtp1(ab, diff = bjts[i], discr = ajts[i], guess = cjt[i], D = D)
    Pts[is.na(Pt)] <- NA
    if (met == "irf")
      f <- f + 0.5 * sum((Pt - Pts)^2 * wt, na.rm = TRUE)
    if (met == "trf")
      f <- f + 0.5 * sum(((rowSums(Pt, na.rm = TRUE) - rowSums(Pts, na.rm = TRUE))^2) * wt)
  }
  return(f)
}

# numerical first derivative of obj with respect to the item parameters
# used for finding the numerical second derivative of obj
# with respect to the equating coefficients and the item parameters
Sgamma <- function(x, par, T, ni, ab, tabl, nummet, itmp, wt, D = D, base) 
{
  tabl$gamma <- x[rownames(tabl)]
  tab <- reshape(tabl, direction = "wide", v.names = "gamma", timevar = "time", idvar = "itms")
  bj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dffclt", ][, -1])
  if (itmp > 1)
    aj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dscrmn", ][, -1]) else aj1T <- matrix(1, ni, T)
  if (itmp == 3)
    cj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Gussng", ][, -1]) else cj1T <- matrix(0, ni, T)
  grad(func = objectivefzRcpp, x = par, T = T, ab = ab, wt = wt, aj1T = aj1T, bj1T = bj1T, cj1T = cj1T, nummet = nummet, itmp = itmp, D = D, base = base)
}



# this functions returns the matrix of second derivatives of obj
# (at the minimum) with respect to the equating coefficients and 
# the item parameters
# calls Rcpp
pABgammaR2C <- function(par, T, ab, wt, aj1T, bj1T, cj1T, nummet, itmp, D, base) 
{
  nb <- sum(!is.na(bj1T))
  posnomi <- bj1T
  posnomi[!is.na(posnomi)] <- 1:nb
  tabnomia <- outer(rownames(aj1T), 1:T, paste, sep = ".")
  tabnomia[is.na(aj1T)] <- NA
  tabnomib <- outer(rownames(bj1T), 1:T, paste, sep = ".")
  tabnomib[is.na(bj1T)] <- NA
  if (itmp == 3) 
  {
    tabnomic <- outer(rownames(cj1T), 1:T, paste, sep = ".")
    tabnomic[is.na(cj1T)] <- NA
  }
  nomia <- tabnomia[!is.na(tabnomia)]
  nomib <- tabnomib[!is.na(tabnomib)]
  if (itmp == 3)
    nomic <- tabnomic[!is.na(tabnomic)]

  tt <- partialABgammaRcpp(par, T, ab, wt, aj1T, bj1T, cj1T, nummet, itmp, D, base, nb, posnomi)

  out_a <- tt[[1]]
  out_b <- tt[[2]]
  if (itmp == 3)
    out_c <- tt[[3]] else out_c <- NULL
  colnames(out_a) <- nomia
  colnames(out_b) <- nomib
  if (itmp == 3)
    colnames(out_c) <- nomic

  out <- cbind(out_c, out_a, out_b)
  out <- out[-c(base, base + T), ]
  return(out)
}



# this function is used to compute the numerical derivatives of a_j^* and b_j^* 
# with respect to the item parameter estimates
# as input requires the item parameter estimates (gamma)
# as output returns a_j^* and b_j^* 
gamma2ab <- function(x, par, T, ni, ab, tabl, nummet, itmp, wt, D = D, base, ini) 
{
  tabl$gamma <- x[rownames(tabl)]
  tab <- reshape(tabl, direction = "wide", v.names = "gamma", timevar = "time", idvar = "itms")
  bj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dffclt", ][, -1])
  if (itmp > 1)
    aj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dscrmn", ][, -1]) else aj1T <- matrix(1, ni, T)
  if (itmp == 3)
    cj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Gussng", ][, -1]) else cj1T <- matrix(0, ni, T)
  out <- nlminb(start = ini, objective = objectivefzRcpp, T = T, ab = ab, wt = wt, aj1T = aj1T, bj1T = bj1T, cj1T = cj1T, nummet = nummet, itmp = itmp, D = 1, base = base, control = list(iter.max = 1e+05, eval.max = 10000, trace = FALSE))
  par <- out$par
  As <- rep(1, T)
  As[-base] <- par[1:(T - 1)]
  Bs <- rep(0, T)
  Bs[-base] <- par[T:(2 * T - 2)]
  bj1Ts <- bj1T
  for (t in 1:T) bj1Ts[, t] <- bj1Ts[, t] * As[t] + Bs[t]
  aj1Ts <- aj1T
  for (t in 1:T) aj1Ts[, t] <- aj1Ts[, t] / As[t]
  bj <- rowMeans(bj1Ts, na.rm = TRUE)
  aj <- rowMeans(aj1Ts, na.rm = TRUE)
  return(c(aj, bj))
}



multiec_irf <- function(mods, base, nq = 30, method, se, start, iter.max, trace) 
{
  if (trace)
    cat("Computation of equating coefficients  .  .  .  . \n")
  if (inherits(start, "mlteqc"))
    start <- c(start$A[-1], start$B[-1])
  itms_all <- lapply(mods, function(x) names(x$coef)) # item labels
  itms <- unique(unlist(itms_all))
  itms <- sort(itms)
  T <- length(mods)
  modsnames <- names(mods)
  tab <- data.frame(itms = itms)
  for (k in 1:length(mods)) tab <- merge(tab, mods[[k]]$coef, by.x = "itms", by.y = 0, all = TRUE, suffixes = c(k - 1, k))
  colnames(tab)[-1] <- modsnames
  tab$itms <- as.character(tab$itms)
  rownames(tab) <- tab$itms
  sel_administered2plus <- rowSums(!is.na(tab[, -1])) >= 2
  items_administered2plus <- tab$itms[sel_administered2plus]
  tab_restr <- tab[sel_administered2plus, ]

  itmp <- 2
  if (sum(substr(tab$itms, 1, 6) == "Dscrmn") == 0)
    itmp = 1
  if (sum(substr(tab$itms, 1, 6) == "Gussng") > 0)
    itmp = 3

  if (is.null(start))
    ini <- c(rep(1, T - 1), rep(0, T - 1))
  if (!is.null(start))
    ini <- start
  names(ini) <- c(modsnames[-base], modsnames[-base])

  gq <- gauss.quad.prob(nq, dist = "normal")
  ab <- gq$nodes
  wt <- gq$weights

  ni <- nrow(tab) / itmp
  ni_restr <- nrow(tab_restr) / itmp
  bj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dffclt", ][, -1])
  if (itmp > 1)
    aj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dscrmn", ][, -1])
  if (itmp == 1) 
  {
    aj1T <- matrix(1, ni, T)
    aj1T[is.na(bj1T)] <- NA
  }
  if (itmp == 3)
    cj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Gussng", ][, -1])
  if (itmp < 3)
    cj1T <- matrix(0, ni, T)
  if (method == "irf")
    nummet <- 1
  if (method == "trf")
    nummet <- 2
  bj1T_restr <- bj1T[rownames(bj1T) %in% items_administered2plus, ]
  if (itmp > 1)
    aj1T_restr <- aj1T[rownames(aj1T) %in% items_administered2plus, ]
  if (itmp == 1) 
  {
    aj1T_restr <- matrix(1, ni_restr, T)
    aj1T[is.na(bj1T)] <- NA
  }
  if (itmp == 3)
    cj1T_restr <- cj1T[rownames(cj1T) %in% items_administered2plus, ]
  if (itmp < 3)
    cj1T_restr <- matrix(0, ni_restr, T)
  # the following uses Rcpp function
  out <- nlminb(start = ini, objective = objectivefzRcpp, gradient = gradRcpp, hessian = hessRcpp, T = T, ab = ab, wt = wt, aj1T = aj1T_restr, bj1T = bj1T_restr, cj1T = cj1T_restr, nummet = nummet, itmp = itmp, D = 1, base = base, control = list(iter.max = iter.max))
  par <- out$par
  As <- rep(1, T)
  names(As) <- modsnames
  if (itmp > 1)
    As[-base] <- par[1:(T - 1)]
  Bs <- rep(0, T)
  names(Bs) <- modsnames
  Bs[-base] <- par[T:(2 * T - 2)]

  bj1Ts <- bj1T
  for (t in 1:T) bj1Ts[, t] <- bj1Ts[, t] * As[t] + Bs[t]
  aj1Ts <- aj1T
  for (t in 1:T) aj1Ts[, t] <- aj1Ts[, t] / As[t]
  bj <- rowMeans(bj1Ts, na.rm = TRUE)
  aj <- rowMeans(aj1Ts, na.rm = TRUE)
  names(bj) <- tab$itms[substr(tab$itms, 1, 6) == "Dffclt"]
  names(aj) <- tab$itms[substr(tab$itms, 1, 6) == "Dscrmn"]
  rownames(aj1T) <- tab$itms[substr(tab$itms, 1, 6) == "Dscrmn"]
  rownames(bj1T) <- tab$itms[substr(tab$itms, 1, 6) == "Dffclt"]
  rownames(cj1T) <- tab$itms[substr(tab$itms, 1, 6) == "Gussng"]

  if (se) 
  {
    if (trace)
      cat("Computation of standard errors ")
    # h1 is the matrix of second derivatives at the minimum (hessian)
    # the following uses Rcpp function
    h1 <- hessRcpp(par = out$par, T = T, ab = ab, wt = wt, aj1T = aj1T_restr, bj1T = bj1T_restr, cj1T = cj1T_restr, nummet = nummet, itmp = itmp, D = 1, base = base)
    colnames(h1) <- c(rep("A", T - 1), rep("B", T - 1))
    if (trace) cat(" . ")

    # pfABg is the matrix of second derivatives with respect to the equating coefficients and the item parameters
    # the following uses Rcpp function
    pfABg <- pABgammaR2C(par = out$par, T = T, ab = ab, wt = wt, aj1T = aj1T_restr, bj1T = bj1T_restr, cj1T = cj1T_restr, nummet = nummet, itmp = itmp, D = 1, base = base)
    if (itmp == 1) 
    {
      pfABg[1:(T - 1), ] <- 0
      pfABg <- pfABg[, grep("Dffclt", colnames(pfABg))]
    }
    if (trace) cat(" . ")

    if (itmp == 1) 
    {
      h1 <- h1[-(1:(T - 1)), -(1:(T - 1))]
      pfABg <- pfABg[-(1:(T - 1)), ]
    }
    # dAB_gamma <- -solve(h1) %*% pfABg
    dAB_gamma <- -solve(h1, pfABg)
    dAB_gamma <- dAB_gamma[, colSums(dAB_gamma) != 0]
    if (trace) cat(" . ")
    sel <- colnames(dAB_gamma)
    vars <- lapply(mods, FUN = function(x) x$var)
    var.names <- list()
    for (i in 1:length(vars)) 
    {
      var.names_i <- paste(rownames(vars[[i]]), i, sep = ".")
      var.names[[i]] <- var.names_i
      rownames(vars[[i]])<-colnames(vars[[i]]) <- var.names_i
    }
    VarAll <- bdiag(vars)
    VarAllNames <- unlist(var.names)
    rownames(VarAll) <- colnames(VarAll) <- VarAllNames
    VarAll_tdAB_gamma <- Matrix::tcrossprod(VarAll[sel, sel], dAB_gamma)
    varAB <- dAB_gamma %*% VarAll_tdAB_gamma
    seAB <- Matrix::diag(varAB)^0.5
    seA <- rep(0, T)
    names(seA) <- modsnames
    if (itmp > 1)
      seA[-base] <- seAB[1:(T - 1)]
    seB <- rep(0, T)
    names(seB) <- modsnames
    if (itmp > 1)
      seB[-base] <- seAB[T:(2 * T - 2)]
    if (itmp == 1)
      seB[-base] <- seAB[1:(T - 1)]

    # standard errors of item parameters on a common scale (a_j* and b_j*)

    invT <- as.matrix(rowSums(!is.na(bj1T)))  # -> u_j
    invT <- invT[, c(rep(1, T))]
    colnames(invT) <- modsnames

    if (itmp > 1) 
    {
      pajA <- as.matrix(aj1T)
      for (t in 1:T) pajA[, t] <- aj1T[, t] / As[t]^2
      pajA[is.na(pajA)] <- 0
      pajA <- -pajA / invT
      pajA <- as.matrix(pajA[, -base])

      tablong <- reshape(tab, direction = "long", varying = list(2:6), idvar = "itms", v.names = "value")
      tablong <- tablong[!is.na(tablong$value), ]
      itms_t <- rownames(tablong)
      dAB_gamma_all <- matrix(0, (T - 1) * 2, length(itms_t))
      colnames(dAB_gamma_all) <- itms_t
      rownames(dAB_gamma_all) <- rownames(dAB_gamma)
      dAB_gamma_all[, colnames(dAB_gamma)] <- dAB_gamma

      grAa <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "A"), select = grepl("^Dscrmn", colnames(dAB_gamma_all)))  # derivatives of equating coefficients A with respect to the discrimination item parameters
      pajgamma_a <- pajA %*% grAa  # this is the result of the sum in Equations (91) and (92) multiplied by -1 / u_j

      colnames_spl <- strsplit(colnames(grAa), split = ".", fixed = TRUE)
      colnames_no_t <- sapply(colnames_spl, FUN = function(x) paste(x[1], x[2], sep = "."))
      sel_sameitem <- outer(rownames(aj1T), colnames_no_t, FUN = "==")

      whicht <- as.numeric(as.character(sapply(colnames_spl, FUN = function(x) x[3])))

      pajgamma_a[sel_sameitem] <- pajgamma_a[sel_sameitem] + (matrix(1 / invT[, 1]) %*% (1 / As[whicht]))[sel_sameitem]  # these are the derivatives in Equations (91) and (92)

      grAb <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "A"), select = grepl("^Dffclt", colnames(dAB_gamma_all)))  # derivatives of equating coefficients A with respect to the difficulty item parameters
      pajgamma_b <- pajA %*% grAb  # this is the result of the sum in Equation (93) multiplied by -1 / u_j

      grAc <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "A"), select = grepl("^Gussng", colnames(dAB_gamma_all)))  # derivatives of equating coefficients A with respect to the difficulty item parameters
      if (itmp == 3)
        pajgamma_c <- pajA %*% grAc  # similar to Equation (92), derivative of a_j^* with respect to c_it
      else pajgamma_c <- c()

      pajgamma <- cbind(pajgamma_a, pajgamma_b, pajgamma_c)

      sel <- colnames(pajgamma)
      VarAll_tpajgamma <- Matrix::tcrossprod(VarAll[sel, sel], pajgamma)
      seaj <- Matrix::diag(pajgamma %*% VarAll_tpajgamma)^0.5
      seaj <- seaj[names(aj)]
    }
    if (itmp == 1)
    {
      tablong <- reshape(tab, direction = "long", varying = list(2:6), idvar = "itms", v.names = "value")
      tablong <- tablong[!is.na(tablong$value), ]
      itms_t <- rownames(tablong)
      dAB_gamma_all <- matrix(0, T - 1, length(itms_t))
      colnames(dAB_gamma_all) <- itms_t
      rownames(dAB_gamma_all) <- rownames(dAB_gamma)
      dAB_gamma_all[, colnames(dAB_gamma)] <- dAB_gamma
      aj <- seaj <- NULL
    }

    if (trace) cat(" . \n")

    pbjA <- as.matrix(bj1T)
    pbjA[is.na(pbjA)] <- 0
    pbjA <- pbjA / invT
    pbjA <- as.matrix(pbjA[, -base])

    pbjB <- pbjA
    pbjB[pbjB != 0] <- 1 / invT[, -base][pbjB != 0]

    grAa <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "A"), select = grepl("^Dscrmn", colnames(dAB_gamma_all)))  # derivatives of equating coefficients A with respect to the discrimination item parameters
    grBa <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "B"), select = grepl("^Dscrmn", colnames(dAB_gamma_all)))  # derivatives of equating coefficients B with respect to the discrimination item parameters
    if (itmp > 1)
      pbjgamma_a <- pbjA %*% grAa + pbjB %*% grBa  # Equation (94)
    if (itmp == 1)
      pbjgamma_a <- pbjB %*% grBa  # Equation (94)

    grAb <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "A"), select = grepl("^Dffclt", colnames(dAB_gamma_all)))  # derivatives of equating coefficients A with respect to the difficulty item parameters
    grBb <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "B"), select = grepl("^Dffclt", colnames(dAB_gamma_all)))  # derivatives of equating coefficients B with respect to the difficulty item parameters
    if (itmp > 1)
      pbjgamma_b <- pbjA %*% grAb + pbjB %*% grBb  # this is the result of the sum in Equations (95) and (96) multiplied by -1 / u_j
    if (itmp == 1)
      pbjgamma_b <- pbjB %*% grBb  # this is the result of the sum in Equations (95) and (96) multiplied by -1 / u_j

    grAc <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "A"), select = grepl("^Gussng", colnames(dAB_gamma_all)))  # derivatives of equating coefficients A with respect to the guessing item parameters
    grBc <- subset(dAB_gamma_all, subset = (rownames(dAB_gamma_all) == "B"), select = grepl("^Gussng", colnames(dAB_gamma_all)))  # derivatives of equating coefficients B with respect to the guessing item parameters
    if (itmp == 3)
      pbjgamma_c <- pbjA %*% grAc + pbjB %*% grBc  # similar to Equation (94), derivative of b_j^* with respect to c_it
    else pbjgamma_c <- c()

    colnames_spl <- strsplit(colnames(grBb), split = ".", fixed = TRUE)
    colnames_no_t <- sapply(colnames_spl, FUN = function(x) paste(x[1], x[2], sep = "."))
    sel_sameitem <- outer(rownames(bj1T_restr), colnames_no_t, FUN = "==")

    whicht <- as.numeric(as.character(sapply(colnames_spl, FUN = function(x) x[3])))

    pbjgamma_b[sel_sameitem] <- pbjgamma_b[sel_sameitem] + (matrix(1 / invT[, 1]) %*% As[whicht])[sel_sameitem]

    pbjgamma <- cbind(pbjgamma_a, pbjgamma_b, pbjgamma_c)

    sel <- colnames(pbjgamma)
    
    VarAll_tpbjgamma <- Matrix::tcrossprod(VarAll[sel, sel], pbjgamma)
    sebj <- Matrix::diag(pbjgamma %*% VarAll_tpbjgamma)^0.5
    
    sebj <- sebj[names(bj)]

    partial <- t(dAB_gamma)  # derivatives of A and B equating coefficients with respect to the item parameters
    if (itmp == 1)
      partial <- cbind(matrix(0, nrow(partial), T - 1), partial)
    namesAB <- c(paste("A", (1:T)[-base], sep = "."), paste("B", (1:T)[-base], sep = "."))
    colnames(partial) <- namesAB
    if (itmp == 1) 
    {
      varAB1 <- matrix(0, 2 * T - 2, 2 * T - 2)
      varAB1[T:(2 * T - 2), T:(2 * T - 2)] <- as.matrix(varAB)
      varAB <- varAB1
    }
    rownames(varAB) <- colnames(varAB) <- namesAB
  } 
  else 
  {
    seA <- rep(NA, T)
    seB <- rep(NA, T)
    varAB <- matrix(NA, T, T)
    seaj <- NULL
    sebj <- NULL
    vars <- NULL
    partial <- NULL
  }
  if (itmp == 1)
    tabs <- bj1Ts
  if (itmp == 2)
    tabs <- rbind(bj1Ts, aj1Ts)
  if (itmp == 3)
    tabs <- rbind(bj1Ts, aj1Ts, cj1T)
  colnames(tabs) <- paste(colnames(tabs), modsnames[base], sep = ".as.")
  tabs <- tabs[, -base]
  tab <- cbind(tab, tabs)
  colnames(tab)[1] <- "Item"
  out <- list(A = As, B = Bs, se.A = seA, se.B = seB, varAB = varAB, as = aj, bs = bj, se.as = seaj, se.bs = sebj, tab = tab, varFull = vars, partial = partial, itmp = itmp, method = method, basename = modsnames[base], convergence = out$convergence)
  class(out) <- "mlteqc"
  return(out)
}






ipf4der <- function(gamma, itms, t, base, aj1T) 
{
  names(gamma) <- itms
  for (tt in 1:max(t)) aj1T[, paste("gamma", tt, sep = ".")] <- gamma[t == tt][rownames(aj1T)]
  mmout <- ipfRcpp(as.matrix(aj1T[, -1]), base, 1e-05)
  As <- mmout[[1]]
  return(As)
}



# this function is used to calculate numerical derivatives
# of a*, A, b*, B with respect to the discrimination and difficulty parameters
der_asAbsB_ab <- function(gamma, tab, base, method1, T, ni) 
{
  tab$gamma <- gamma
  itms <- tab$itempar

  tabDscrmn <- tab[substr(tab$itempar, 1, 6) == "Dscrmn", ]  # table with discrimination parameters
  tabDffclt <- tab[substr(tab$itempar, 1, 6) == "Dffclt", ]  # table with difficulty parameters

  # MULTIPLE MEAN-GEOMETRIC MEAN (performed anyway)
  # step 1: estimation of A equating coefficients
  if (nrow(tabDscrmn) > 0) 
  {
    if (any(tabDscrmn$gamma < 0))
      warning("Negative discrimination parameter. Converted to positive value.")
    tabDscrmn$gamma <- abs(tabDscrmn$gamma)  # negative discriminations converted to positive values
    tabDscrmn$t <- as.factor(tabDscrmn$t)
    tabDscrmn$t <- relevel(tabDscrmn$t, ref = base)  #set the A equating coefficient of the reference administration to 1
    reg1 <- lm(log(gamma) ~ factor(itempar) + t - 1, data = tabDscrmn, x = TRUE)  # estimation of the A equating coefficients with the multiple mean-geometric mean method (Haberman, 2009)
    As <- rep(1, T)
    As[-base] <- exp(reg1$coef)[(ni + 1):(ni + T - 1)]  # store A equating coefficients
    aj <- exp(reg1$coef[1:ni])  # store discrimination parameters on a common scale
    names(aj) <- substr(names(aj), 16, 1000)
  }
  # MULTIPLE MEAN-MEAN
  # step 1: estimation of A equating coefficients
  if (nrow(tabDscrmn) > 0 & method1 == "mean-mean") 
  {
    aj1T <- reshape(tabDscrmn[, c("itempar", "gamma", "t")], direction = "wide", v.names = "gamma", timevar = "t", idvar = "itempar")  # discrimination parameters in wide format
    rownames(aj1T) <- aj1T$itempar
    mmout <- ipfRcpp(as.matrix(aj1T[, -1]), base, 1e-07)
    As <- mmout[[1]]
    aj_tmp <- mmout[[2]]
    names(aj_tmp) <- aj1T$itempar
    aj <- aj_tmp[names(aj)]
  }

  # if there are no discrimination parameters set all A equating coefficients to 1
  if (nrow(tabDscrmn) == 0) 
  {
    As <- rep(1, T)
    aj <- rep(1, ni)
    names(aj) <- names(seaj) <- paste("Dscrmn", substr(itms, 8, 100), sep = ".")
  }

  # step 2: estimation of B equating coefficients (for both multiple mean-mean and multiple mean geometric mean methods)
  vettA <- As[tabDffclt$t]
  tabDffclt$gammaA <- tabDffclt$gamma * vettA  # response variable for the second regression model
  tabDffclt$t <- as.factor(tabDffclt$t)
  tabDffclt$t <- relevel(tabDffclt$t, ref = base)  #set the B equating coefficient of the reference administration to 0
  reg2 <- lm(gammaA ~ factor(itempar) + t - 1, data = tabDffclt, x = TRUE)  # estimation of the B equating coefficients
  Bs <- rep(0, T)
  Bs[-base] <- -reg2$coef[(ni + 1):(ni + T - 1)]  # store B equating coefficients
  bj <- reg2$coef[1:ni]  # store difficulty parameters on a common scale
  names(bj) <- substr(names(bj), 16, 1000)

  return(c(aj, As[-base], bj, Bs[-base]))
}


# multiple equating coefficients
multiec <- function(mods, base = 1, method = "mean-mean", se = TRUE, nq = 30, start = NULL, iter.max = 100000, obsinf = TRUE, trace = TRUE) 
{
  stopifnot(
    "'mods' must an output of function 'modIRT'" = inherits(mods, "modIRT"),
    "'base' must always be an integer value between 1 and the length of mods" = round(base) - base == 0 & base >= 1 & base <= length(mods),
    "'method' must be one of 'mean-mean', 'mean-gmean', 'irf', 'trf' or 'lik'"  = method %in% c("mean-mean", "mean-gmean", "irf", "trf", "lik"),
    "'se' must be logical" = is.logical(se),
    "'nq' must be an integer value between 3 and 100" = (round(nq) - nq == 0 & nq >= 3 & nq <= 100),
    "'start' must be an output of function 'multiec' or a vector with dimension 2*T-2, where T is the number of forms" = (is.null(start) | inherits(start, "mlteqc") | (is.vector(start) & length(start) == length(mods) * 2 - 2)),
    "'iter.max' must be numeric value" = is.numeric(iter.max),
    "'obsinf' must be logical" = is.logical(obsinf),
    "'trace' must be logical" = is.logical(trace)
  )
  
  if (method == "mean-mean" | method == "mean-gmean")
    out <- multiec_moments(mods = mods, base = base, method = method, se = se, trace = trace)
  if (method == "irf" | method == "trf")
    out <- multiec_irf(mods = mods, base = base, method = method, se = se, nq = nq, start = start, iter.max = iter.max, trace = trace)
  if (method == "lik")
    out <- multiec_lik(mods = mods, base = base, se = se, start = start, iter.max = iter.max, obsinf = obsinf, trace = trace)
  return(out)
}


# multiple equating coefficients with methods based on moments
multiec_moments <- function(mods, base, method, se, trace) 
{
  if (trace)
    cat("Computation of equating coefficients  .  .  .  . \n")
  itms_all <- lapply(mods, function(x) names(x$coef)) # item labels
  itms <- unique(unlist(itms_all))
  itms <- sort(itms)
  itms <- itms[substr(itms, 1, 6) != "Gussng"]
  T <- length(mods)  # number of administrations
  modsnames <- names(mods)
  # tab: table with all item parameters. Columns: itempar (label of item parameters), gamma (item parameter value), t (administration)
  tab1 <- list()
  for (k in 1:length(mods))
    tab1[[k]] <- data.frame(itempar = names(mods[[k]]$coef), gamma = mods[[k]]$coef, t = k, stringsAsFactors = FALSE)
  tab <- do.call("rbind", tab1)
  rownames(tab) <- paste(tab$itempar, tab$t, sep = ".")
  tabDscrmn <- tab[substr(tab$itempar, 1, 6) == "Dscrmn", ]  # table with discrimination parameters
  tabDffclt <- tab[substr(tab$itempar, 1, 6) == "Dffclt", ]  # table with difficulty parameters
  ni <- length(grep("Dffclt", itms))  # numer of items

  # MULTIPLE MEAN-GEOMETRIC MEAN (performed anyway)
  # step 1: estimation of A equating coefficients
  if (nrow(tabDscrmn) > 0) 
  {
    if (any(tabDscrmn$gamma < 0))
      warning("Negative discrimination parameter. Converted to positive value.")
    tabDscrmn$gamma <- abs(tabDscrmn$gamma)  # negative discriminations converted to positive values
    tabDscrmn$t <- as.factor(tabDscrmn$t)
    tabDscrmn$t <- relevel(tabDscrmn$t, ref = base)  #set the A equating coefficient of the reference administration to 1
    reg1 <- lm(log(gamma) ~ factor(itempar) + t - 1, data = tabDscrmn, x = TRUE)  # estimation of the A equating coefficients with the multiple mean-geometric mean method (Haberman, 2009)
    As <- rep(1, T)
    names(As) <- modsnames
    As[-base] <- exp(reg1$coef)[(ni + 1):(ni + T - 1)]  # store A equating coefficients
    aj <- exp(reg1$coef[1:ni])  # store discrimination parameters on a common scale
    names(aj) <- substr(names(aj), 16, 1000)
  }

  # MULTIPLE MEAN-MEAN
  # step 1: estimation of A equating coefficients
  if (nrow(tabDscrmn) > 0 & method == "mean-mean") 
  {
    aj1T <- reshape(tabDscrmn[, c("itempar", "gamma", "t")], direction = "wide", v.names = "gamma", timevar = "t", idvar = "itempar")  # discrimination parameters in wide format
    rownames(aj1T) <- aj1T$itempar
    mmout <- ipfRcpp(as.matrix(aj1T[, -1]), base, 1e-07)
    As <- mmout[[1]]
    aj_tmp <- mmout[[2]]
    names(As) <- modsnames
    names(aj_tmp) <- aj1T$itempar
    aj <- aj_tmp[names(aj)]
  }

  # if there are no discrimination parameters set all A equating coefficients to 1 and standard errors to 0
  if (nrow(tabDscrmn) == 0) 
  {
    As <- rep(1, T)
    names(As) <- modsnames
    seA <- rep(0, T)
    names(seA) <- modsnames
    aj <- rep(1, ni)
    seaj <- rep(0, ni)
    names(aj) <- names(seaj) <- paste("Dscrmn", substr(itms, 8, 100), sep = ".")
  }

  # step 2: estimation of B equating coefficients (for both multiple mean-mean and multiple mean geometric mean methods)
  vettA <- As[tabDffclt$t]
  tabDffclt$gammaA <- tabDffclt$gamma * vettA  # response variable for the second regression model
  tabDffclt$t <- as.factor(tabDffclt$t)
  tabDffclt$t <- relevel(tabDffclt$t, ref = base)  #set the B equating coefficient of the reference administration to 0
  reg2 <- lm(gammaA ~ factor(itempar) + t - 1, data = tabDffclt, x = TRUE)  # estimation of the B equating coefficients
  Bs <- rep(0, T)
  names(Bs) <- modsnames
  Bs[-base] <- -reg2$coef[(ni + 1):(ni + T - 1)]  # store B equating coefficients
  bj <- reg2$coef[1:ni]  # store difficulty parameters on a common scale
  names(bj) <- substr(names(bj), 16, 1000)

  if (se) 
  {
    # ================================
    # computation of standard errors
    # ================================

    if (trace)
      cat("Computation of standard errors ")

    tabDscrmn$t <- as.numeric(as.character(tabDscrmn$t))

    if (nrow(tabDscrmn) > 0) 
    {
      vars <- lapply(mods, FUN = function(x) x$var)
      var.names <- list()
      for (i in 1:length(vars)) 
      {
        var.names_i <- paste(rownames(vars[[i]]), i, sep = ".")
        var.names[[i]] <- var.names_i
        rownames(vars[[i]])<-colnames(vars[[i]]) <- var.names_i
      }
      VarAll <- bdiag(vars)
      VarAllNames <- unlist(var.names)
      rownames(VarAll) <- colnames(VarAll) <- VarAllNames

      X1 <- reg1$x  # design matrix
      X2 <- reg2$x  # design matrix
      X2[, (ni + 1):(ni + T - 1)] <- -X2[, (ni + 1):(ni + T - 1)]
      colnames(X1)[1:ni] <- substr(colnames(X1)[1:ni], 16, 1000)
      colnames(X2)[1:ni] <- substr(colnames(X2)[1:ni], 16, 1000)
      sel1 <- rownames(tabDscrmn)
      sel2 <- rownames(tabDffclt)
      sel <- c(sel1, sel2)

      if (trace) cat(" . ")

      if (method == "mean-gmean") 
      {
        invXX1 <- summary(reg1)$cov.unscaled
        tmp2 <- matD(t(X1), 1 / tabDscrmn$gamma)
        tmp1 <- Dmat(exp(reg1$coef), invXX1)
        pAa <- tmp1 %*% tmp2
        rownames(pAa) <- rownames(t(X1))
        colnames(pAa) <- sel1

        if (trace) cat(" . ")

      }
      if (method == "mean-mean") 
      {
        pAa0 <- jacobian(func = ipf4der, x = tabDscrmn$gamma, itms = tabDscrmn$itempar, t = tabDscrmn$t, base = base, aj1T = aj1T)

        if (trace) cat(" . ")

        rownames(pAa0) <- modsnames
        colnames(pAa0) <- rownames(tabDscrmn)
        sumajs <- tapply(tabDscrmn$gamma, tabDscrmn$itempar, sum)
        sumAs <- tapply(As[tabDscrmn$t], tabDscrmn$itempar, sum)
        sumpAa0 <- aggregate(pAa0[tabDscrmn$t, ], by = list(Group.1 = tabDscrmn$itempar), FUN = sum)
        rn <- sumpAa0$Group.1
        sumpAa0 <- sumpAa0[, -1]
        tmp1<-apply(sumpAa0, 2, sumajs / sumAs^2, FUN="*")
        tmp2 <- as.matrix(1 / sumAs)
        tmp2 <- tmp2[, rep(1, ncol(tmp1))]
        tmp2 <- tmp2 * t(X1[, 1:ni])
        colnames(tmp2) <- colnames(tmp1)
        pAa1 <- tmp2 - tmp1
        pAa <- rbind(pAa1, pAa0[-base, ])  # see Battauz (2016) Equations (41) and (42)

      }

      if (trace) cat(" . ")

      pAb <- matrix(0, ncol(X1), nrow(X1))
      colnames(pAb) <- sel2
      tmp <- Dmat(tabDffclt$gamma, -X2[, (ni + 1):(ni + T - 1)])
      pBa2 <- t(X2) %*% tmp
      invXX2 <- summary(reg2)$cov.unscaled
      colnames(invXX2)[1:ni] <- rownames(invXX2)[1:ni] <- colnames(X2)[1:ni]
      invXX2[1:ni, -(1:ni)] <- -invXX2[1:ni, -(1:ni)]
      invXX2[-(1:ni), 1:ni] <- -invXX2[-(1:ni), 1:ni]
      pBa <- (invXX2 %*% (pBa2)) %*% pAa[(ni + 1):(ni + T - 1), ]  # see Battauz (2016) Equation (39)
      invXX2_tX2 <- Matrix::tcrossprod(invXX2, X2)
      pBb <- matD(invXX2_tX2, vettA)  # see Battauz (2016) Equation (40)
      tmp1 <- cbind(pAa, pAb)
      tmp2 <- cbind(pBa, pBb)
      colnames(tmp2) <- colnames(tmp1)
      part <- rbind(tmp1, tmp2)
      if (trace) cat(" . \n")
      VarAll_tpart <- Matrix::tcrossprod(VarAll[sel, sel], part)
      varABgamma <- part %*% VarAll_tpart # covariance matrix of equating coefficients and item parameters on a common scale
      seABgamma <- Matrix::diag(varABgamma)^0.5
      seA <- rep(0, T)
      names(seA) <- modsnames
      seA[-base] <- seABgamma[(ni + 1):(ni + T - 1)]  # standard errors of A equating coefficients
      seB <- rep(0, T)
      names(seB) <- modsnames
      seB[-base] <- seABgamma[(ni + ncol(X1) + 1):(ni + ncol(X1) + T - 1)]  # standard errors of B equating coefficients
      varAB <- varABgamma[c((ni + 1):(ni + T - 1), (ni + ncol(X1) + 1):(ni + ncol(X1) + T - 1)), c((ni + 1):(ni + T - 1), (ni + ncol(X1) + 1):(ni + ncol(X1) + T - 1))]  # covariance matrix of equating coefficients
      seaj <- seABgamma[1:ni]  # standard errors of discrimination parameters aj on a common scale
      sebj <- seABgamma[(ni + T):(2 * ni + T - 1)]  # standard errors of difficulty parameters bj on a common scale

      partial <- t(part[c((ni + 1):(ni + T - 1), (ni + ncol(X1) + 1):(ni + ncol(X1) + T - 1)), ])  # derivatives of A and B equating coefficients with respect to the item parameters
      namesAB <- c(paste("A", (1:T)[-base], sep = "."), paste("B", (1:T)[-base], sep = "."))
      colnames(partial) <- namesAB
      rownames(varAB) <- colnames(varAB) <- namesAB
    }

    if (nrow(tabDscrmn) == 0) 
    {
      X2 <- reg2$x
      mat <- Matrix::tcrossprod(summary(reg2)$cov.unscaled, X2)
      if (trace) cat(" . ")
      colnames(mat) <- sel <- paste(tabDffclt$itempar, tabDffclt$t, sep = ".")
      vars <- lapply(mods, FUN = function(x) x$var)
      var.names <- list()
      for (i in 1:length(vars)) 
      {
        var.names_i <- paste(rownames(vars[[i]]), i, sep = ".")
        var.names[[i]] <- var.names_i
        rownames(vars[[i]])<-colnames(vars[[i]]) <- var.names_i
      }
      VarAll <- bdiag(vars)
      VarAllNames <- unlist(var.names)
      rownames(VarAll) <- colnames(VarAll) <- VarAllNames
      if (trace)  cat(" . ")
      sel <- paste(tabDffclt$itempar, tabDffclt$t, sep = ".")
      VarAll_tmat <- Matrix::tcrossprod(VarAll[sel, sel], mat)
      varB <- mat %*% VarAll_tmat
      if (trace) cat(" . ")
      seBt <- sqrt(Matrix::diag(varB))
      seB <- rep(0, T)
      names(seB) <- modsnames
      seB[-base] <- seBt[(ni + 1):(ni + T - 1)]
      sebj <- seBt[1:ni]
      names(sebj) <- substr(names(sebj), 16, 1000)
      if (trace) cat(" . \n")
      partial <- matrix(0, nrow(tab), (T - 1) * 2)
      rownames(partial) <- paste(tab$itempar, tab$t, sep = ".")
      tmp5 <- t(mat[(ni + 1):(ni + T - 1), ])
      partial[rownames(tmp5), T:(2 * T - 2)] <- tmp5
      namesAB <- c(paste("A", (1:T)[-base], sep = "."), paste("B", (1:T)[-base], sep = "."))
      colnames(partial) <- namesAB
      varAB <- matrix(0, (T - 1) * 2, (T - 1) * 2)
      varAB[T:(T * 2 - 2), T:(T * 2 - 2)] <- as.matrix(varB[(ni + 1):(ni + T - 1), (ni + 1):(ni + T - 1)])
      rownames(varAB) <- colnames(varAB) <- namesAB
    }
  } 
  else
  {
    seA <- rep(NA, T)
    seB <- rep(NA, T)
    varAB <- matrix(NA, T, T)
    seaj <- NULL
    sebj <- NULL
    vars <- NULL
    partial <- NULL
  }

  tabwide <- reshape(tab[, c("itempar", "gamma", "t")], direction = "wide", v.names = "gamma", timevar = "t", idvar = "itempar")  # parameters in wide format
  colnames(tabwide) <- c("Item", modsnames)
  tabwide <- tabwide[order(tabwide$Item), ]
  tabwides <- tabwide[, -1]
  for (t in 1:T) 
  {
    selDffclt <- grep("Dffclt", tabwide$Item)
    selDscrmn <- grep("Dscrmn", tabwide$Item)
    tabwides[selDffclt, t] <- tabwides[selDffclt, t] * As[t] + Bs[t]  # conversion of difficulty parameters to the scale of the base form
    tabwides[selDscrmn, t] <- tabwides[selDscrmn, t] / As[t]  # conversion of discrimination parameters to the scale of the base form
  }
  colnames(tabwides) <- paste(colnames(tabwides), modsnames[base], sep = ".as.")
  tabwides <- tabwides[, -base]
  tabwide <- cbind(tabwide, tabwides)
  itmp <- 2
  if (sum(substr(tabwide$Item, 1, 6) == "Dscrmn") == 0)
    itmp = 1
  if (sum(substr(tabwide$Item, 1, 6) == "Gussng") > 0)
    itmp = 3

  out <- list(A = As, B = Bs, se.A = seA, se.B = seB, varAB = varAB, as = aj, bs = bj, se.as = seaj, se.bs = sebj, tab = tabwide, varFull = vars, partial = partial, itmp = itmp, method = method, basename = modsnames[base])
  class(out) <- "mlteqc"
  return(out)
}


print.mlteqc <- function(x, ...) 
{
  cat("Multiple equating coefficients \n")
  cat("Method: ")
  cat(x$method, "\n")
}


summary.mlteqc <- function(object, ...) 
{
  ct <- data.frame(EQ = rep(c("A", "B"), each = length(object$A)), Form = c(names(object$A), names(object$B)), Estimate = c(object$A, object$B), StdErr = c(object$se.A, object$se.B))
  cat("Equating coefficients:\n")
  print(ct, digits = 5, row.names = FALSE)
}




# product of a matrix mat and a diagonal matrix with diagonal D: mat%*%diag(D)
matD <- function(mat, D) 
{
  nr <- nrow(mat)
  for (i in 1:nr) mat[i, ] <- mat[i, ] * D
  rownames(mat) <- colnames(mat) <- NULL
  return(mat)
}

# product of a diagonal matrix with diagonal D and a matrix mat: diag(D)%*%mat
Dmat <- function(D, mat) 
{
  nc <- ncol(mat)
  for (i in 1:nc) mat[, i] <- D * mat[, i]
  rownames(mat) <- colnames(mat) <- NULL
  return(mat)
}





itm.mlteqc <- function(x, ...) x$tab


eqc.mlteqc <- function(x, ...) 
{
  A1 <- x$A
  B1 <- x$B
  out <- data.frame(A = A1, B = B1)
  rownames(out) <- NULL
  return(out)
}

item.common <- function(x, ...) UseMethod("item.common")

item.common.mlteqc <- function(x, ...)
{
  as <- x$as
  bs <- x$bs
  if (!is.null(x$se.as)) se.as <- x$se.as
  else se.as <- NA
  if (!is.null(x$se.bs)) se.bs <- x$se.bs
  else se.bs <- NA
  as.df <- data.frame(Item = names(as), Estimate = as, StdErr = se.as)
  bs.df <- data.frame(Item = names(bs), Estimate = bs, StdErr = se.bs)
  out <- rbind(as.df, bs.df)
  rownames(out) <- NULL
  out
}


score.mlteqc <- function(obj, method = "TSE", D = 1, scores = NULL, se = TRUE, nq = 30, w = 0.5, theta = NULL, weights = NULL, ...) 
{
  stopifnot(
    "obj must be an output of fucntion 'multiec'" = inherits(obj, "mlteqc"),
    "'method' must be 'TSE' or 'OSE'"  = method %in% c("TSE", "OSE"),
    "'D' must be a numeric value" = is.numeric(D),
    "'scores' should include integer values" = ifelse(is.null(scores), TRUE, !any(round(scores) != scores)), 
    "'se' must be logical" = is.logical(se),
    "'nq' must be an integer value between 3 and 100" = (round(nq) - nq == 0 & nq >= 3 & nq <= 100),
    "'w' must be a numeric value between 0 and 1" = is.null(w) | (is.numeric(w) & w >= 0 & w <= 1),
    "'theta' must be a numeric vector with length > 1" = ifelse(is.null(theta), TRUE, is.numeric(theta) & length(theta) > 1),
    "'weights' must be a numeric vector with the same length of theta" = ifelse(is.null(weights), TRUE, is.numeric(weights) & length(weights) == length(theta))
  )
  
  itmpar <- itm(obj)
  basename <- obj$basename
  modsnames <- names(obj$A)
  T <- length(obj$A)
  base <- (1:T)[modsnames == basename]
  out <- NULL
  printmessages <- TRUE
  count <- 0
  for (t in 1:T) 
  {
    if (base != t) 
    {
      count <- count + 1
      if (count == 2) printmessages <- FALSE
      sel <- c("Item", modsnames[t], basename, paste(modsnames[t], basename, sep = ".as."))
      sel1 <- c(paste("A", t, sep = "."), paste("B", t, sep = "."))
      itm_prepare <- equateIRT:::score_prepare(itmpar[, sel], suff1 = paste(".", t, sep = ""), suff2 = paste(".", base, sep = ""))
      out_t <- equateIRT:::score_compute(method = method, itm_prepare = itm_prepare, D = D, varFull = obj$varFull, partial = obj$partial[, sel1], varAB = obj$varAB[sel1, sel1], itmp = obj$itmp, A = obj$A[t], B = obj$B[t], scores = scores, se = se, nq = nq, w = w, theta = theta, weights = weights, names = sel[3:4], printmessages = printmessages)
      if (is.null(out)) 
      {
        k <- ncol(out_t)
        if (se)
          colnames(out_t)[k] <- paste(colnames(out_t)[k], colnames(out_t)[k - 1], sep = "_")
        out <- out_t
      } 
      else 
      {
        if (method == "TSE")
          out_t <- subset(out_t, select = 3:ncol(out_t))
        if (method == "OSE")
          out_t <- subset(out_t, select = 2:ncol(out_t))
        if (se)
          colnames(out_t)[2] <- paste(colnames(out_t)[2], colnames(out_t)[1], sep = "_")
        out <- cbind(out, out_t)
      }
    }
  }
  return(out)
}






# ProfLik<-function(par,itmpar,itmvar,num.forms,base,DffcltNum,DscrmnNum)
# {
#   coefs<-NULL
#   A<-rep(1,num.forms)
#   A[-base]<-par[1:(num.forms-1)]
#   B<-rep(0,num.forms)
#   B[-base]<-par[((num.forms-1)+1):(2*(num.forms-1))]
#   itmpar$A<-A[itmpar$t]
#   itmpar$B<-B[itmpar$t]
#   itmpar$Y<-9999
#   itmpar[DscrmnNum,'Y']<-itmpar[DscrmnNum,coef / A]
#   itmpar[DffcltNum,'Y']<-itmpar[DffcltNum,coef*A+B]
#   Y<-itmpar$Y
#   X<-model.matrix(~factor(items)-1,data=itmpar)
#   row.names(X)<-itmpar$items.t
#   X_list<-list()
#   for (f in 1:num.forms) X_list[[f]]<-X[itmpar$t==f,]
#   Y_list<-list()
#   for (f in 1:num.forms) Y_list[[f]]<-Y[itmpar$t==f]
#   omega<-itmvar
#   for (t in 1:num.forms)
#   {
#     selds<-grep('Dscrmn',rownames(omega[[t]]))
#     seldf<-grep('Dffclt',rownames(omega[[t]]))
#     omega[[t]][selds,selds]<-omega[[t]][selds,selds] / A[t]^2
#     omega[[t]][seldf,seldf]<-omega[[t]][seldf,seldf]*A[t]^2
#   }
#   V.chol <- try(lapply(omega, Matrix::chol))
#   V.chol.inv <- try(lapply(V.chol, Matrix::solve))
#   if (!isa(V.chol.inv,'try-error'))
#   {
#     t.V.chol.inv.X<-mapply(FUN=Matrix::crossprod,V.chol.inv,X_list,SIMPLIFY=F)
#     t.V.chol.inv.X<-do.call('rbind', t.V.chol.inv.X)
#     tX.V.chol.inv_t.V.chol.inv.X<-Matrix::crossprod(t.V.chol.inv.X)
#     # max(abs(tX.V.chol.inv_t.V.chol.inv.X-tX.Oinv.X))
#     t.V.chol.inv.Y<-mapply(FUN=Matrix::crossprod,V.chol.inv,Y_list,SIMPLIFY=F)
#     t.V.chol.inv.Y<-do.call('rbind', t.V.chol.inv.Y)
#     tX.V.chol.inv_t.V.chol.inv.Y<-Matrix::crossprod(t.V.chol.inv.X,t.V.chol.inv.Y)
#     beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X,tX.V.chol.inv_t.V.chol.inv.Y)
#     # max(abs(beta-beta1))
#     
#     rownames(beta)<-substr(colnames(X),14,100)
#     itmpar$coefs<-beta[itmpar$items,]
#     itmpar$mean<-9999
#     itmpar[DscrmnNum,mean:=coefs*A]
#     itmpar[DffcltNum,mean:=(coefs-B) / A]
#     logden<-itmpar[,dmvnorm(x=coef,mean=mean,sigma=itmvar[[t]],log=TRUE),by=t]
#     out<-sum(logden$V1)
#   }
#   else out<- 1e100000000
#   -out
# }
# 
# 
# 
# ProfLik_1PL<-function(par,itmpar,itmvar,num.forms,base)
# {
#   coefs<-NULL
#   B<-rep(0,num.forms)
#   B[-base]<-par
#   itmpar$B<-B[itmpar$t]
#   itmpar$Y<-9999
#   itmpar$Y<-itmpar[,coef+B]
#   Y<-itmpar$Y
#   X<-model.matrix(~factor(items)-1,data=itmpar)
#   X_list<-list()
#   for (f in 1:num.forms) X_list[[f]]<-X[itmpar$t==f,]
#   Y_list<-list()
#   for (f in 1:num.forms) Y_list[[f]]<-Y[itmpar$t==f]
#   omega<-itmvar
#   V.chol <- try(lapply(omega, Matrix::chol))
#   V.chol.inv <- try(lapply(V.chol, Matrix::solve))
#   if (!isa(V.chol.inv,'try-error'))
#   {
#     t.V.chol.inv.X<-mapply(FUN=Matrix::crossprod,V.chol.inv,X_list,SIMPLIFY=F)
#     t.V.chol.inv.X<-do.call('rbind', t.V.chol.inv.X)
#     tX.V.chol.inv_t.V.chol.inv.X<-Matrix::crossprod(t.V.chol.inv.X)
#     # max(abs(tX.V.chol.inv_t.V.chol.inv.X-tX.Oinv.X))
#     t.V.chol.inv.Y<-mapply(FUN=Matrix::crossprod,V.chol.inv,Y_list,SIMPLIFY=F)
#     t.V.chol.inv.Y<-do.call('rbind', t.V.chol.inv.Y)
#     tX.V.chol.inv_t.V.chol.inv.Y<-Matrix::crossprod(t.V.chol.inv.X,t.V.chol.inv.Y)
#     beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X,tX.V.chol.inv_t.V.chol.inv.Y)
#     # max(abs(beta-beta1))
#     
#     rownames(beta)<-substr(colnames(X),14,100)
#     itmpar$coefs<-beta[itmpar$items,]
#     itmpar$mean<-9999
#     itmpar[,mean:=(coefs-B)]
#     logden<-itmpar[,dmvnorm(x=coef,mean=mean,sigma=itmvar[[t]],log=TRUE),by=t]
#     out<-sum(logden$V1)
#   }
#   else out<- 1e100000000
#   -out
# }
# 



# estimates of the item parameters on a common scale using the estimates from all the forms
getBetas <- function(par, itmpar, itmvar, num.forms, base, DffcltNum, DscrmnNum, X_list) 
{
  A <- rep(1, num.forms)
  A[-base] <- par[1:(num.forms - 1)]
  B <- rep(0, num.forms)
  B[-base] <- par[((num.forms - 1) + 1):(2 * (num.forms - 1))]
  itmpar$A <- A[itmpar$t]
  itmpar$B <- B[itmpar$t]
  itmpar$Y <- 9999
  itmpar[DscrmnNum, "Y"] <- itmpar[DscrmnNum, coef / A]
  itmpar[DffcltNum, "Y"] <- itmpar[DffcltNum, coef * A + B]
  Y <- itmpar$Y
  Y_list <- list()
  for (f in 1:num.forms) Y_list[[f]] <- Y[itmpar$t == f]

  omega <- itmvar
  for (t in 1:num.forms) 
  {
    selds <- grep("Dscrmn", rownames(omega[[t]]))
    seldf <- grep("Dffclt", rownames(omega[[t]]))
    omega[[t]][selds, selds] <- omega[[t]][selds, selds] / A[t]^2
    omega[[t]][seldf, seldf] <- omega[[t]][seldf, seldf] * A[t]^2
  }
  V.chol <- try(lapply(omega, Matrix::chol))
  V.chol.inv <- try(lapply(V.chol, Matrix::solve))
  V.inv <- lapply(V.chol.inv, Matrix::tcrossprod)

  if (!isa(V.chol.inv, "try-error")) 
  {
    t.V.chol.inv.X <- mapply(FUN = Matrix::crossprod, V.chol.inv, X_list, SIMPLIFY = FALSE)
    t.V.chol.inv.X <- do.call("rbind", t.V.chol.inv.X)

    # the following is necessary to compute standard errors of as and bs
    V.inv.X <- mapply(FUN = Matrix::crossprod, V.inv, X_list, SIMPLIFY = FALSE)
    V.inv.X <- do.call("rbind", V.inv.X)

    tX.V.chol.inv_t.V.chol.inv.X <- Matrix::crossprod(t.V.chol.inv.X)
    t.V.chol.inv.Y <- mapply(FUN = Matrix::crossprod, V.chol.inv, Y_list, SIMPLIFY = FALSE)
    t.V.chol.inv.Y <- do.call("rbind", t.V.chol.inv.Y)
    tX.V.chol.inv_t.V.chol.inv.Y <- Matrix::crossprod(t.V.chol.inv.X, t.V.chol.inv.Y)
    beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X, tX.V.chol.inv_t.V.chol.inv.Y)
  } else beta <- NA
  list(beta = beta, tX.Vinv.X = tX.V.chol.inv_t.V.chol.inv.X, V.inv.X = V.inv.X)
}



getBetas_1PL <- function(par, itmpar, itmvar, num.forms, base, X_list) 
{
  B <- rep(0, num.forms)
  B[-base] <- par
  itmpar$B <- B[itmpar$t]
  itmpar$Y <- 9999
  itmpar$Y <- itmpar[, coef + B]
  Y <- itmpar$Y
  Y_list <- list()
  for (f in 1:num.forms) Y_list[[f]] <- Y[itmpar$t == f]

  omega <- itmvar
  V.chol <- try(lapply(omega, Matrix::chol))
  V.chol.inv <- try(lapply(V.chol, Matrix::solve))
  V.inv <- lapply(V.chol.inv, Matrix::tcrossprod)

  if (!isa(V.chol.inv, "try-error")) 
  {
    t.V.chol.inv.X <- mapply(FUN = Matrix::crossprod, V.chol.inv, X_list, SIMPLIFY = FALSE)
    t.V.chol.inv.X <- do.call("rbind", t.V.chol.inv.X)

    # the following is necessary to compute standard errors of as and bs
    V.inv.X <- mapply(FUN = Matrix::crossprod, V.inv, X_list, SIMPLIFY = FALSE)
    V.inv.X <- do.call("rbind", V.inv.X)

    tX.V.chol.inv_t.V.chol.inv.X <- Matrix::crossprod(t.V.chol.inv.X)
    t.V.chol.inv.Y <- mapply(FUN = Matrix::crossprod, V.chol.inv, Y_list, SIMPLIFY = FALSE)
    t.V.chol.inv.Y <- do.call("rbind", t.V.chol.inv.Y)
    tX.V.chol.inv_t.V.chol.inv.Y <- Matrix::crossprod(t.V.chol.inv.X, t.V.chol.inv.Y)
    beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X, tX.V.chol.inv_t.V.chol.inv.Y)
  } else beta <- NA
  list(beta = beta, tX.Vinv.X = tX.V.chol.inv_t.V.chol.inv.X, V.inv.X = V.inv.X)

}





getitmpar <- function(mods, t) 
{
  ct <- mods$coefficients
  itemstype <- substr(names(ct), 1, 6)
  noGussng <- itemstype != "Gussng"
  ct <- ct[noGussng]
  itemstype <- itemstype[noGussng]
  itemslab <- substr(names(ct), 8, 100)
  data.frame(items = names(ct), itemstype = itemstype, itemslab = itemslab, coef = ct, t = t)
}

delGussng <- function(x) 
{
  which_guess <- grep("Gussng", rownames(x))
  if (length(which_guess) > 0)
    x <- x[-which_guess, -which_guess]
  x
}

multiec_lik <- function(mods, base, se, obsinf, start, iter.max, trace) 
{
  itemstype <- items <- Y <- items.t <- NULL
  if (trace) cat("Computation of equating coefficients  .  .  . ")
  num.forms <- length(mods)

  itmpar <- mapply(FUN = getitmpar, mods = mods, t = 1:num.forms, SIMPLIFY = FALSE)
  itmpar <- rbindlist(itmpar)
  itmpar[, items.t := paste(items, t, sep = ".")]

  DscrmnNum <- grep("Dscrmn", itmpar$items)
  DffcltNum <- grep("Dffclt", itmpar$items)

  itmvar <- lapply(mods, FUN = function(x) x$var)
  itmvar <- lapply(itmvar, delGussng)
  for (i in 1:length(itmvar)) rownames(itmvar[[i]]) <- colnames(itmvar[[i]]) <- paste(rownames(itmvar[[i]]), i, sep = ".")

  itmp <- 2
  if (sum(substr(itmpar$items, 1, 6) == "Dscrmn") == 0)
    itmp = 1
  if (sum(substr(itmpar$itmes, 1, 6) == "Gussng") > 0)
    itmp = 3

  if (itmp >= 2) 
  {
    if (inherits(start, "mlteqc"))
      start <- c(start$A[-base], start$B[-base])
    if (is.null(start)) 
    {
      capture.output(irf <- multiec(mods, se = FALSE, method = "irf", base = base))
      start <- c(irf$A[-base], irf$B[-base])
    }
  }
  if (itmp == 1) 
  {
    if (inherits(start, "mlteqc"))
      start <- start$B[-base]
    if (is.null(start)) 
    {
      irf <- multiec(mods, se = FALSE, method = "irf", base = base)
      start <- irf$B[-base]
    }
  }

  ini <- start

  se.A <- se.B <- rep(NA, num.forms)
  varAB <- NULL

  X <- model.matrix(~factor(items) - 1, data = itmpar)
  row.names(X) <- itmpar$items.t
  colnames(X) <- substr(colnames(X), 14, 100)
  X_list <- list()
  for (f in 1:num.forms) X_list[[f]] <- X[itmpar$t == f, ]
  pos <- match(itmpar$items, colnames(X))

  if (itmp >= 2)
    opt <- try(optim(par = ini, fn = profLikRcpp, coef = itmpar$coef, t = itmpar$t - 1, X_list = X_list, itmvar = itmvar, numforms = num.forms, notbase = (1:num.forms)[-base] - 1, DffcltNum = DffcltNum - 1, DscrmnNum = DscrmnNum - 1, pos = pos - 1, method = "BFGS", hessian = se, control = list(maxit = iter.max)))
  if (itmp == 1)
    opt <- try(optim(par = ini, fn = profLikRcpp_1PL, coef = itmpar$coef, t = itmpar$t - 1, X_list = X_list, itmvar = itmvar, numforms = num.forms, notbase = (1:num.forms)[-base] - 1, pos = pos - 1, method = "BFGS", hessian = se, control = list(maxit = iter.max)))
  if (trace) cat(" .   \n")
  if (!isa(opt, "try-error")) 
  {
    A <- rep(1, num.forms)
    B <- rep(0, num.forms)
    if (itmp >= 2) 
    {
      A[-base] <- opt$par[1:(num.forms - 1)]
      B[-base] <- opt$par[num.forms:(2 * (num.forms - 1))]
    }
    if (itmp == 1) 
      B[-base] <- opt$par
    modsnames <- names(mods)
    names(A) <- names(B) <- modsnames
    conv <- opt$convergence

    if (itmp >= 2) 
    {
      outBetas <- getBetas(par = opt$par, itmpar = itmpar, itmvar = itmvar, num.forms = num.forms, base = base, DffcltNum, DscrmnNum, X_list = X_list)
      betas <- outBetas$beta
      as <- betas[grep("Dscrmn", rownames(betas)), ]
      bs <- betas[grep("Dffclt", rownames(betas)), ]
    }
    if (itmp == 1) 
    {
      outBetas <- getBetas_1PL(par = opt$par, itmpar = itmpar, itmvar = itmvar, num.forms = num.forms, base = base, X_list = X_list)
      betas <- outBetas$beta
      rownames(betas) <- substr(colnames(X), 14, 100)
      bs <- betas[grep("Dffclt", rownames(betas)), ]
      as <- NULL
    }

    itmpar$A <- A[itmpar$t]
    itmpar$B <- B[itmpar$t]
    itmpar$Y <- 9999
    itmpar[itemstype == "Dscrmn", Y := coef / A]
    itmpar[itemstype == "Dffclt", Y := coef * A + B]

    tab <- data.table::dcast(itmpar, items ~ t, value.var = c("coef", "Y"))
    tab <- as.data.frame(tab)
    sel <- which(colnames(tab) == paste("Y", base, sep = "_"))
    tab <- tab[, -sel]  # delete base form converted
    colnames(tab) <- c("Item", modsnames, paste(modsnames[-base], modsnames[base], sep = ".as."))

    partial <- NULL

    if (se) 
    {
      if (trace) cat("Computation of standard errors ")
      if (obsinf) 
      {
        varAB <- Matrix::solve(opt$hessian)
        if (trace) cat(" .  .  .  . \n")
      } 
      else 
      {
        if (trace) cat(" . ")
        derS_AB <- opt$hessian
        common <- as.character(itmpar$items[duplicated(itmpar$items)])
        itmpar_common <- itmpar[items %in% common]
        gamma <- itmpar_common$coef
        names(gamma) <- itmpar_common$items.t
        if (itmp >= 2)
          derS_gamma <- jacobian(func = derAB, x = gamma, par = opt$par, itmpar = itmpar, itmvar = itmvar, num.forms = num.forms, method = "simple", base = base, DffcltNum = DffcltNum, DscrmnNum = DscrmnNum, X_list = X_list, pos = pos)
        if (itmp == 1)
          derS_gamma <- jacobian(func = derAB_1PL, x = gamma, par = opt$par, itmpar = itmpar, itmvar = itmvar, num.forms = num.forms, method = "simple", base = base, X_list = X_list, pos = pos)
        if (trace) cat(" . ")
        colnames(derS_gamma) <- itmpar_common$items.t
        # derAB_gamma <- -Matrix::solve(derS_AB) %*% derS_gamma
        derAB_gamma <- -Matrix::solve(derS_AB, derS_gamma)

        partial <- t(derAB_gamma)  # derivatives of A and B equating coefficients with respect to the item parameters
        if (itmp == 1)
          partial <- cbind(matrix(0, nrow(partial), num.forms - 1), partial)
        namesAB <- c(paste("A", (1:num.forms)[-base], sep = "."), paste("B", (1:num.forms)[-base], sep = "."))
        colnames(partial) <- namesAB

        var_gamma <- bdiag(itmvar)
        
        VarNames_list <- sapply(itmvar, function(x) rownames(x))
        VarNames <- as.vector(VarNames_list)
        rownames(var_gamma) <- colnames(var_gamma) <- VarNames

        sel <- colnames(derAB_gamma)
        if (trace) cat(" . ")
        vargamma_tderAB_gamma <- Matrix::tcrossprod(var_gamma[sel, sel], derAB_gamma)
        varAB <- derAB_gamma %*% vargamma_tderAB_gamma
        if (trace) cat(" . \n")
      }
      se.A <- se.B <- rep(0, num.forms)
      names(se.A) <- names(se.B) <- modsnames
      if (itmp >= 2) 
      {
        se.A[-base] <- Matrix::diag(varAB)[1:(num.forms - 1)]^0.5
        se.B[-base] <- Matrix::diag(varAB)[(num.forms):(2 * num.forms - 2)]^0.5
        namesAB <- c(paste("A", (1:num.forms)[-base], sep = "."), paste("B", (1:num.forms)[-base], sep = "."))
        colnames(varAB) <- rownames(varAB) <- namesAB
      }
      if (itmp == 1) 
      {
        se.B[-base] <- Matrix::diag(varAB)^0.5
        tmp <- matrix(0, (2 * num.forms - 2), (2 * num.forms - 2))
        tmp[num.forms:(2 * num.forms - 2), num.forms:(2 * num.forms - 2)] <- as.matrix(varAB)
        varAB <- tmp
        namesAB <- c(paste("A", (1:num.forms)[-base], sep = "."), paste("B", (1:num.forms)[-base], sep = "."))
        colnames(varAB) <- rownames(varAB) <- namesAB
      }

      # computation of se of as and bs
      tX.Vinv.X <- outBetas$tX.Vinv.X
      V.inv.X <- outBetas$V.inv.X
      der_beta_gamma <- Matrix::tcrossprod(solve(tX.Vinv.X), V.inv.X)

      if (obsinf) 
      {
        var_gamma <- bdiag(itmvar)
        VarNames_list <- sapply(itmvar, function(x) rownames(x))
        VarNames <- as.vector(VarNames_list)
        rownames(var_gamma) <- colnames(var_gamma) <- VarNames
      }
      sel <- colnames(der_beta_gamma)
      var_gamma_tder_beta_gamma <- Matrix::tcrossprod(var_gamma[sel, sel], der_beta_gamma)
      var_ab <- der_beta_gamma %*% var_gamma_tder_beta_gamma
      se_ab <- Matrix::diag(var_ab)^0.5
      se.as <- se_ab[grep("Dscrmn", names(se_ab))]
      se.bs <- se_ab[grep("Dffclt", names(se_ab))]
    }
  } 
  else
  {
    A <- B <- rep(NA, length(se.A))
    conv <- "optimization failed"
    as <- bs <- NULL
  }
  out <- list(A = A, B = B, se.A = se.A, se.B = se.B, varAB = varAB, as = as, bs = bs, se.as = se.as, se.bs = se.bs, tab = tab, varFull = itmvar, partial = partial, itmp = itmp, method = "lik", basename = modsnames[base], convergence = conv)
  class(out) <- "mlteqc"
  return(out)
}

# this function is used for finding the second derivatives of the profile log-likelihood
# with respect to the equating coefficients and the item parameters
# x = item parameters
derAB <- function(x, par, itmpar, itmvar, num.forms, base, DffcltNum, DscrmnNum, X_list, pos) 
{
  itmpar$coef <- x[itmpar$items.t]
  grad(func = profLikRcpp, x = par, coef = itmpar$coef, t = itmpar$t - 1, X_list = X_list, itmvar = itmvar, numforms = num.forms, notbase = (1:num.forms)[-base] - 1, DffcltNum = DffcltNum - 1, DscrmnNum = DscrmnNum - 1, pos = pos - 1)
}

derAB_1PL <- function(x, par, itmpar, itmvar, num.forms, base, X_list, pos) 
{
  itmpar$coef <- x[itmpar$items.t]
  grad(func = profLikRcpp_1PL, x = par, coef = itmpar$coef, t = itmpar$t - 1, X_list = X_list, itmvar = itmvar, numforms = num.forms, notbase = (1:num.forms)[-base] - 1, pos = pos - 1)
}



plot.mlteqc <- function(x, form, ask = prod(par("mfcol")) < x$itmp*2 && dev.interactive(), ...)
{
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  basef <- x$basename
  if (is.numeric(form)) form <- names(x$A)[form]
  form2 <- paste(form, basef, sep = ".as.")
  tab <- itm(x)[, c("Item", basef, form, form2)]
  type <- as.factor(substr(tab$Item, 1, 6))
  tabDffclt <- tab[type == "Dffclt",]
  if (!any(complete.cases(tabDffclt))) message("Form ", x$basename, " and Form ", form, " don't have common items.")
  rnDff <- range(tabDffclt[, -1], na.rm = TRUE)
  plot(tabDffclt[, c(form, basef)], xlim = rnDff, ylim = rnDff, main = "Difficulty parameters")
  abline(0, 1)
  abline(x$B[form], x$A[form], col = 5, lwd = 2)
  plot(tabDffclt[, c(form2, basef)], xlim = rnDff, ylim = rnDff, main = "Difficulty parameters converted")
  abline(0, 1)
  if (x$itmp>1) 
  {
    tabDscrmn <- tab[type == "Dscrmn",]
    rnDsc <- range(tabDscrmn[, -1], na.rm = TRUE)
    plot(tabDscrmn[, c(form, basef)], xlim = rnDsc, ylim = rnDsc, main = "Discrimination parameters")
    abline(0, 1)
    abline(0, 1/x$A[form], col = 5, lwd = 2)
    plot(tabDscrmn[, c(form2, basef)], xlim = rnDsc, ylim = rnDsc, main = "Discrimination parameters converted")
    abline(0, 1)
  }
}



