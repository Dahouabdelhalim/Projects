testsetup <- function (liminput, ...) 
{
    calcterm <- function(term, usecomp = TRUE) {
        calc <- term$val
        if (!is.na(term$par1)) 
            calc <- calc * parval[term$par1]
        if (!is.na(term$par2)) 
            calc <- calc * parval[term$par2]
        if (!is.na(term$par3)) 
            calc <- calc * parval[term$par3]
        if (!is.na(term$par4)) 
            calc <- calc * parval[term$par4]
        if (!is.na(term$comp))
        	  if (usecomp) calc <- calc * compval[term$comp]
        if (!is.na(term$external)) 
            if (usecomp) calc <- calc * externval[term$extern]
        return(calc)
    }
    calcvalues <- function(pars, npar) {
        if (npar == 0) 
            return(NULL)
        parval <- rep(0, length = npar)
        for (i in 1:npar) {
            term <- pars[pars$nr == i, ]
            for (ii in 1:nrow(term)) {
                parval[i] <- parval[i] + calcterm(term[ii, ])
            }
        }
        return(parval)
    }
    calcvariable <- function(vars, ii, vec = rep(0, nunkn + 1)) {
        term <- vars[vars$nr == ii, ]
        for (i in 1:nrow(term)) {
            ct <- calcterm(term[i, ])
            iv <- term$var[i]
            ifl <- term[i, unkncol]
            if (!is.na(iv)) 
                vec[1:nunkn] <- vec[1:nunkn] + varmat[iv, 1:nunkn] * 
                  ct
            if (!is.na(ifl)) {
                vec[ifl] <- vec[ifl] + ct
            }
            if (is.na(iv) & is.na(ifl)) 
                vec[nunkn + 1] <- vec[nunkn + 1] - ct
        }
        return(vec)
    }
    ncomp <- length(liminput$compnames)
    nextern <- length(liminput$externnames)
    ifelse(is.null(liminput$rate), nrate <- 0, nrate <- ncomp)
    ifelse(is.null(liminput$marker), nmarker <- 0, nmarker <- nrow(liminput$marker))
    npar <- length(liminput$parnames)
    nvars <- length(liminput$varnames)
    flows <- liminput$flows
    reactions <- liminput$reactions
    ispos <- TRUE
    ifelse(is.null(reactions), nreac <- 0, nreac <- max(reactions$nr))
    nflow <- nrow(flows)
    neqs <- nineqs <- 0
    if (is.null(nflow)) 
        nflow <- 0
    if (liminput$Type == "flow") 
        Unknown <- "F"
    else if (liminput$Type == "reaction") 
        Unknown <- "R"
    else Unknown <- "S"
    if (Unknown == "F") {
        nunkn <- nflow
        if (nunkn == 0) 
            stop("cannot create inverse matrices: number of flows=0")
        unkncol <- 9
        Unknownnames <- liminput$flows$name
        if (any(!liminput$posreac)) 
            ispos <- FALSE
    }
    else if (Unknown == "R") {
        nunkn <- max(reactions$nr)
        if (nunkn == 0) 
            stop("cannot create inverse matrices: number of unknowns =0")
        unkncol <- 12
        Unknownnames <- unique(liminput$reactions$name)
        if (any(!liminput$posreac)) 
            ispos <- FALSE
    }
    else {
        nunkn <- ncomp
        unkncol <- 10
        if (nunkn == 0) 
            stop("cannot create inverse matrices: number of components=0")
        Unknownnames <- liminput$compnames
    }
    parval <- rep(NA, npar)
    compval <- rep(NA, ncomp)
    if (!is.null(liminput$pars)) {
        while (any(is.na(parval))) parval <- calcvalues(liminput$pars, 
            npar)
    }
    ifelse(Unknown %in% c("F", "R"), compval <- calcvalues(liminput$comp, 
        ncomp), compval <- rep(1, ncomp))
    externval <- calcvalues(liminput$extern, nextern)
    rateval <- calcvalues(liminput$rate, nrate)
    markerval <- calcvalues(liminput$marker, nmarker)
    if (length(rateval) == 0) 
        rateval <- rep(0, ncomp)
    rateval[is.na(rateval)] <- 0
    VarA <- VarB <- NULL
    if (!is.null(liminput$vars)) {
        varmat <- matrix(data = 0, nrow = nvars, ncol = nunkn + 
            1)
        for (i in 1:nvars) varmat[i, ] <- calcvariable(liminput$vars, 
            i)
        VarA <- varmat[, 1:nunkn]
        VarB <- varmat[, nunkn + 1]
        if (nvars > 1) 
            ii <- which(rowSums(abs(VarA)) == 0)
        else ii <- NULL
        if (length(ii) > 0) {
            print(liminput$varnames[ii])
            stop(paste("cannot proceed: the above variable(s) are  not a true variable, but a parameter:"))
        }
    }
    eqmat <- NULL
    if (!is.null(liminput$equations)) {
        neqs <- max(liminput$equations$nr)
        eqmat <- matrix(data = 0, nrow = neqs, ncol = nunkn + 
            1)
        for (i in 1:neqs) eqmat[i, ] <- calcvariable(liminput$equations, 
            i)
    }
    ineqmat <- NULL
    if (!is.null(liminput$constraints)) {
        nineqs <- max(liminput$constraints$nr)
        ineqmat <- matrix(data = 0, nrow = nineqs, ncol = nunkn + 
            1)
        for (i in 1:nineqs) ineqmat[i, ] <- calcvariable(liminput$constraints, 
            i)
    }
    cost <- NULL
    if (!is.null(liminput$cost)) {
        ncost <- max(liminput$cost$nr)
        costmat <- matrix(data = 0, nrow = ncost, ncol = nunkn + 
            1)
        for (i in 1:ncost) costmat[i, ] <- calcvariable(liminput$cost, 
            i)
        cost <- costmat[, 1:nunkn]
    }
    profit <- NULL
    if (!is.null(liminput$profit)) {
        nprofit <- max(liminput$profit$nr)
        profitmat <- matrix(data = 0, nrow = nprofit, ncol = nunkn + 
            1)
        for (i in 1:nprofit) profitmat[i, ] <- calcvariable(liminput$profit, 
            i)
        profit <- profitmat[, 1:nunkn]
    }
    A <- B <- G <- H <- NULL
    if (nflow > 0) {
        A <- matrix(data = 0, nrow = ncomp, ncol = nunkn)
        for (i in 1:nflow) {
            if (flows$from[i] > 0) 
                A[flows$from[i], i] <- A[flows$from[i], i] - 
                  1
            if (flows$to[i] > 0) 
                A[flows$to[i], i] <- A[flows$to[i], i] + 1
        }
        B <- rateval
    }
    else if (nreac > 0) {
        A <- matrix(data = 0, nrow = ncomp, ncol = nunkn)
        for (i in 1:nrow(reactions)) {
            st <- reactions[i, "comp"]
            rr <- reactions[i, "nr"]
            A[st, rr] <- A[st, rr] + calcterm(reactions[i, ], usecomp = FALSE)
        }
        B <- rateval
    }
    if (!is.null(eqmat)) {
        A <- rbind(A, eqmat[, 1:nunkn])
        B <- c(B, eqmat[, nunkn + 1])
    }
    if (nflow > 0 || nreac > 0) {
        Aext <- cbind(A[1:ncomp, ], rateval)
        Aext[Aext != 0] <- 1
        RS <- rowSums(Aext)
        if (any(RS <= 1)) {
            warning("following component(s) is possibly an external:")
            warning(paste(" ", liminput$compnames[which(RS <= 
                1)]))
        }
        G <- diag(nunkn)
        G <- G[liminput$posreac, ]
        H <- rep(0, sum(liminput$posreac))
    }
    if (!is.null(ineqmat)) {
        G <- rbind(ineqmat[, 1:nunkn], G)
        H <- c(ineqmat[, nunkn + 1], H)
    }
    if (nflow > 0) {
        nel <- ncomp + nextern
        elnames <- c(liminput$compnames, liminput$externnames)
        Flowmatrix <- matrix(data = 0, nrow = nel, ncol = nel, 
            dimnames = list(elnames, elnames))
        flowsft <- as.matrix(flows[, 1:2])
        flowsft[flowsft < 0] <- -flowsft[flowsft < 0] + ncomp
        Flowmatrix[flowsft] <- 1:nflow
    }
    else Flowmatrix <- NULL
    if (nmarker == 0) {
        markers = NULL
    }
    else markers = data.frame(name = c(liminput$compnames, liminput$externnames), 
        val = markerval)
    res <- list(file = liminput$file, NUnknowns = nunkn, NEquations = neqs, 
        NConstraints = nineqs, NComponents = ncomp, NExternal = nrow(liminput$extern), 
        NVariables = nvars, A = A, B = B, G = G, H = H, Cost = cost, 
        Profit = profit, Flowmatrix = Flowmatrix, VarA = VarA, 
        VarB = VarB, Parameters = data.frame(name = liminput$parnames, 
            val = parval), Components = data.frame(name = liminput$compnames, 
            val = compval), Externals = data.frame(name = liminput$externnames, 
            val = externval), rates = data.frame(name = liminput$compnames, 
            val = rateval), markers = markers, Variables = liminput$varnames, 
        costnames = unique(liminput$cost$name), profitnames = unique(liminput$profit$name), 
        eqnames = unique(liminput$equations$name), ineqnames = unique(liminput$constraints$name), 
        Unknowns = Unknownnames, ispos = ispos)
    class(res) <- "lim"
    return(res)
}
