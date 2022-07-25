#new
"localLsector" <-
  function(X, ..., rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL,
           begin=0, end=360, units=c("degrees", "radians"), domain = NULL)
  {
    localKsector(X, wantL=TRUE, rmax = rmax,
           correction=correction, verbose=verbose, rvalue=rvalue,
           begin=begin, end=end, units=units, domain = domain)
  }

"localLinhomsector" <-
  function(X, lambda=NULL, ..., rmax = NULL, correction="Ripley",
           verbose=TRUE, rvalue=NULL,
           sigma=NULL, varcov=NULL, update=TRUE, leaveoneout=TRUE,
           begin=0, end=360, units=c("degrees", "radians"), domain = NULL)
  {
    localKinhomsector(X, lambda=lambda, wantL=TRUE, ..., rmax = rmax, 
                correction=correction, verbose=verbose, rvalue=rvalue,
                sigma=sigma, varcov=varcov,
                update=update, leaveoneout=leaveoneout,
                begin=begin, end=end, units=units, domain = domain)
  }

"localKsector" <-
  function(X, ..., rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL,
           begin=0, end=360, units=c("degrees", "radians"), domain = NULL)
  {
    verifyclass(X, "ppp")
    localKsectorengine(X, ..., rmax = rmax,
                 correction=correction, verbose=verbose, rvalue=rvalue,
                 begin=begin, end=end, units=units, domain = domain)
  }

"localKinhomsector" <-
  function(X, lambda=NULL, ...,
           rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL,
           sigma=NULL, varcov=NULL, update=TRUE, leaveoneout=TRUE, 
           begin=0, end=360, units=c("degrees", "radians"), domain = NULL)
  {
    verifyclass(X, "ppp")
    
    a <- resolve.lambda(X, lambda, ...,
                        sigma=sigma, varcov=varcov,
                        update=update, leaveoneout=leaveoneout)
    result <- localKsectorengine(X, lambda=a$lambda, ..., rmax = rmax,
                           correction=correction,
                           verbose=verbose, rvalue=rvalue,
                           begin=begin, end=end, units=units, domain = domain)
    if(a$danger)
      attr(result, "dangerous") <- a$dangerous
    return(result)
  }
#new

"localKsectorengine" <-
  function(X, ..., wantL=FALSE, lambda=NULL, rmax = NULL,
           correction="Ripley", verbose=TRUE, rvalue=NULL, 
           begin=0, end=360, units=c("degrees", "radians"), domain = NULL)
  {
    npts <- npoints(X)
    W <- X$window
    areaW <- area.owin(W)
    lambda.ave <- npts/areaW
    lambda1.ave <- (npts - 1)/areaW
    
    weighted <- !is.null(lambda)
    
    if(is.null(rvalue)) 
      rmaxdefault <- rmax %orifnull% rmax.rule("K", W, lambda.ave)
    else {
      stopifnot(is.numeric(rvalue))
      stopifnot(length(rvalue) == 1)
      stopifnot(rvalue >= 0)
      rmaxdefault <- rvalue
    }
    breaks <- handle.r.b.args(NULL, NULL, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
    
    #new
    if(!is.null(domain)) {
      domain <- as.owin(domain)
      stopifnot(is.subset.owin(domain, Window(X)))
      areaW <- area(domain)
    }
    
    units <- match.arg(units)
    switch(units,
           radians = {
             if(missing(end)) end <- 2 * pi
             check.1.real(begin)
             check.1.real(end)
             check.in.range(begin, c(-pi, 2*pi))
             check.in.range(end, c(0, 2*pi))
             stopifnot(begin < end)
             stopifnot((end - begin) <= 2 * pi)
             BEGIN <- begin
             END   <- end
             Bname <- simplenumber(begin/pi, "pi") %orifnull% signif(begin, 3)
             Ename <- simplenumber(end/pi, "pi") %orifnull% signif(end, 3)
           },
           degrees = {
             check.1.real(begin)
             check.1.real(end)
             check.in.range(begin, c(-90, 360))
             check.in.range(end, c(0, 360))
             stopifnot(begin < end)
             stopifnot((end - begin) <= 360)
             if(verbose && (end - begin) <= 2 * pi)
               warning("Very small interval in degrees: did you mean radians?")
             BEGIN <- pi* (begin/180)
             END   <- pi * (end/180)
             Bname <- signif(begin, 3)
             Ename <- signif(end, 3)
           })
    #new
    
    correction.given <- !missing(correction)
    correction <- pickoption("correction", correction,
                             c(none="none",
                               isotropic="isotropic",
                               Ripley="isotropic",
                               trans="translate",
                               translate="translate",
                               translation="translate",
                               best="best"),
                             multi=FALSE)
    
    correction <- implemented.for.K(correction, W$type, correction.given)
    
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
    
    # identify all close pairs
    #print(r)
    rmax <- max(r)
    #print(rmax)
    close <- as.data.frame(closepairs(X, rmax)) # Maximum distance between pairs of points to be counted as close pairs.
    
    #print(nrow(close))
    
    #new
    if(!is.null(domain)) {
      ## restrict to pairs with first point in 'domain'
      indom <- with(close, inside.owin(xi, yi, domain))
      close <- close[indom, , drop=FALSE]
    }
    
    ## select pairs in angular range
    ang <- with(close, atan2(dy, dx)) %% (2*pi)
    if(BEGIN >= 0) {
      ## 0 <= begin < end
      ok <- (BEGIN <= ang) & (ang <= END)
    } else {
      ## begin < 0 <= end
      ok <- (ang >= 2 * pi + BEGIN) | (ang <= END)
    }
    close <- close[ok, , drop=FALSE]
    #new
    
    DIJ <- close$d
    XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
    I <- close$i
    if(weighted) {
      J <- close$j
      lambdaJ <- lambda[J]
      weightJ <- 1/lambdaJ
    } 
    
    # initialise
    df <- as.data.frame(matrix(NA, length(r), npts))
    labl <- desc <- character(npts)
    
    if(verbose) state <- list()
    
    switch(correction,
           none={
             # uncorrected! For demonstration purposes only!
             for(i in 1:npts) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val,
                           if(weighted) weightJ[ii] else NULL)  # no edge weights
               df[,i] <- cumsum(wh)
               icode <- numalign(i, npts)
               names(df)[i] <- paste("un", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("uncorrected estimate of %s",
                                "for point", icode)
               if(verbose) state <- progressreport(i, npts, state=state)
             }
             if(!weighted) df <- df/lambda1.ave
           },
           translate={
             # Translation correction
             XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
             edgewt <- edge.Trans(XI, XJ, paired=TRUE)
             if(weighted)
               edgewt <- edgewt * weightJ
             for(i in 1:npts) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
               Ktrans <- cumsum(wh)
               df[,i] <- Ktrans
               icode <- numalign(i, npts)
               names(df)[i] <- paste("trans", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("translation-corrected estimate of %s",
                                "for point", icode)
               if(verbose) state <- progressreport(i, npts, state=state)
             }
             if(!weighted) df <- df/lambda1.ave
             h <- diameter.owin(W)/2
             df[r >= h, ] <- NA
           },
           isotropic={
             # Ripley isotropic correction
             edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
             if(weighted)
               edgewt <- edgewt * weightJ
             for(i in 1:npts) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
               Kiso <- cumsum(wh)
               df[,i] <- Kiso
               icode <- numalign(i, npts)
               names(df)[i] <- paste("iso", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("Ripley isotropic correction estimate of %s", 
                                "for point", icode)
               if(verbose) state <- progressreport(i, npts, state=state)
             }
             if(!weighted) df <- df/lambda1.ave
             h <- diameter.owin(W)/2
             df[r >= h, ] <- NA
           })
    # transform values if L required
    if(wantL){
      print("wantL=TRUE")
      df <- sqrt(df/pi)
    }
    
    # return vector of values at r=rvalue, if desired
    if(!is.null(rvalue)) {
      nr <- length(r)
      if(r[nr] != rvalue)
        stop("Internal error - rvalue not attained")
      return(as.numeric(df[nr,]))
    }
    # function value table required
    # add r and theo
    if(!wantL) {
      df <- cbind(df, data.frame(r=r, theo=pi * r^2))
      if(!weighted) {
        fnam <- c("K", "loc")
        yexp <- ylab <- quote(K[loc](r))
      } else {
        fnam <- c("K", "list(inhom,loc)")
        ylab <- quote(K[inhom,loc](r))
        yexp <- quote(K[list(inhom,loc)](r))
      }
    } else {
      df <- cbind(df, data.frame(r=r, theo=r))
      if(!weighted) {
        fnam <- c("L", "loc")
        yexp <- ylab <- quote(L[loc](r))
      } else {
        fnam <- c("L", "list(inhom,loc)")
        ylab <- quote(L[inhom,loc](r))
        yexp <- quote(L[list(inhom,loc)](r))
      }
    }
    desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
    labl <- c(labl, c("r", "{%s[%s]^{pois}}(r)"))
    # create fv object
    K <- fv(df, "r", ylab, "theo", , alim, labl, desc,
            fname=fnam, yexp=yexp)
    # default is to display them all
    formula(K) <- . ~ r
    unitname(K) <- unitname(X)
    attr(K, "correction") <- correction
    return(K)
  }
