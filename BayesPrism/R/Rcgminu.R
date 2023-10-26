Rcgminu <- function(par, fn, gr, control = list(), ...) {
    ## An R version of the conjugate gradient minimization
    ## using the Dai-Yuan ideas
    #  This version is for unconstrained functions.
    #
    # Input:
    #  par  = a vector containing the starting point
    # fn = objective function (assumed to be sufficeintly
    #   differentiable)
    #  gr = gradient of objective function
    #  control = list of control parameters
    #           maxit = a limit on the number of iterations (default 500)
    #           maximize = TRUE to maximize the function (default FALSE)
    #           trace = 0 (default) for no output,
    #                  >0 for output (bigger => more output)
    # eps=1.0e-7 (default) for use in computing numerical
    #   gradient approximations.
    # dowarn=TRUE by default. Set FALSE to suppress warnings.
    #
    # Output:
    #    A list with components:
    #
    #     par: The best set of parameters found.
    #
    #   value: The value of 'fn' corresponding to 'par'.
    #
    # counts: A two-element integer vector giving the number of
    #   calls to
    # 'fn' and 'gr' respectively. This excludes those calls
    #   needed
    # to compute the Hessian, if requested, and any calls to
    #   'fn'
    # to compute a finite-difference approximation to the
    #   gradient.
    #
    # convergence: An integer code. '0' indicates successful
    #   convergence.
    #          Error codes are
    #          '0' converged
    # '1' indicates that the function evaluation count
    #   'maxfeval'
    #               was reached.
    #          '2' indicates initial point is infeasible
    #
    # message: A character string giving any additional
    #   information returned
    #          by the optimizer, or 'NULL'.
    #
    #
    #  Author:  John C Nash
    #  Date:  April 2, 2009; revised July 28, 2009
    #################################################################
    # control defaults -- idea from spg
    
    #edited by Tinyi add maxRst
    ctrl <- list(maxit = 500, maximize = FALSE, trace = 0, eps = 1e-07, 
        dowarn = TRUE, tol=0, maxNA=500)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    npar<-length(par)
    if (ctrl$tol == 0) tol <- npar * (npar * .Machine$double.eps)  # for gradient test.  
    # Note -- integer overflow if npar*npar*.Machine$double.eps
    else tol<-ctrl$tol
    maxit <- ctrl$maxit  # limit on function evaluations
    maximize <- ctrl$maximize  # TRUE to maximize the function
    trace <- ctrl$trace  # 0 for no output, >0 for output (bigger => more output)
    #edited by Tinyi
    maxNA <- ctrl$maxNA

    
    
    if (trace > 2) 
        cat("trace = ", trace, "\n")
    eps <- ctrl$eps
    fargs <- list(...)  # the ... arguments that are extra function / gradient data
    grNULL <- is.null(gr)
    dowarn <- ctrl$dowarn  #
    #############################################
    if (maximize) {
       warning("Rcgmin no longer supports maximize 111121 -- see documentation")
       msg<-"Rcgmin no longer supports maximize 111121"
       ans <- list(par, NA, c(0, 0), 9999, msg)
       return(ans)
    }
    #############################################
    # gr MUST be provided
    if (grNULL) {
    ##   require(numDeriv) # in NAMESPACE
       if (control$dowarn) 
          warning("A NULL gradient function is being replaced by numDeriv 'grad()'for Rcgmin")
       if (ctrl$trace > 1) {
           cat("Using following function in numDeriv grad()\n")
           print(fn)
       }
       mygr<-function(prm, func=fn, ...){
           gv<-grad(func=func, x=prm, ...)
       }
  #############################################
    } else { mygr <- gr }
  ############# end test gr ####################
  ## Set working parameters (See CNM Alg 22)
    if (trace > 0) {
        cat("Rcgminu -- J C Nash 2009 - unconstrained version CG min\n")
        cat("an R implementation of Alg 22 with Yuan/Dai modification\n")
    }
    bvec <- par  # copy the parameter vector
    n <- length(bvec)  # number of elements in par vector
    maxfeval <- round(sqrt(n + 1) * maxit)  # change 091219
    ig <- 0  # count gradient evaluations
    ifn <- 1  # count function evaluations (we always make 1 try below)
    stepredn <- 0.15  # Step reduction in line search
    acctol <- 1e-04  # acceptable point tolerance
    reltest <- 100  # relative equality test
    accpoint <- as.logical(FALSE)  # so far do not have an acceptable point
    cyclimit <- min(2.5 * n, 10 + sqrt(n))  #!! upper bound on when we restart CG cycle
    #!! getting rid of limit makes it work on negstart BUT inefficient
    # This does not appear to be in Y H Dai & Y Yuan, Annals of
    #   Operations Research 103, 33–47, 2001 aor01.pdf
    # in Alg 22 pascal, we can set this as user. Do we wish to allow that?
    ##    tol <- n * (n * .Machine$double.eps)  # # for gradient test.  
    ## Note -- integer overflow if n*n*d.eps
    fargs <- list(...)  # function arguments
    if (trace > 2) {
        cat("Extra function arguments:")
        print(fargs)
    }
    # Initial function value -- may NOT be at initial point
    #   specified by user.
    if (trace > 2) {
        cat("Try function at initial point:")
        print(bvec)
    }
    f <- try(fn(bvec, ...), silent = TRUE)  # Compute the function at initial point.
    if (trace > 0) {
        cat("Initial function value=", f, "\n")
    }
    if (class(f) == "try-error") {
        msg <- "Initial point is infeasible."
        if (trace > 0) 
            cat(msg, "\n")
        ans <- list(par, NA, c(ifn, 0), 2, msg)
        names(ans) <- c("par", "value", "counts", "convergence", 
            "message")
        return(ans)
    }
    fmin <- f
    if (trace > 0) 
        cat("Initial fn=", f, "\n")
    if (trace > 2) 
        print(bvec)
    # Start the minimization process
    keepgoing <- TRUE
    msg <- "not finished"  # in case we exit somehow
    oldstep <- 0.8  #!! 2/3 #!!?? Why this choice?
    ####################################################################
    fdiff <- NA  # initially no decrease
    cycle <- 0  # !! cycle loop counter
    naD.num <- 0 #count consecutive NAs
    while (keepgoing && (naD.num < maxNA)) {
        # main loop -- must remember to break out of it!!
        t <- as.vector(rep(0, n))  # zero step vector
        c <- t  # zero 'last' gradient
        
        while (keepgoing && (cycle < cyclimit) && (naD.num < maxNA)) {
            ## cycle loop
            cycle <- cycle + 1
            
            fdiff.old <- fdiff
            
            if (trace > 0) 
                cat(ifn, " ", ig, " ", cycle, " ", fmin, "  last decrease=", 
                  fdiff, "\n")
            if (trace > 2) {
                print(bvec)
                cat("\n")
            }
            if (ifn > maxfeval) {
                msg <- paste("Too many function evaluations (> ", 
                  maxfeval, ") ", sep = "")
                if (trace > 0) 
                  cat(msg, "\n")
                ans <- list(par, fmin, c(ifn, ig), 1, msg)  # 1 indicates not converged in function limit
                names(ans) <- c("par", "value", "counts", "convergence", 
                  "message")
                return(ans)
            }
            par <- bvec  # save best parameters
            ig <- ig + 1
            if (ig > maxit) {
                msg <- paste("Too many gradient evaluations (> ", 
                  maxit, ") ", sep = "")
                if (trace > 0) 
                  cat(msg, "\n")
                ans <- list(par, fmin, c(ifn, ig), 1, msg)  # 1 indicates not converged in function or gradient limit
                names(ans) <- c("par", "value", "counts", "convergence", 
                  "message")
                return(ans)
            }
            g <- mygr(bvec, ...)
            g1 <- sum(g * (g - c))  # gradient * grad-difference
            g2 <- sum(t * (g - c))  # oldsearch * grad-difference
            gradsqr <- sum(g * g)
            if (trace > 1) {
                cat("Gradsqr = ", gradsqr, " g1, g2 ", g1, " ", 
                  g2, " fmin=", fmin, "\n")
            }
            c <- g  # save last gradient
            g3 <- 1  # !! Default to 1 to ensure it is defined -- t==0 on first cycle
            if (gradsqr > tol * (abs(fmin) + reltest)) {
                if (g2 > 0) {
                  betaDY <- gradsqr/g2
                  betaHS <- g1/g2
                  g3 <- max(0, min(betaHS, betaDY))  # g3 is our new 'beta' !! Dai/Yuan 2001, (4.2)
                }
            }
            else {
                msg <- paste("Very small gradient -- gradsqr =", 
                  gradsqr, sep = " ")
                if (trace > 0) 
                  cat(msg, "\n")
                keepgoing <- FALSE  # done loops -- should we break ??
                break  # to leave inner loop
            }
            if (trace > 2) 
                cat("Betak = g3 = ", g3, "\n")
            if (g3 == 0 || cycle >= cyclimit) {
                # we are resetting to gradient in this case
                if (trace > 0) {
                  if (cycle < cyclimit) cat("Yuan/Dai cycle reset\n")  
                  else cat("Cycle limit reached -- reset\n")
                }
                fdiff <- NA
                cycle <- 0
                break  #!!
                #!! oldstep<-1 # !!
                #!! don't reset stepsize ## oldstep<-1 #!! reset
                #!! break # to quit inner loop
            }
            else {
                # drop through if not Yuan/Dai cycle reset
                t <- t * g3 - g  # t starts at zero, later is step vector
                gradproj <- sum(t * g)  # gradient projection
                if (trace > 1) 
                  cat("Gradproj =", gradproj, "\n")
                # ?? Why do we not check gradproj size??
                ########################################################
                ####                  Line search                   ####
                OKpoint <- FALSE
                if (trace > 2) 
                  cat("Start linesearch with oldstep=", oldstep, 
                    "\n")
                steplength <- oldstep * 1.5  #!! try a bit bigger
                f <- fmin
                changed <- TRUE  # Need to set so loop will start
                while ((f >= fmin) && changed) {
                  bvec <- par + steplength * t
                  changed <- (!identical((bvec + reltest), (par + reltest)))
                  if (changed) {
                    # compute newstep, if possible
                    f <- fn(bvec, ...)  # Because we need the value for linesearch, don't use try()
                    # instead preferring to fail out, which will hopefully be
                    #   unlikely.
                    ifn <- ifn + 1
                    if (is.na(f) || (!is.finite(f))) {
                      warning("Rcgmin - undefined function")
                      f <- .Machine$double.xmax
                    }
                    if (f < fmin) {
                      f1 <- f  # Hold onto value
                    }
                    else {
                      savestep<-steplength
                      steplength <- steplength * stepredn
                      if (steplength >=savestep) changed<-FALSE
                      if (trace > 0) 
                        cat("*")
                    }
                  }
                }  # end while
                changed1 <- changed  # Change in parameters occured in step reduction
                if (changed1) 
                  {
                    ## ?? should we check for reduction? or is this done in if
                    #   (newstep >0) ?
                    newstep <- 2 * (f - fmin - gradproj * steplength)  # JN 081219 change
                    if (newstep > 0) {
                      newstep = -(gradproj * steplength * steplength/newstep)
                    }
                    bvec <- par + newstep * t
                    changed <- (!identical((bvec + reltest), 
                      (par + reltest)))
                    if (changed) {
                      f <- fn(bvec, ...)
                      ifn <- ifn + 1
                    }
                    if (trace > 2) 
                      cat("fmin, f1, f: ", fmin, f1, f, "\n")
                    if (isTRUE (f < min(fmin, f1))) {
                      # success
                      OKpoint <- TRUE
                      accpoint <- (f <= fmin + gradproj * newstep * 
                        acctol)
                      fdiff <- (fmin - f)  # check decrease
                      fmin <- f
                      oldstep <- newstep  # !! save it
                    }
                    else {
                      if (f1 < fmin) {
                        bvec <- par + steplength * t  # reset best point
                        accpoint <- (f1 <= fmin + gradproj * 
                          steplength * acctol)
                        OKpoint <- TRUE  # Because f1 < fmin
                        fdiff <- (fmin - f1)  # check decrease
                        fmin <- f1
                        oldstep <- steplength  #!! save it
                      }
                      else {
                        # no reduction
                        fdiff <- NA
                        accpoint <- FALSE
                      }  # f1<?fmin
                    }  # f < min(f1, fmin)
                    if (trace > 1) 
                      cat("accpoint = ", accpoint, " OKpoint = ", 
                        OKpoint, "\n")
                    if (!accpoint) {
                      msg <- "No acceptable point -- exit loop"
                      if (trace > 0) 
                        cat("\n", msg, "\n")
                      keepgoing <- FALSE
                      break  #!!
                    }
                  }  # changed1
                else {
                  # not changed on step redn
                  if (cycle == 1) {
                    msg <- " Converged -- no progress on new CG cycle"
                    if (trace > 0) 
                      cat("\n", msg, "\n")
                    keekpgoing <- FALSE
                    break  #!!
                  }
                }  # end else
            }  # end of test on Yuan/Dai condition
            #### End line search ####
        
        	#added by Tinyi
        	if( (is.na(fdiff) && is.na(fdiff.old)) ) naD.num <- naD.num+1 #consecutive NA decrease value
        	else 	naD.num<-0 # reset na count
        
        }  # end of inner loop (cycle)
        if (oldstep < acctol) {
            oldstep <- acctol
        }
        #   steplength
        if (oldstep > 1) {
            oldstep <- 1
        }
        if (trace > 1) 
            cat("End inner loop, cycle =", cycle, "\n")
    }  # end of outer loop
    
    if(naD.num == maxNA)  msg <- paste("maxNA", naD.num,"reached")
    else msg <- "Rcgmin seems to have converged"
    
    if (trace > 0) 
        cat(msg, "\n")
    #  par: The best set of parameters found.
    #  value: The value of 'fn' corresponding to 'par'.
    #  counts: number of calls to 'fn' and 'gr' (2 elements)
    # convergence: An integer code. '0' indicates successful
    #   convergence.
    #  message: A character string or 'NULL'.
#    if (maximize) 
#        fmin <- -fmin
    ans <- list(par, fmin, c(ifn, ig), 0, msg)
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message")
    return(ans)
}  ## end of Rcgminu
