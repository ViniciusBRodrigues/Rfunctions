rdiagnostic <- function(object,...,nmodels=NULL) {
  
  ## Warning about the nmodels argument
  if(is.null(nmodels)!=TRUE) {
    warning("nmodels' argument is not necessary anymore, ?rdiagnostic for more information")
  }
  
  ## Make a list of models to be tested
  objects <- list(object,...)
  
  ## Define if this is plotted
  par(mfrow=c(length(objects),2),pty="s",cex.axis=1.5,cex.lab=1.5)
  
  for(model in objects) {
    ##Para Normal
    if(model$family$family=="gaussian") {
      main <- "Gaussian Model"
      X <- model.matrix(model)
      n <- nrow(X)
      p <- ncol(X)
      H <- X%*%solve(t(X)%*%X)%*%t(X)
      h <- diag(H)
      si <- lm.influence(model)$sigma
      r <- resid(model)
      tsi <- r/(si*sqrt(1-h))
      ident <- diag(n)
      epsilon <- matrix(0,n,100)
      e <- matrix(0,n,100)
      e1 <- numeric(n)
      e2 <- numeric(n)
      for(i in 1:100) {
        epsilon[,i] <- rnorm(n,0,1)
        e[,i] <- (ident - H)%*%epsilon[,i]
        u <- diag(ident - H)
        e[,i] <- e[,i]/sqrt(u)
        e[,i] <- sort(e[,i])
      }
      for(i in 1:n) {
        eo <- sort(e[i,])
        e1[i] <- eo[5]
        e2[i] <- eo[95]
      }
      med <- apply(e,1,mean)
      faixa <- range(tsi,e1,e2)
      qqnorm(tsi,xlab="Theoretical Quantiles",
             ylab="Std. deviance resid.", main=main, ylim=faixa, pch=1)
      par(new=T)
      qqnorm(e1,axes=F,xlab="",ylab="",ylim=faixa,main="",type="l",lty=1)
      par(new=T)
      qqnorm(e2,axes=F,xlab="",ylab="",ylim=faixa,main="", type="l",lty=1)
      par(new=T)
      qqnorm(med,axes=F,xlab="",ylab="",ylim=faixa,main="",type="l",lty=2)
      mtext("Normal Q-Q plot", 3, 0.25)
    }
    else {
      ##Para Poisson 
      if(model$family$family=="poisson") {
        main <- "Poisson Model"
        X <- model.matrix(model)
        n <- nrow(X)
        p <- ncol(X)
        w <- model$weights
        W <- diag(w)
        H <- solve(t(X)%*%W%*%X)
        H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
        h <- diag(H)
        td <- resid(model,type="deviance")/sqrt(1-h)
        e <- matrix(0,n,100)
        for(i in 1:100) {
          nresp <- rpois(n, fitted(model))
          fit <- glm(nresp ~ X , family=poisson)
          w <- fit$weights
          W <- diag(w)
          H <- solve(t(X)%*%W%*%X)
          H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
          h <- diag(H)
          e[,i] <- sort(resid(fit,type="deviance")/sqrt(1-h))
        }
        e1 <- numeric(n)
        e2 <- numeric(n)
        for(i in 1:n) {
          eo <- sort(e[i,])
          e1[i] <- eo[5]
          e2[i] <- eo[95]
        }
        med <- apply(e,1,mean)
        faixa <- range(td,e1,e2)
        qqnorm(td,xlab="Theoretical Quantiles",
               ylab="Std. deviance resid.", ylim=faixa,main=main,pch=1)
        par(new=T)
        qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,main="",lty=1)
        par(new=T)
        qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,main="",lty=1)
        par(new=T)
        qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,main="",lty=2)
        mtext("Normal Q-Q plot", 3, 0.25)
      }
      else {
        ##Para Binomial
        if(model$family$family=="binomial") {
          main <- "Binomial Model"
          X <- model.matrix(model)
          n <- nrow(X)
          p <- ncol(X)
          w <- model$weights
          W <- diag(w)
          H <- solve(t(X)%*%W%*%X)
          H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
          h <- diag(H)
          td <- resid(model,type="deviance")/sqrt(1-h)
          e <- matrix(0,n,100)
          for(i in 1:100) {
            dif <- runif(n) - fitted(model)
            dif[ dif >= 0 ] <- 0
            dif[ dif < 0] <- 1
            nresp <- dif
            fit <- glm(nresp ~ X, family=binomial)
            w <- fit$weights
            W <- diag(w)
            H <- solve(t(X)%*%W%*%X)
            H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
            h <- diag(H)
            e[,i] <- sort(resid(fit,type="deviance")/sqrt(1-h))
          }
          e1 <- numeric(n)
          e2 <- numeric(n)
          for(i in 1:n) {
            eo <- sort(e[i,])
            e1[i] <- eo[5]
            e2[i] <- eo[95]
          }
          med <- apply(e,1,mean)
          faixa <- range(td,e1,e2)
          qqnorm(td,xlab="Theoretical Quantiles",
                 ylab="Std. deviance resid.", ylim=faixa, main=main,pch=1)
          par(new=T)
          qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,main="",lty=1)
          par(new=T)
          qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,main="",lty=1)
          par(new=T)
          qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,main="",lty=2)
          mtext("Normal Q-Q plot", 3, 0.25)
        }
        else {
          ## Para Gamma
          if(model$family$family=="Gamma") {
            main <- "Gamma Model"
            X <- model.matrix(model)
            n <- nrow(X)
            p <- ncol(X)
            w <- model$weights
            W <- diag(w)
            H <- solve(t(X)%*%W%*%X)
            H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
            h <- diag(H)
            ro <- resid(model,type="response")
            fi <- (n-p)/sum((ro/(fitted(model)))^ 2)
            td <- resid(model,type="deviance")*sqrt(fi/(1-h))
            e <- matrix(0,n,100)
            for(i in 1:100) {
              resp <- rgamma(n,fi)
              resp <- (fitted(model)/fi)*resp
              fit <- glm(resp ~ X, family=Gamma(link=log))
              w <- fit$weights
              W <- diag(w)
              H <- solve(t(X)%*%W%*%X)
              H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
              h <- diag(H)
              ro <- resid(fit,type="response")
              phi <- (n-p)/sum((ro/(fitted(fit)))^ 2)
              e[,i] <- sort(resid(fit,type="deviance")*sqrt(phi/(1-h)))
            }
            e1 <- numeric(n)
            e2 <- numeric(n)
            for(i in 1:n) {
              eo <- sort(e[i,])
              e1[i] <- eo[5]
              e2[i] <- eo[95]
            }
            med <- apply(e,1,mean)
            faixa <- range(td,e1,e2)
            qqnorm(td,xlab="Theoretical Quantiles",
                   ylab="Std. deviance resid.", ylim=faixa,main=main, pch=1)
            par(new=T)
            qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,main="",lty=1)
            par(new=T)
            qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,main="",lty=1)
            par(new=T)
            qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,main="",lty=2)
            mtext("Normal Q-Q plot", 3, 0.25)
          }
          else {
            if(strsplit(model$family[[1]],split="\\(")[[1]][1]=="Negative Binomial") {
              main <- "Binomial negative Model"
              ##Para Binomial Negativa
              library(MASS)
              X <- model.matrix(model)
              n <- nrow(X)
              p <- ncol(X)
              fi <- model$theta
              w <- fi*fitted(model)/(fi + fitted(model))
              W <- diag(w)
              H <- solve(t(X)%*%W%*%X)
              H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
              h <- diag(H)
              td <- resid(model,type="deviance")/sqrt(1-h)
              fi <- model$theta
              e <- matrix(0,n,100)
              for(i in 1:100) {
                resp <- rnegbin(n, fitted(model),fi)
                fit <- glm.nb(resp ~ X,maxit=1000)
                w <- fit$weights
                W <- diag(w)
                H <- solve(t(X)%*%W%*%X)
                H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
                h <- diag(H)
                e[,i] <- sort(resid(fit,type="deviance")/sqrt(1-h))
              }
              e1 <- numeric(n)
              e2 <- numeric(n)
              for(i in 1:n) {
                eo <- sort(e[i,])
                e1[i] <- eo[5]
                e2[i] <- eo[95]
              }
              med <- apply(e,1,mean)
              faixa <- range(td,e1,e2)
              qqnorm(td,xlab="Theoretical Quantiles",
                     ylab="Std. deviance resid.", ylim=faixa,main=main, pch=1)
              par(new=T)
              qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,main="",lty=1)
              par(new=T)
              qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,main="",lty=1)
              par(new=T)
              qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,main="",lty=2)
              mtext("Normal Q-Q plot", 3, 0.25)
            }
            else {
              stop("Unknow Distribution")
            }
          }
        }
      }
    }
    plot(model,which=1:1,main=main)
  }
  par(mfrow=c(1,1),pty="m",cex.axis=1,cex.lab=1)
}