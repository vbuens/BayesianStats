###########################################################################
# VISUALISING THE BETA AND GAMMA DISTRIBUTIONS
#
# These functions plot the beta and gamma distribution and allow the user
# to vary the parameters with Tcl/Tk sliders.
#
# To run with default values, just type "vis.beta()" or "vis.gamma()".
# Inital values of parameters can be specified if desired.
#
# Note that R will freeze until the slider window is closed. The slider window
# may disappear behind other windows: look for the "TK" icon on the task bar.
#
# The functions return the last parameter values selected by the user.
#
# Code by Mike Meredith, updated 27 May 2007
# Inspired by functions by Greg Snow in the 'TeachingDemos' package.
# More updates by Bob Burn, 10 July 2007 and by Sandro Leidi, August 2014
###########################################################################

vis.beta <- function (Alpha=2, Beta=2) 
{
 ## Preliminaries:
    if (!exists("slider.env")) 
        slider.env <<- new.env()
    library(tcltk)
    assign("Alpha", tclVar(Alpha), env = slider.env)
    assign("Beta", tclVar(Beta), env = slider.env)
 ## This draws the figure:
    Beta.refresh <- function(...) {
        Alpha <- as.numeric(evalq(tclvalue(Alpha), env = slider.env))
        Beta <- as.numeric(evalq(tclvalue(Beta), env = slider.env))
        mu <- Alpha / (Alpha+Beta)
        md <- ifelse(Alpha+Beta != 2, max(0,min(1,(Alpha-1)/(Alpha+Beta-2))), 0.5)
        SD<- sqrt((Alpha * Beta) / ((Alpha+Beta)^2 * (Alpha + Beta + 1)))
        ylim <- c(0,4)
        plot(0:1, ylim, type='n',main="", xlab = expression(theta), ylab = "Probability density")
        curve(dbeta(x, Alpha, Beta), 0, 1, col='red', lwd=2, add=TRUE)
        segments(mu,0,mu,dbeta(mu, Alpha, Beta), lty=2, col='blue')
        segments(md,0,md,dbeta(md, Alpha, Beta), lty=2, col='green')
        title(main='Beta Distribution')
        x.pos <- ifelse(md>0.5,0,0.75)
        x.pos <- ifelse((Alpha<1)|(Beta<1),0.4,x.pos)
        text(0.4,4,bquote(alpha == .(Alpha)),pos=4)
        text(0.4,3.8,bquote(beta == .(Beta)),pos=4)
        text(x.pos,2,paste('Mean = ', round(mu,3)),pos=4,col='blue')
        text(x.pos,1.75,paste('Mode = ', round(md,3)),pos=4,col='green')
        text(x.pos,1.5,paste('S.D.   = ', round(SD,3)),pos=4)
    }
 ## Initial display of figure
    Beta.refresh()
    bringToTop()

 ## Set up and run the widget:
    m <- tktoplevel()
    tkwm.title(m, "Beta Distribution")
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Alpha", width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = Beta.refresh, from = 0.1, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = T), 
        side = "left")
    assign("sc", sc, env = slider.env)
    evalq(tkconfigure(sc, variable = Alpha), env = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Beta", width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = Beta.refresh, from = 0.1, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = T), 
        side = "left")
    assign("sc", sc, env = slider.env)
    evalq(tkconfigure(sc, variable = Beta), env = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tkbutton(m, text = "Refresh", command = Beta.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
    cat("Waiting for TK slider window to close...") ; flush.console()
    tkwait.window(m)
    cat("okay.\n")

 ## When window is closed, return final values:
    Alpha <- as.numeric(evalq(tclvalue(Alpha), env = slider.env))
    Beta <- as.numeric(evalq(tclvalue(Beta), env = slider.env))
    output <- c(Alpha, Beta)
    names(output) <- c("Alpha", "Beta")
    rm("slider.env", pos=1)
    return(output)
}

###########################################################################

vis.gamma <- function (Alpha=2, Beta=2) 
{
 ## Preliminaries:
    if (!exists("slider.env")) 
        slider.env <<- new.env()
    library(tcltk)
    assign("Alpha", tclVar(Alpha), env = slider.env)
    assign("Beta", tclVar(Beta), env = slider.env)
 ## This draws the figure:
    Gamma.refresh <- function(...) {
        Alpha <- as.numeric(evalq(tclvalue(Alpha), env = slider.env))
        Beta <- as.numeric(evalq(tclvalue(Beta), env = slider.env))
        mu <- Alpha / Beta
        md <- max(0,(Alpha-1) / Beta)
        SD<- sqrt(Alpha) / Beta
        p25<-qgamma(0.025,Alpha,Beta)
        p975<-qgamma(0.975,Alpha,Beta)
        x.max <- 5
        x.pos <- x.max-1.5
        y.seq <- seq(0,2,length=length(0:x.max))
        y.pos <- max(y.seq)-0.2
        plot(0:x.max, y.seq, type='n', xlab = expression(theta), ylab = "Probability density")
        curve(dgamma(x, shape=Alpha, rate=Beta), 0, x.max, col='red', lwd=2, add=TRUE)
        segments(mu,0,mu,dgamma(mu, shape=Alpha, rate=Beta), lty=2, col='blue')
        segments(md,0,md,dgamma(md, shape=Alpha, rate=Beta), lty=2, col='green')
        segments(p25,0,p25,dgamma(p25, shape=Alpha, rate=Beta), lty=1, col='brown')
        segments(p975,0,p975,dgamma(p975, shape=Alpha, rate=Beta), lty=1, col='brown')
        title(main='Gamma Distribution')
        text(x.pos,y.pos,bquote(alpha == .(Alpha)),pos=4)
        text(x.pos,y.pos-0.1,bquote(beta == .(Beta)),pos=4)
        text(x.pos,y.pos-0.3,paste('Mean = ', round(mu,3)),pos=4,col='blue')
        text(x.pos,y.pos-0.4,paste('Mode = ', round(md,3)),pos=4,col='green')
        text(x.pos,y.pos-0.5,paste('S.D.   = ', round(SD,3)),pos=4)
        text(x.pos,y.pos-0.6,paste('P2.5   = ', round(p25,3)),pos=4,col='brown')
        text(x.pos,y.pos-0.7,paste('P97.5   = ', round(p975,3)),pos=4,col='brown')
    }
 ## Initial display of figure
    Gamma.refresh()
    bringToTop()

 ## Set up and run the widget:
    m <- tktoplevel()
    tkwm.title(m, "Gamma Distribution")
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Alpha", width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = Gamma.refresh, from = 0.1, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = T), 
        side = "left")
    assign("sc", sc, env = slider.env)
    evalq(tkconfigure(sc, variable = Alpha), env = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Beta", width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = Gamma.refresh, from = 0.1, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = T), 
        side = "left")
    assign("sc", sc, env = slider.env)
    evalq(tkconfigure(sc, variable = Beta), env = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tkbutton(m, text = "Refresh", command = Gamma.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
    cat("Waiting for TK slider window to close...") ; flush.console()
    tkwait.window(m)
    cat("okay.\n")

 ## When window is closed, return final values:
    Alpha <- as.numeric(evalq(tclvalue(Alpha), env = slider.env))
    Beta <- as.numeric(evalq(tclvalue(Beta), env = slider.env))
    output <- c(Alpha, Beta)
    names(output) <- c("Alpha", "Beta")
    rm("slider.env", pos=1)
    return(output)
}

###########################################################################
