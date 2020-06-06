
# R programs to combine asymmetric "systematic" errors
# replacing java code
# R Barlow  

combineasym <- function(sp,sm,model=1){ 
	# sp and sm are vectors of all the sigma+, sigma- values
	#      to be combined
	#    returns a named vector 
 
	pi <- 2*acos(0)  #   pi may be set up, but this makes sure
	rtwopi <- sqrt(2*pi)
	
	#    calculate cumulants
	if(model==1){
           d <- sp-sm     
           mu <- sum(d)/rtwopi
           V <- sum((sp^2+sm^2)/2 - d^2/(2*pi))
	   gamma<- sum(2*(sp^3-sm^3)-1.5*d*(sp^2+sm^2)+d^3/pi)/rtwopi
	 } else {  # for model 2
           alpha=(sp-sm)/2
           mu <- sum(alpha)
           s <- (sp+sm)/2
           V <- sum(s**2+2*alpha**2)
           gamma <- sum(6*s^2*alpha+8*alpha**3)
         }

	# find S and D that give this V and gamma by iteration
	D=0
	for (i in 1:10){  #  10 seems to be plenty: increase if necessary
           oldD <- D
           if(model==1){
            S <- 2*V+D^2/pi
	    D <- 2*(gamma*rtwopi-D^3*(1/pi-1))/(3*S)
	    if(abs(D-oldD) < 1.E-5) 	{ # judged to have converged
	    	ssum <- sqrt(2*S-D^2)
	        return(c(splus=(ssum+D)/2,sminus=(ssum-D)/2,shift=mu-D/rtwopi))
                }  
            } else {  # model 2
              D <- gamma/(6*V-4*D**2)
              if(abs(D-oldD) < 1.E-5) { # converged
                s <- sqrt(V-2*D^2)
                return(c(splus=s+D,sminus=s-D,shift=mu-D))
               }
            }
	    }  # end of iteration-loop
	print("Sorry! Failed to converge in combineasym. Try increasing number ")
}

# test job

print(format(combineasym(c(1,1.2),c(1,.8)),digits=2))
print(format(combineasym(c(1.2,1.2),c(.8,.8)),digits=2))
print(format(combineasym(c(1.5,1.2),c(.5,.8)),digits=2))
print(format(combineasym(c(1.5,1.5),c(.5,.5)),digits=2))
print(format(combineasym(c(1,1.2),c(1,.8),model=2),digits=3))
print(format(combineasym(c(1.2,1.2),c(.8,.8),2),digits=3))
print(format(combineasym(c(1.5,1.2),c(.5,.8),2),digits=3))
print(format(combineasym(c(1.5,1.5),c(.5,.5),model=2),digits=3))
