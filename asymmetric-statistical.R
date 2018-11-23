#
#  combining results with asymmetric statistical errors 
#
# usage:  combine(results,sigmaplus,sigmaminus, useVariance, tolerance)
#
# first three arguments (required)  are vectors (must have same length) 
#     giving values, positive and negative errors
#
# 4th argument (optional)is boolean, switching between linear-sigma and linear-variance forms. Default is linear-sigma
#
# 5th argument (optional) is required accuracy of solution (as a fraction of overall range)). Default 1E-6. Unlikely to need to change this.

# the result is a self-describing list 

s=sprime=V=Vprime=0  # Globals. Their use saves parameters. 

combineds <- function(x,values){ # likelihood using linear sigma
y=rep(0,length.out=length(x))
for (i in 1:length(values)) y <- y-0.5*(x-values[i])**2/(s[i]+sprime[i]*(x-values[i]))**2 
return(y)
}
#
combinedv <- function(x,values){ # likelihood using linear variance
#              V can go negative, which must be protected against
#               but only far from the solution so it is just cosmetic
y=rep(0,length.out=length(x))
for (i in 1:length(values)){ 
  thisV=V[i]+Vprime[i]*(x-values[i])
  y[thisV>0]=y[thisV>0]-0.5*(x[thisV>0]-values[i])**2/thisV[thisV>0]
  y[thisV<=0]=-999 
}
return(y)
}
#
# two functions used for finding the Delta ln L =-.5 points
wanteds <- function (x,values,maxy) { return (combineds(x,values)+0.5-maxy)}
wantedv<- function (x,values,maxy) { return (combinedv(x,values)+0.5-maxy)}
#
combine <- function(results,sp,sm,useVariance=FALSE,TOL=1.E-6){  # the main function
#
  xlo=min(results-2*sm)
  xhi=max(results+2*sp)
  small=TOL*(xhi-xlo)    #  used for solving numerically 
  
# set up and draw empty picture. 
  x=seq(xlo,xhi,length.out=200)
  plot(0,0,type='n',xlim=c(xlo,xhi),ylim=c(-2,0),xlab='x',ylab='ln L')
  lines(c(xlo,xhi),c(0,0))
  lines(c(xlo,xhi),c(0,0)-.5)

# calcuate quantities from Equations 16 and 18 of the paper
  s <<- 2*sp*sm/(sp+sm)
  sprime <<- (sp-sm)/(sp+sm) 
  V <<- sm*sp
  Vprime <<- sp-sm
#
# draw a curve for each value
for (i in 1:length(results)){
  if(useVariance){
# linear V can go negative at large x. Does not affect results but can give crazy plots 
    thisV=V[i]+Vprime[i]*(x-results[i])
    y=-0.5*(x[thisV>0]-results[i])**2/thisV[thisV>0]
    lines(x[thisV>0],y,col='green')
  } else {
    y=-0.5*(x-results[i])**2/(s[i]+sprime[i]*(x-results[i]))**2
    lines(x,y,col='blue')
    }
  }
#
# iterate to find maximum of the sum
guess=mean(results)
test=100*small # ensure first iteration - R does not have do {} while()!
ntries=0 # collapse ungracefuly if exceeded
while(test>small){
  if(useVariance) {
     w=V/(V+Vprime*(guess-results))**2
     newguess=sum(w*(results-Vprime*(guess-results)**2*0.5/V)) /sum(w)
   } else {
     w=s/(s+sprime*(guess-results))**3
     newguess=sum(results*w)/sum(w)
   }
  test=abs(newguess-guess)
  guess=newguess
  ntries=ntries+1
  if(ntries > 100) stop("Problem - does not converge")
  }
#
lines(guess*c(1,1),c(-.05,.05),col='red') # draw solution on plot
answer=list()
answer[["mean"]]=guess
#
if(useVariance) {yyy=combinedv(x,results)} else {yyy=combineds(x,results)}
lines(x,yyy-max(yyy),col='red')  # draw the combined curve
answer[["chisquared"]] <- -2*max(yyy)
#
# now find the Delta ln L = -1/2 points using the R function  uniroot
if(useVariance) {u= uniroot(wantedv,c(guess,xhi),results,max(yyy))} else { u=uniroot(wanteds,c(guess,xhi),results,max(yyy))}
answer[["sigmaplus"]] <- u$root-guess
lines(u$root*c(1,1),c(-.45,-.55),col='red')
if(useVariance){u=uniroot(wantedv,c(xlo,guess),results,max(yyy))} else { u=uniroot(wanteds,c(xlo,guess),results,max(yyy))}
answer[["sigmaminus"]] <- guess- u$root
lines(u$root*c(1,1),c(-.45,-.55),col='red')
return(answer)
}
#
# here's an example combining three values
#
values <- c(1.9,2.4,3.1)
sigmam <- c(.5,.8,.4)
sigmap <- c(.7,.6,.5)
par(mfrow=c(1,2))
result <- combine(values,sigmap,sigmam)
print(result)
result <- combine(values,sigmap,sigmam,TRUE)
print(result)
#