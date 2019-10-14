rm(list=ls())
library(nlme)
library(IPEC)

# Massa fresca do fruto 
massa<- c(0.55,2.94,9.99,19.71,39.20,49.94,64.46,79.22,82.09,93.12,96.02,97.65)
daa <- c(40,60,80,100,120,140,160,180,200,220,240,260)

log <- gnls (massa~alpha/(1+ exp(k*(gamma-daa))), star=list(alpha = 100,gamma=59, k=0.04))
summary(log)
coef(log)
intervals(log)
confint.default(log)

library(nlstools)
log1 <- nls (massa~alpha/(1+ exp(k*(gamma-daa))), star=list(alpha = 100,gamma=59, k=0.04))
confint2(log1, level = 0.95, method = c("asymptotic"))
confint2(log1, level = 0.95, method = c( "profile"))



# Medidas de acordo com o SAS 

logistic <- function(theta, X){
  x    <- X
  alpha <- theta[1]
  gamma <- theta[2]
  k <- theta[3]
  alpha/(1+ exp(k*(gamma-x)))
}

ini.val <- c(100, 59, 0.04)

ajustando       <- fitIPEC( logistic, x=daa, y=massa, ini.val=ini.val, xlim=NULL, ylim=NULL,fig.opt=FALSE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
parA     <- ajustando$par
parA
viciodebox<- biasIPEC( logistic, theta=parA, x=daa, y=massa, tol=1e-20 ) 
viciodebox

curvatura <- curvIPEC(logistic, theta=parA, x=daa, y=massa, alpha=0.05, method.args=list(eps=1e-4, d=0.11, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2)) 
curvatura

assimetria <- skewIPEC( logistic, theta=parA, x=daa, y=massa, tol= 1e-20 )
assimetria


The_NLIN_Procedure <- as.data.frame(matrix(NA, nrow=3, ncol=7))

rownames(The_NLIN_Procedure) <- c("alpha","gamma","k")
colnames(The_NLIN_Procedure) <- c("Estimate", "Approx_Std_Error","Confidence Limits LI","Confidence Limits LS", "Skewness", "Bias","PercentBias")
The_NLIN_Procedure[,"Estimate"] <- c(round(coef(log),4))
The_NLIN_Procedure[,"Approx_Std_Error"] <- c(sqrt(diag(vcov(log))))
The_NLIN_Procedure[,"Confidence Limits LI"] <- c(intervals(log)[1]$coef[,1])
The_NLIN_Procedure[,"Confidence Limits LS"] <- c(intervals(log)[1]$coef[,3])

The_NLIN_Procedure[,"Skewness"] <- c(assimetria)

The_NLIN_Procedure[,"Bias"] <- c(viciodebox$bias)
The_NLIN_Procedure[,"PercentBias"] <- c(viciodebox$percent.bias)


The_NLIN_Procedure


# Global Nonlinearity Measures

Global_Nonlinearity_Measures<- as.data.frame(matrix(NA, nrow=3, ncol=1))
rownames(Global_Nonlinearity_Measures) <- c("RMS Intrinsic Curvature","RMS Parameter-Effects Curvature","Curvature Critical Value")
colnames(Global_Nonlinearity_Measures)<-c("")
Global_Nonlinearity_Measures["RMS Intrinsic Curvature",]<-c(curvatura$rms.ic)

Global_Nonlinearity_Measures["RMS Parameter-Effects Curvature",]<-c(curvatura$rms.pec)


Global_Nonlinearity_Measures["Curvature Critical Value",]<-c(curvatura$critical.c)

Global_Nonlinearity_Measures


# Gompertz


Gompertz <- function(theta, X){
  x    <- X
  alpha <- theta[1]
  gamma <- theta[2]
  k <- theta[3]
  alpha*exp(-exp(k*(gamma-X)))
}

ini.val <- c(100, 48, 0.04)

ajustando       <- fitIPEC( Gompertz, x=daa, y=massa, ini.val=ini.val,fig.opt=TRUE )
parA     <- ajustando$par
parA
viciodebox<- biasIPEC( Gompertz, theta=parA, x=daa, y=massa, tol=1e-20 ) 
viciodebox

curvatura <- curvIPEC(Gompertz, theta=parA, x=daa, y=massa, alpha=0.05, method.args=list(eps=1e-4, d=0.11, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2)) 
curvatura

assimetria <- skewIPEC( Gompertz, theta=parA, x=daa, y=massa, tol= 1e-20 )
assimetria


The_NLIN_Procedure <- as.data.frame(matrix(NA, nrow=3, ncol=7))

rownames(The_NLIN_Procedure) <- c("alpha","gamma","k")
colnames(The_NLIN_Procedure) <- c("Estimate", "Approx_Std_Error","Confidence Limits LI","Confidence Limits LS", "Skewness", "Bias","PercentBias")
The_NLIN_Procedure[,"Estimate"] <- c(round(coef(log),4))
The_NLIN_Procedure[,"Approx_Std_Error"] <- c(sqrt(diag(vcov(log))))
The_NLIN_Procedure[,"Confidence Limits LI"] <- c(intervals(log)[1]$coef[,1])
The_NLIN_Procedure[,"Confidence Limits LS"] <- c(intervals(log)[1]$coef[,3])

The_NLIN_Procedure[,"Skewness"] <- c(assimetria)

The_NLIN_Procedure[,"Bias"] <- c(viciodebox$bias)
The_NLIN_Procedure[,"PercentBias"] <- c(viciodebox$percent.bias)


The_NLIN_Procedure


# Global Nonlinearity Measures

Global_Nonlinearity_Measures<- as.data.frame(matrix(NA, nrow=3, ncol=1))
rownames(Global_Nonlinearity_Measures) <- c("RMS Intrinsic Curvature","RMS Parameter-Effects Curvature","Curvature Critical Value")
colnames(Global_Nonlinearity_Measures)<-c("")
Global_Nonlinearity_Measures["RMS Intrinsic Curvature",]<-c(curvatura$rms.ic)

Global_Nonlinearity_Measures["RMS Parameter-Effects Curvature",]<-c(curvatura$rms.pec)


Global_Nonlinearity_Measures["Curvature Critical Value",]<-c(curvatura$critical.c)

Global_Nonlinearity_Measures

# Como obter a do R 

n<- length(massa)
p <- length(coef(log))
qf(0.95,p,n-p)
curvatura$rms.ic*sqrt(qf(0.95,p, n-p))
curvatura$rms.pec*sqrt(qf(0.95,p, n-p))
1/(sqrt(qf(0.95,p, n-p)))
#1/(2*sqrt(qf(0.95,p, n-p)))
2*(1/(sqrt(qf(0.95,p,n-p))))




alpha<- coef(log)[1]
gamma<- coef(log)[2]
k<- coef(log)[3]
x <- daa

m<-cbind(daa,coef(log))

dfunção<- deriv(~alpha/(1+ exp(k*(gamma-x))),"x")
x=daa
plot(x,attr(eval(dfunção),"gradient"))


library(ggpubr)

library(ggExtra)

dados <- data.frame(massa,daa)
p <- ggscatter(dados, x = "daa", y = "massa")
ggMarginal(p, type = "boxplot",margins = c("y"))


D(expression(alpha/(1+ exp(k*(gamma-x)))), "x")
tcc<- alpha * (exp(k * (gamma - x)) * k)/(1 + exp(k * (gamma - x)))^2
cbind(x,tcc)

par(mar = c(5,5,3,5))
plot(massa~daa, type="n", xlab="Dias após emergência", ylab="Massa em (g)")
lines(daa, fitted(log))
par(new = T)
plot(x,tcc, type="l",axes=F, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4)
mtext(side = 4, line = 2, 'Taxa de crescimento')


par(mar = c(5,5,3,5))
plot(massa~daa, type="n", xlab="Dias após emergência", ylab="Massa em (g)", axes=F)
axis(1, at=c(0,daa, 280))
axis(2, at=c(-10,seq(0.55,105, 5)))
lines(daa, fitted(log), lwd="2", col="red")
abline(v=seq(0,280,10),h=seq(0,280,10),lty=2, col="gray")
par(new = T)
plot(x,tcc, type="b",axes=F, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4, at=c(-10,seq(0.05,0.90,0.05)))
mtext(side = 4, line = 2, 'Taxa de crescimento')
legend("topleft",
       legend=c("Logistico", "TCC"),
       lty=c(2,1), col=c("red3", "black"), bty="n")


