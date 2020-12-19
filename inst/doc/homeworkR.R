## -----------------------------------------------------------------------------
library(lattice)
n <- seq(20, 50, 10)
x <- rnorm(sum(n))
y <- factor(rep(n, n), labels=paste("n =", n))
densityplot(~ x | y,
panel = function(x, ...) {
panel.densityplot(x, col="red", ...)
panel.mathdensity(dmath=dnorm,
args=list(mean=mean(x), sd=sd(x)),
col="darkblue")
})

## -----------------------------------------------------------------------------
x1<-c(0.10, 0.11, 0.12, 0.13, 0.14, 0.15,
0.16, 0.17, 0.18, 0.20)
x2<-c(1,2,3,4,5,5,4,3,2,1)
y<-c(42.0, 43.5, 49.0, 45.5, 45.0, 47.5,
49.0, 53.0, 53.0, 55.0)
lm.sol<-lm(y ~ x1+x2)
knitr::kable(summary(lm.sol)$coef)

## -----------------------------------------------------------------------------
f1=function(n,a,b){
#n: sample size; Pareto(a,b)
u=runif(n)
x=b*(1-u)^(-1/a) 
hist(x, prob = TRUE, main = expression(f(x)==(a*b^a)/x^(a+1)))
y=seq(2, 100, .01)
lines(y, (a*b^a)/y^(a+1),col='red')
}
f1(1000,2,2)


## -----------------------------------------------------------------------------
f2=function(n){
#n:sample size
k=0;
y=numeric(n);
while(k<n){
  u1=runif(1,-1,1)
  u2=runif(1,-1,1)
  u3=runif(1,-1,1)
  k=k+1
  if ((abs(u3)>=abs(u2))&&(abs(u3)>=abs(u1))){
    y[k]=u2;}
  else{
    y[k]=u3;}
}
hist(y,prob = TRUE, main = expression(f(y)==-0.25*x^3+0.75*x+0.5))
}
f2(5000)

## -----------------------------------------------------------------------------
n=1000;
r=4;
beta=2;
lameda=rgamma(n,r,beta);
x=rexp(n,lameda)
hist(x, prob = TRUE)
y=seq(0, 15, .01)
lines(y, r*beta^r*(beta+y)^(-r-1),col='red')

## -----------------------------------------------------------------------------
sample_size=10000;
t=runif(sample_size,0,pi/3)#generate 10000 number from U(0,pi/3)
int_hat=mean(sin(t))*pi/3 #estimate the value of the integral
int_exact=0.5;
print(c(int_hat,int_exact))#outout the estimate and the exact value
print(abs(int_hat-int_exact))#the absolute value of the difference 

## -----------------------------------------------------------------------------
sample_size=10000;
t=runif(sample_size/2,0,1)#generate 5000 number from U(0,1) for antithetic variate approach
v1=runif(sample_size/2,0,1)#generate another 5000 number from U(0,1) for simple MC
u1=c(t,v1)#random number for simple MC method
e=exp(1)
P1=e^u1
P2=(e^u1+e^(1-u1))/2
int_hat1=mean(P1) #estimate the value of the integral with simple MC method
int_hat2=mean(P2)
int_exact=e-1;
print(c(int_hat1,int_hat2,int_exact))#output the two estimate and the exact value
print(c(var(P1),var(P2)))
var1=-0.5*e^2+2*e-1.5
var2=(-3*e^2+10*e-5)/4
reduction=(var1-var2)/var1
reductionhat=(var(P1)-var(P2))/var(P1)
print(c(var1,var2))#Var(e^u),Var((e^u+e^(1-u))/2)
print(c(reduction,reductionhat))

## -----------------------------------------------------------------------------
#plot the function graph for comparison
g=function(x) (x^2/(sqrt(2*pi)))*exp(-0.5*x^2)
f1=function(x) 0.5*exp(-0.5*(x-1))
f2=function(x) (4/pi)/(1+x^2)
plot(g,xlim=c(1,10),ylim=c(0,0.8),col='black')
plot(f1,xlim=c(1,10),ylim=c(0,0.8),col='blue',add=T)
plot(f2,xlim=c(1,10),ylim=c(0,0.8),col='red',add=T)


## -----------------------------------------------------------------------------
#x~f Y=X+1~exp(0.5)
#thus we just need to generate number from exp(0.5)
n=10000;
ratio=numeric(n);
for (i in 1:n){
  u=rexp(1,0.5)+1;
  ratio[i]=g(u)/f1(u);
}
result1=mean(ratio)
print(result1)#the result
print(var(ratio))#the variance using f1

## -----------------------------------------------------------------------------
n=10000;
ratio=numeric(n);
for (i in 1:n){
u0=runif(1);
u=tan((pi/4)*(u0+1))##inverse transform method
ratio[i]=g(u)/f2(u);
}
result2=mean(ratio)
print(result2)#the result
print(var(ratio))#the variance using f2

## -----------------------------------------------------------------------------
M=1000; k=5 #5 layers
r=M/k
N=50
T=numeric(k)
g=function(x) exp(-x)/(1+x^2)#def g(x)
c=function(j) 1/(exp(-(j-1)/5)-exp(-j/5)) #def the coefficient c, which changes with j
ff=function(x,j) 1/(exp(-(j-1)/5)-exp(-j/5))*exp(-x) #def f(x)
est=numeric(N)
for(i in 1:N){#repeat calculation for N times
  for(j in 1:k){#use important function 
    u=runif(r)
    x0=(j-1)/k
    x=-log(exp(-x0)-u/c(j))
    T[j]=mean(g(x)/ff(x,j))}
  est[i]=sum(T)
}
print(mean(est))
print(sd(est))# estimate the sd of the result

## -----------------------------------------------------------------------------
set.seed(54321)
n=1000
N=1000#repeat the process to estimate the confidence level
k=0;
for (i in 1:N){
X=rlnorm(n,3,3)#random number from LN(3,3^2)
t=qt(0.975,n-1)
left=mean(log(X))-t*sd(log(X))/sqrt(n)
right=mean(log(X))+t*sd(log(X))/sqrt(n)
if((left<=3)&&(right>=3))
    k=k+1
}
print(k/N)#the confidence level

## -----------------------------------------------------------------------------
set.seed(54321)
n=20
N=1000#repeat times
k=0
af=0.05#alpha
left=numeric(N);
right=numeric(N);
for (i in 1:N){
X=rchisq(n,2)
t=qt(0.975,n-1)
left[i]=mean(X)-t*sd(X)/sqrt(n)
right[i]=mean(X)+t*sd(X)/sqrt(n)
if((left[i]<=2)&&(right[i]>=2))#Monte Carlo experiment to estimate the coverage probability
    k=k+1
}
print(k/N) # the coverage probability
print(sd(left))
print(sd(right))# the sd of the CI

## -----------------------------------------------------------------------------
set.seed(54321)
n=20
N=1000
a=0.05
sigma=2
UCL=numeric(N)
for (i in 1:N){
Y=rnorm(n, mean=0, sd=2)
UCL[i]=(n-1)*var(Y)/qchisq(a,n-1)}
s=sum(sigma<=UCL)
print(s/N)
print(sd(UCL))

## -----------------------------------------------------------------------------
rm(list = ls())
#define a function to compute the sample skewness coeff
sk=function(x) {
xx=mean(x)
e3=mean((x - xx)^3)
e2=mean((x - xx)^2)
return( e3 / e2^1.5 )
}
n=c(20,100,500);
#different sample size
l=length(n)
m=10000
a=2#beta distribition beta(a,a)
nu=5#t-distribution t(5)
alpha=0.05#confidence level
w=qnorm(.975, 0, sqrt(6/n))#alternative hypothesis
test1=test2=numeric(m)

result.beta=result.t=numeric(l)
for (i in 1:l){
  for (j in 1:m){
    x=rbeta(n[i],a,a)
    y=rt(n[i],nu)
    test1[j]= as.integer(abs(sk(x)) >= w[i])
    test2[j]= as.integer(abs(sk(y)) >= w[i])
  }
result.beta[i]=mean(test1)#beta distribution
result.t[i]=mean(test2)#t distribution
}
print(result.beta)
print(result.t)

## -----------------------------------------------------------------------------
#count5test function
#reject H0:return 1; not reject H0:return 0
count5test <- function(a, b) {
A=a-mean(a)
B=b-mean(b)
px=sum(A > max(B)) + sum(A < min(B))
py=sum(B > max(A)) + sum(B < min(A))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(px, py)) > 5))
}

#comparison of count5test and F-test
sigma1=1
sigma2=1.5
ahat=0.055
m=10000
n=c(20,100,1000)
#test for small, medium, and large sample sizes
power1=power2=numeric(3)
#power1--count5test; power2--F-test
c=f=numeric(m)
for(i in 1:3){
  for(j in 1:m){
# generate samples under H1 to estimate power
x=rnorm(n[i], 0, sigma1)
y=rnorm(n[i], 0, sigma2)
c[j]=count5test(x, y)
f.sol=var.test(x,y,ratio=1,conf.level = ahat)#F-test
f[j]=(f.sol$p.value<ahat)#compare the p-value
}
power1[i]=mean(c) 
power2[i]=mean(f) 
}
print(power1)#count5test
print(power2)#F-test

## -----------------------------------------------------------------------------
#Example 6.8
library(MASS)
skw<-function(dt){
  n=nrow(dt)
  c=ncol(dt)
  central=dt
  for(i in 1:c){
    central[,i]<-dt[,i]-mean(dt[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(54867)
mu=c(0,0,0)
sigma=matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n=c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    dt=mvrnorm(n[i],mu,sigma) 
    skw(dt)
  }))
}
a

## -----------------------------------------------------------------------------
set.seed(54321)
msk=function(x){
n=nrow(x)#sample size of x and y
d=ncol(x)
sigma2=cov(x)#the MLE of the sample covariance
sigma=solve(sigma2)
xbar=apply(x,2,mean)#average
xx=x-xbar
sum=0
for (j in 1:n){
  for (i in 1:n){
sum=sum+(t(xx[i,])%*%sigma%*%(xx[j,]))^3
  }
}
sk=sum/n^2
return(sk)  
}

####main function
library(MASS)
n=30#sample size
m=100
alpha=0.1
d=2
eps=c(seq(0, .15, .01), seq(.15, 1, .05))
le=length(eps)
result=numeric(le)#power
temp=numeric(m)
x=matrix(nrow=n,ncol=d)
free=d*(d+1)*(d+2)/6#degrees of freedom
reject=6*qchisq(.975,free)/n
for (i in 1:le){
for (j in 1:m){
sigma=sample(c(1, 100), replace = TRUE,
size =n, prob = c(1-eps[i], eps[i]))#replacement sampling to produce random vector
for (k in 1:n)
{x[k,]=mvrnorm(1, numeric(d), sigma[k]*diag(d))}
temp[j]=as.integer(abs(msk(x)) >= reject)
}
result[i]=mean(temp) #save the power of different epsilon 
}
print(result)
######plot the figure#####
par(pin=c(1.5,1))
plot(eps, result, type = "b",
xlab = bquote(eps), ylim = c(0,1))
abline(h = .1, lty = 3)
sd <- sqrt(result * (1-result) / m) #figure out the sd
lines(eps, result+sd, lty = 3)
lines(eps, result-sd, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap) 
data1=as.matrix(law)
cor_sample=cor(law$LSAT, law$GPA)
print(cor_sample)#sample correlation
sample_size=nrow(data1)#sample size
cor_jack=numeric(sample_size)
for (i in 1:sample_size){
cor_jack[i]=cor(data1[-i,1],data1[-i,2])  
#calculate each correlation  
}
jbias=(sample_size-1)*(mean(cor_jack)-cor_sample)
jse=sqrt((sample_size-1)*mean((cor_jack-cor_sample)^2))
round(c(bias_jack=jbias,se_jack=jse),4)

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
data2=as.matrix(aircondit)
#define a function to calculate the mean of the sample
bootbar=function(x,i){
  mean(x[i])
} 
size=nrow(data2)#sample size
result=boot(data2, statistic =bootbar, R=1000)
result
#figure out different bootstrap intervals
boot.ci(result, type=c("norm","basic","perc", "bca"))

## -----------------------------------------------------------------------------
library(bootstrap)
Sigmahat=cov(scor);
n=nrow(scor)#sample size
ev1=eigen(Sigmahat)
thetahat=ev1$values[1]/sum(ev1$values)
round(c(thetahat=thetahat),4)
#jackknife
thetajack=numeric(n)
for (i in 1:n){
Sigmajack=cov(scor[-i,])#calculate the sample cov
evjack=eigen(Sigmajack)
thetajack[i]=evjack$values[1]/sum(evjack$values)
biasjack=(n-1)*(mean(thetajack)-thetahat)}#jackknife bias
sejack=sqrt((n-1)*mean((thetajack-mean(thetajack))^2))#jackknife se
round(c(jackknife_bias=biasjack,jackknife_se=sejack),4)

## -----------------------------------------------------------------------------
library(DAAG)
library(lattice)
dt=ironslag
#leave-two-out cross validation
n=length(dt$magnetic)
err1=err2=err3=err4=matrix(nrow=n*(n-1)/2,ncol=2)
k=1;#save error
for (i in 1:(n-1)){#choose the first test sample
  for (j in (i+1):n){#choose the second test sample
    sp1=dt$magnetic[c(-i,-j)]
    sp2=dt$chemical[c(-i,-j)]
    #Linear
    J1 <- lm(sp1 ~ sp2)
   j11 <- J1$coef[1] + J1$coef[2] * dt$chemical[i]
   j12 <- J1$coef[1] + J1$coef[2] * dt$chemical[j]
    err1[k,1] <- dt$magnetic[i] - j11
    err1[k,2] <- dt$magnetic[j] - j12
    # Quadratic
    J2 <- lm(sp1 ~ sp2 + I(sp2^2))
    j21 <- J2$coef[1] + J2$coef[2] * dt$chemical[i] +
      J2$coef[3] * dt$chemical[i]^2
    j22 <- J2$coef[1] + J2$coef[2] * dt$chemical[j] +
      J2$coef[3] * dt$chemical[j]^2
    err2[k,1] <- dt$magnetic[i] - j21
    err2[k,2] <- dt$magnetic[j] - j22
    #Exponential
    J3 <- lm(log(sp1) ~ sp2)
    j31 <- J3$coef[1] + J3$coef[2] * dt$chemical[i]
    j31 <- exp(j31)
    j32 <- J3$coef[1] + J3$coef[2] * dt$chemical[j]
    j32 <- exp(j32)
    err3[k,1] <- dt$magnetic[i] - j31
    err3[k,2] <- dt$magnetic[j] - j32
    #Log-Log
    J4 <- lm(log(sp1) ~ log(sp2))
    j41 <- J4$coef[1] + J4$coef[2] * log(dt$chemical[i])
    j41 <- exp(j41)
    j42 <- J4$coef[1] + J4$coef[2] * log(dt$chemical[j])
    j42 <- exp(j42)
    err4[k,1] <- dt$magnetic[i] - j41
    err4[k,2] <- dt$magnetic[j] - j42
    
    k=k+1
  }
}

#estimate the prediction error
c(mean(err1^2), mean(err2^2), mean(err3^2), mean(err4^2))

## -----------------------------------------------------------------------------
count5test = function(x, y) {
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y))
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}

# Count Five test permutation
c5p= function(z) {
n = length(z)
x = z[1:(n/2)]
y = z[-(1:(n/2))]
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y)) 
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0) 
return(as.integer(max(c(outx, outy)) > 5))
}
permutation = function(z,R) {
  n = length(z)
  out = numeric(R)
  for (r in 1: R){
      p = sample(1:n ,n ,replace = FALSE)
      out[r] = c5p(z[p])
  }
  sum(out)/R
}              
#main
n1 = 20
n2 = 50
mu1 = mu2 = 0
sigma1 = sigma2 = 1
m = 1e3
alphahat1 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean
y = y - mean(y)
count5test(x, y)
}))
alphahat2 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
permutation(z,1000) 
})<0.05)

round(c(count5test=alphahat1,count5test_permutation=alphahat2),4)

## -----------------------------------------------------------------------------
library(Ball)
library(RANN)
library(energy)
library(boot)
#def stat function
stat=function(z, ix, sizes,k) {
  n1=sizes[1]; n2=sizes[2]; n=n1 + n2
  if(is.vector(z)) z=data.frame(z,0);
  z=z[ix, ];
  NN=nn2(data=z, k=k+1) 
  tp1=NN$nn.idx[1:n1,-1]
  tp2=NN$nn.idx[(n1+1):n,-1]
  i1=sum(tp1 < n1 + .5); i2= sum(tp2 > n1+.5)
  (i1 + i2) / (k * n)
} 
#def knn funcion
eqdist.nn=function(z,sizes,k,R){
  boot.obj=boot(data=z,statistic=stat,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value=mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
} 
#####simulation##### 
m=200;
knn=3; 
di=2; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*di,0,1.5),ncol=di);
  y=cbind(rnorm(n2),rnorm(n2));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)

## -----------------------------------------------------------------------------
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*di,0,2),ncol=di);
  y=cbind(rnorm(n2));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)

## -----------------------------------------------------------------------------
#####simulation##### 
m=200;
knn=3; 
di=2; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*di,-0.5,1.5),ncol=di);
  y=cbind(rnorm(n2,0,1),rnorm(n2,0,1));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)

## -----------------------------------------------------------------------------
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*di,-0.5,2),ncol=di);
  y=cbind(rnorm(n2,0,1));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)


## -----------------------------------------------------------------------------
###########t-distribution###########
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rt(n1*di,df=1),ncol=di);
  y=cbind(rt(n2,df=10));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3) 

## -----------------------------------------------------------------------------
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rt(n1*di,df=1),ncol=di);
  y=cbind(rnorm(n2));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3) 


## -----------------------------------------------------------------------------
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999;
#bimodel distribution
mu1=c(-1,0)
sig1=c(1.5,1)
mu2=c(-0.5,1)
sig2=c(4,2)
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  lab1=sample(c(1,2), replace = TRUE,
              size = n1, prob = c(0.5,0.5))
  lab2=sample(c(1,2), replace = TRUE,
              size = n2, prob = c(0.5,0.5))
  x=matrix(rnorm(n1*di,mu1[lab1],sig1[lab1]),ncol=di);
  y=cbind(rnorm(n2,mu2[lab2],sig1[lab2]));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3) 


## -----------------------------------------------------------------------------
m=100;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=30
n2=30; 
R=999;
mu1=c(-1,0)
sig1=c(1.5,1)
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  lab1=sample(c(1,2), replace = TRUE,
              size = n1, prob = c(0.3,0.7))
  x=matrix(rnorm(n1*di,mu1[lab1],sig1[lab1]),ncol=di);
  y=cbind(rt(n2,1));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)

## -----------------------------------------------------------------------------
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=100
n2=10; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*di,-0.5,2),ncol=di);
  y=cbind(rnorm(n2,0,1));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)


## -----------------------------------------------------------------------------
###########t-distribution###########
m=200;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=100
n2=10; 
R=999; 
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rt(n1*di,df=1),ncol=di);
  y=cbind(rt(n2,df=10));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3) 

## -----------------------------------------------------------------------------
m=100;
knn=3; 
di=1; #dimension
set.seed(31415)
n1=100
n2=10; 
R=999;
mu1=c(-1,0)
sig1=c(1.5,1)
mu2=c(-0.5,1)
sig2=c(4,2)
n=n1+n2; 
N=c(n1,n2)
#save p-value
vp=matrix(NA,m,3)
for(i in 1:m){
  lab1=sample(c(1,2), replace = TRUE,
              size = n1, prob = c(0.3,0.7))
  lab2=sample(c(1,2), replace = TRUE,
              size = n2, prob = c(0.2,0.8))
  x=matrix(rnorm(n1*di,mu1[lab1],sig1[lab1]),ncol=di);
  y=cbind(rnorm(n2,mu2[lab2],sig1[lab2]));
  z=rbind(x,y)
  ##knn, erergy, ball
  vp[i,1]=eqdist.nn(z,N,knn,R)$p.value
  vp[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  vp[i,3]= bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
} 
alpha=0.1;
pwr1=colMeans(vp<alpha) #calculate power under 3 different tests
round(c(nn=pwr1[1],energy=pwr1[2],
        ball=pwr1[3]),3)

## -----------------------------------------------------------------------------
#define a function to calculate the pdf of Laplace distribution
La=function(x){
  0.5*exp(-abs(x))
}
#define a function to generate the chain
#x0: initial value;
#N: the length of the chain
LRWM=function(sigma, x0, N) { 
x=numeric(N)
x[1]=x0
u=runif(N)
acp=0 #acceptance rate
for (i in 2:N) { 
  y=rnorm(1, x[i-1], sigma) 
  if (u[i] <= (La(y) / La(x[i-1])))
  {x[i]=y
   acp=acp+1 
   } 
  else {
    x[i]=x[i-1] 
  } 
  }
return(list(x=x, acp=acp)) 
}
N=2000
sigma=c(0.05,0.5,1,2,9,16)
x0=0 
rw1=LRWM(sigma[1], x0, N) 
rw2=LRWM(sigma[2], x0, N) 
rw3=LRWM(sigma[3], x0, N) 
rw4=LRWM(sigma[4], x0, N)
rw5=LRWM(sigma[5], x0, N)
rw6=LRWM(sigma[6], x0, N)
round(c(rw1$acp/N, rw2$acp/N, rw3$acp/N, rw4$acp/N,rw5$acp/N,rw6$acp/N),4)
par(pin=c(1.5,0.75))
plot(1:2000,rw1$x,type='l',main="sigma=0.05",ylab="x")
plot(1:2000,rw2$x,type='l',main="sigma=0.5",ylab="x")
plot(1:2000,rw3$x,type='l',main="sigma=1",ylab="x")
plot(1:2000,rw4$x,type='l',main="sigma=2",ylab="x")
plot(1:2000,rw5$x,type='l',main="sigma=9",ylab="x")
plot(1:2000,rw6$x,type='l',main="sigma=16",ylab="x")

## -----------------------------------------------------------------------------
#define a function to calculate the diagnostic statstics
Gelman.Rubin=function(psi) { 
psi=as.matrix(psi)
n=ncol(psi)
k=nrow(psi)

psi.means=rowMeans(psi)
B=n * var(psi.means)
psi.w <- apply(psi, 1, "var")
W <- mean(psi.w)
v.hat=W*(n-1)/n + (B/n) 
r.hat=v.hat / W
return(r.hat) 
}
#choose sigma
sig=0.9
k=4
n=15000
b=1000
set.seed(12354)
#choose initial values
x0=c(-10, -5, 5, 10)
chain=matrix(nrow=k,ncol=n)
#generate chains
for (i in 1:k){
  chain[i,]=LRWM(sig,x0[i],n)$x
}
#calculate the statistics
psi=t(apply(chain, 1, cumsum)) 
for (i in 1:nrow(psi)) 
  psi[i,]= psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains 
par(pin=c(1.5,0.75))
for (i in 1:k) 
  plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat=rep(0, n)
for (j in (b+1):n) 
 rhat[j]=Gelman.Rubin(psi[,1:j]) 
par(pin=c(1.5,0.75))
plot(rhat[(b+1):n], type="l", xlab="", ylab="R") 
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
#choose sigma 
set.seed(12354)
sig=2
k=4
n=15000
b=1000
#choose initial values
x0=c(-10, -5, 5, 10)
chain=matrix(nrow=k,ncol=n)
#generate chains
for (i in 1:k){
  chain[i,]=LRWM(sig,x0[i],n)$x
}
#calculate the statistics
psi=t(apply(chain, 1, cumsum)) 
for (i in 1:nrow(psi)) 
  psi[i,]= psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains 
par(pin=c(1.5,0.75)) 
for (i in 1:k) 
  plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat=rep(0, n)
for (j in (b+1):n) 
  rhat[j]=Gelman.Rubin(psi[,1:j]) 
par(pin=c(1.5,0.75))
plot(rhat[(b+1):n], type="l", xlab="", ylab="R") 
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
k=c(4:25,100,500,1000)
S = function(a,k){
 ck = sqrt(a^2*k/(k+1-a^2))
 pt(ck,df=k,lower.tail=FALSE)
}

f = function(a,k){S(a,k)-S(a,k-1)}
#curve(f(x),xlim = c(0,sqrt(k)))
a <- seq(0, 4, by=0.01)
par(pin=c(1.5,0.75))
plot(a, f(a, k[23]), lty=1, col=1, type="l", xlim=c(0, 4), xlab="a", ylab="f(a|k)", main="f(a) with different k")
lines(a, f(a, k[24]), xlim = c(0, 4), lty=2, col=2)
lines(a, f(a, k[25]), xlim = c(0, 4), lty=3, col=3)

## -----------------------------------------------------------------------------
solve = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}
root = matrix(0,2,length(k))
for (i in 1:length(k)){
  root[2,i]=round(solve(k[i]),4)
}
root[1,] = k
rownames(root) = c('k','A(k)')
root

## -----------------------------------------------------------------------------
nA=444;
nB=132;
noo=361;
nAB=63;
#calculate the expectation of nAA
EnAA=function(p){
  nA*p^2/(p^2+2*p*r)
}
#calculate the expectation of nBB
EnBB=function(q){
  nB*q^2/(q^2+2*q*r)
}
omle=function(p,q){
nA*log(p*(p+2*r)) +nB*log(q*(q+2*r))+2*noo*log(r)+nAB*log(2*p*q) 
}
obj=function(x){#def the opposite MLE
  p=x[1];
  q=x[2];
  r=1-p-q;
  y=-(2*noo*log(r)+nAB*log(2*p*q)+nA*log(2*p*r)+nB*log(2*q*r)+(2*log(p)-log(2*p*r))*EnAA(pi)+(2*log(q)-log(2*q*r))*EnBB(qi))#we need to change the values of pi and qi in the algorithm
}
#set the initial p and q
#the initial value of p and q are saved in p[2] and q[2] 
p0=0.2;
q0=0.2;
R=20#maximum number of cycles
P=numeric(R+1);
Q=numeric(R+1);
OMLE=numeric(R);
P[1]=p0;
Q[1]=q0;
P[2]=p0+0.01;
Q[2]=q0+0.01;
i=2;
while ((abs(P[i]-P[i-1])>1e-6)||(abs(Q[i]-Q[i-1])>1e-6)){
  pi=P[i];
  qi=Q[i];
  r=1-pi-qi;
sol=optim(par=c(p0,q0),fn=obj) #solve the minimum
#save the MLE at each step
P[i+1]=sol$par[1]  
Q[i+1]=sol$par[2]
OMLE[i-1]=omle(P[i+1],Q[i+1])
i=i+1;
}
iter=i-1
P[3:i-1]
Q[3:i-1]
OMLE[1:iter-1]

## -----------------------------------------------------------------------------
formulas=list(
mpg ~ disp, 
mpg ~ I(1 / disp), 
mpg ~ disp + wt, 
mpg ~ I(1 / disp) + wt
)
n=length(formulas)
for (i in 1:n){
print(formulas[[i]])  
print(lm(formulas[[i]],data=mtcars))
}

## -----------------------------------------------------------------------------
lr=function(f){ 
  lm(f, data=mtcars)}
lapply(formulas,lr)

## -----------------------------------------------------------------------------
trials=replicate(100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE)
sapply(trials,function(pv) pv$p.value)

## -----------------------------------------------------------------------------
xs=list(runif(10), runif(5)) 
ws=list(rpois(10, 2) + 1, rpois(5, 2) + 1)
vapply(Map(weighted.mean,xs,ws),function(x) x,FUN.VALUE=c(wm1=0))

## -----------------------------------------------------------------------------
library(StatComp20087)

sigma=c(0.05,1,9,16)
x0=0
N=2000
rw1=LARMC(sigma[1], x0, N) 
rw2=LARMC(sigma[2], x0, N)
rw3=LARMC(sigma[3], x0, N) 
rw4=LARMC(sigma[4], x0, N)
par(pin=c(1.5,0.75));
plot(1:2000,rw1,type='l',main="sigma=0.05",ylab="x"); plot(1:2000,rw2,type='l',main="sigma=1",ylab="x"); plot(1:2000,rw3,type='l',main="sigma=9",ylab="x"); plot(1:2000,rw4,type='l',main="sigma=16",ylab="x")

## -----------------------------------------------------------------------------
La=function(x){
  0.5*exp(-abs(x))
}
#define a function to generate the chain
#x0: initial value;
#N: the length of the chain
LRWMR=function(sigma, x0, N) { 
x=numeric(N)
x[1]=x0
u=runif(N)
for (i in 2:N) { 
  y=rnorm(1, x[i-1], sigma) 
  if (u[i] <= (La(y) / La(x[i-1])))
  {x[i]=y
   } 
  else {
    x[i]=x[i-1] 
  } 
  }
return(x) 
}
par(pin=c(1.5,0.75))
sigma=0.05;
X1=LRWMR(sigma, x0, N);
X2=LARMC(sigma, x0, N);
qqplot(X1,X2,main="sigma=0.05")
sigma=1;
X3=LRWMR(sigma, x0, N);
X4=LARMC(sigma, x0, N);
qqplot(X3,X4,main="sigma=1")
sigma=9;
X5=LRWMR(sigma, x0, N);
X6=LARMC(sigma, x0, N);
qqplot(X5,X6,main="sigma=9")
sigma=16;
X7=LRWMR(sigma, x0, N);
X8=LARMC(sigma, x0, N);
qqplot(X7,X8,main="sigma=16")

## -----------------------------------------------------------------------------
library(microbenchmark)
t=microbenchmark(LRWMR(sigma, x0, N),LARMC(sigma, x0, N))
summary(t)

