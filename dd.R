rm(list=ls())

# 1 -----------------------------------------------------------------------
##a
fun1 = function(x){ out=x^3-2*x-5; return(out)}
bi.fun = function(a,b,x,iter.max){
  for(i in 1:iter.max){
    if(fun1(a)*fun1(x)>0){a=x;b=b}else{a=a;b=x}
    nx=(a+b)/2
    if(abs(nx-x)<eps){return(nx)}else{x=nx}
  }
}

#interval
a=0;b=3;eps=1e-7;iter.max=1e+3 
a=0;b=1.5;eps=1e-7;iter.max=1e+3
a=0.5;b=1;eps=1e-7;iter.max=1e+3
x=(a+b)/2
bi.fun(a,b,x,iter.max)

##b
x = seq(0,3,by=0.05)
plot(x,fun1(x))

a=0.5;b=1;eps=1e-7;iter.max=1e+3
x=(a+b)/2
for(i in 1:iter.max){
  print(i) #iteration step
  if(fun1(a)*fun1(x)>0){a=x;b=b}else{a=a;b=x}
  nx=(a+b)/2
  if(abs(nx-x)<eps){return(nx)}else{x=nx}
  abline(v=nx)
}

# 2 -----------------------------------------------------------------------
##a
fun2=function(x){ sin(x)-((x^2)/10) }
a = -5; b = 5; num = 100
gd.fun = function(a,b,num){
  x.vec = seq(a,b,length.out=num)
  f.vec = rep(0,num)
  for(x.pos in 1:length(x.vec)){
    f.vec[x.pos] = fun2(x.vec[x.pos])
  }
  min.pos = which.min(f.vec)
  max.pos = which.max(f.vec)
  x.min = x.vec[min.pos]
  f.min = f.vec[min.pos]
  x.max = x.vec[max.pos]
  f.max = f.vec[max.pos]
  return(list(x.min=x.min,f.min=f.min,x.max=x.max,f.max=f.max))
}
gd.fun(a,b,num)

##b
x = seq(-5,5,length.out=num)
plot(x,fun2(x))
abline(v=gd.fun(a=-5,b=5,num=100)$x.min,lty=2)
abline(v=gd.fun(a=-5,b=5,num=100)$x.max,lty=2)

# 3 -----------------------------------------------------------------------
n=100;
w.vec = seq(0,4,length=n)
x.vec = rnorm(n,mean=0,sd=1)

##a
f.fun3a = function(x.vec,w.vec,mu) {sum((w.vec)*abs(x.vec-mu))}

golden = function(f,dn,up,eps=1e-7){
  golden.ratio = 2/(sqrt(5)+1)
  mu1 = dn + golden.ratio*(up-dn)
  mu2 = up - golden.ratio*(up-dn)
  f1 = f(x.vec,w.vec,mu1)
  f2 = f(x.vec,w.vec,mu2)
  
while(abs(up-dn)>eps){
  if(f2>f1){
    dn = mu2; mu2 = mu1; f2 = f1
    mu1 = dn+golden.ratio*(up-dn)
    f1 = f(x.vec,w.vec,mu1)
  }
  else{
    up = mu1; mu1 = mu2; f1 = f2
    mu2 = up-golden.ratio*(up-dn)
    f2 = f(x.vec,w.vec,mu2)
  }
} 
  (up+dn)/2
}

golden(f.fun3a,-5,5,eps=1e-5)

##b
f.fun3b = function(x.vec,w.vec,mu){sum((w.vec)*abs(x.vec-mu))}
golden(f.fun3b,-5,5,eps=1e-5)

##c
f.fun3c = function(alp,x.vec,w.vec,mu){
  alp*f.fun3a(x.vec,w.vec,mu)+(1-alp)*f.fun3b(x.vec,w.vec,mu)}
golden = function(f,dn,up,eps=1e-7,alp){
  golden.ratio = 2/(sqrt(5)+1)
  mu1 = dn + golden.ratio*(up-dn)
  mu2 = up - golden.ratio*(up-dn)
  f1 = f.fun3c(alp,x.vec,w.vec,mu1)
  f2 = f.fun3c(alp,x.vec,w.vec,mu2)
  
  while(abs(up-dn)>eps){
    if(f2>f1){
      dn = mu2; mu2 = mu1; f2 = f1
      mu1 = dn+golden.ratio*(up-dn)
      f1 = f.fun3c(alp,x.vec,w.vec,mu1)
    }
    else{
      up = mu1; mu1 = mu2; f1 = f2
      mu2 = up-golden.ratio*(up-dn)
      f2 = f.fun3c(alp,x.vec,w.vec,mu2)
    }
  } 
  (dn+up)/2
}
golden(f.fun3C,-5,5,eps=1e-5,alp=2)

# 4 -----------------------------------------------------------------------
rm(list=ls())

loss.po.fun = function(y.vec,x.mat,b.vec) {
  xb.vec = drop(x.mat %*% b.vec)+1
  ret = -sum(y.vec*xb.vec) + sum(exp(xb.vec))
  return(ret)
}
grad.po.fun = function(y.vec, x.mat, b.vec) {
  exb.vec = exp(drop(x.mat %*% b.vec)+1)
  ret = t(x.mat) %*% (-y.vec + exb.vec)
  return(ret)
}
hess.po.fun = function(y.vec, x.mat, b.vec) {
  exb.vec = exp(drop(x.mat %*% b.vec)+1)
  ret = t(x.mat) %*% diag(exb.vec) %*% x.mat 
  return(ret)
}

set.seed(2020)
n = 30; p = 4
x.mat = matrix(rnorm(n*p),ncol=p,nrow=n)
tb.vec = b.vec = 1/(1:p)
mu.vec = exp((x.mat%*%tb.vec)+1)
y.vec = rpois(n=n,lambda=mu.vec)

##a
nb1.vec = nb2.vec = b.vec
gr = 2/(sqrt(5)+1)
eps = 1e-3; iter.max = 1e+5

cd.fun = function(y.vec,x.mat,b.vec,eps,iter.max){
  
  for(iter in 1:iter.max){
    ob.vec = b.vec
    for(j in 1:p){
      a = -10; b = 10
      c = gr*a + (1-gr)*b
      d = gr*b + (1-gr)*a
      nb1.vec[j] = c; nb2.vec[j] = d
      f1 = loss.po.fun(y.vec,x.mat,nb1.vec)
      f2 = loss.po.fun(y.vec,x.mat,nb2.vec)
      
      # golden 
      while(abs(b-a)>eps){
        if(f1>f2){
          a = c;
          c = gr*a + (1-gr)*b
          d = gr*b + (1-gr)*a
          nb1.vec[j] = c; nb2.vec[j] = d
          f1 = loss.po.fun(y.vec,x.mat,nb1.vec)
          f2 = loss.po.fun(y.vec,x.mat,nb2.vec)
        }
        else{
          b = d
          c = gr*a + (1-gr)*b
          d = gr*b + (1-gr)*a
          nb1.vec[j] = c; nb2.vec[j] = d
          f1 = loss.po.fun(y.vec,x.mat,nb1.vec)
          f2 = loss.po.fun(y.vec,x.mat,nb2.vec)
        }
      } 
      b.vec[j] = (a+b)/2
    }
    if(sum(abs(b.vec-ob.vec))<eps) break
    if(iter>=iter.max) break
  }
  return(b.vec)
}

cd.fun(y.vec,x.mat,b.vec,eps,iter.max)

##b
gd.fun = function(y.vec,x.mat,b.vec,eps,iter.max){
  loss.vec = rep(0,iter.max)
  for(iter in 1:iter.max){
    loss.vec[iter] = loss.po.fun(y.vec,x.mat,b.vec)
    alp = 2
    for(i in 1:100){
      d = grad.po.fun(y.vec,x.mat,b.vec)/sum(abs(grad.po.fun(y.vec,x.mat,b.vec)))
      nb.vec = b.vec - alp*d
      if(loss.po.fun(y.vec,x.mat,nb.vec)-loss.po.fun(y.vec,x.mat,b.vec)< -eps) break
      alp = alp*0.6      
    }
    if(sum(abs(nb.vec-b.vec))<eps) break
    b.vec = nb.vec
  }
  return(drop(b.vec))
}
gd.fun(y.vec,x.mat,b.vec,eps=1e-5,iter.max=1e+3)
  
##c
nr.fun=function(y.vec,x.mat,b.vec,eps,iter.max){
  for(iter in 1:iter.max){
    nb.vec = b.vec-drop(solve(hess.po.fun(y.vec,x.mat,b.vec))%*%
                          grad.po.fun(y.vec,x.mat,b.vec))
    if(sum(abs(nb.vec-b.vec))<eps) break
    b.vec = nb.vec
  }
  return(b.vec)
}
nr.fun(y.vec,x.mat,b.vec,eps=1e-5,iter.max=1e+3)

