rdata.hmm.Cartesian<-function(NUM, pii, A, f0, f1, f2)
{

theta<-rep(0, NUM)


x1<-rep(0, NUM)
x2<-rep(0, NUM)


## generating the states
 # initial state
theta[1]<-sample(1:4, 1, replace=FALSE, prob=pii)
 # other states
for (i in 2:NUM)
{
  if (theta[i-1]==1)
  {
     theta[i]<-sample(1:4, 1, replace=FALSE, prob=A[1, ])
  }
  else if (theta[i-1]==2)
  {
     theta[i]<-sample(1:4, 1, replace=FALSE, prob=A[2, ])
  }
  else if (theta[i-1]==3)
  {
     theta[i]<-sample(1:4, 1, replace=FALSE, prob=A[3, ])
  }
  else
  {
     theta[i]<-sample(1:4, 1, replace=FALSE, prob=A[4, ]) 
  }   
}

## generating the observations
for (i in 1:NUM)
{
  if (theta[i]==1)
  {
    x1[i]<-rnorm(1, mean=f0[1], sd=f0[2])
    x2[i]<-rnorm(1, mean=f0[1], sd=f0[2])   
  }
  else if (theta[i]==2)
  { 
    x1[i]<-rnorm(1, mean=f1[1], sd=f1[2])
    x2[i]<-rnorm(1, mean=f0[1], sd=f0[2])
  }
  else if (theta[i]==3)
  { 
    x1[i]<-rnorm(1, mean=f0[1], sd=f0[2])
    x2[i]<-rnorm(1, mean=f2[1], sd=f2[2])
  }
  else
  {
    x1[i]<-rnorm(1, mean=f1[1], sd=f1[2])
    x2[i]<-rnorm(1, mean=f2[1], sd=f2[2])
  }
}
data<-list(s=theta, x1=x1, x2=x2)
return (data)

}

