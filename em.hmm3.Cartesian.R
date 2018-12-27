em.hmm3.Cartesian<-function(x1, x2, x3, maxiter=200)
{

#####################################################################################

NUM<-length(x1)
# precision tolerance level
ptol<-1e-4
niter<-0

### initializing model parameters

pii.new<-c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125)

A.new<-matrix(c(c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125), 
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
                c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125)), 8, 8, byrow=TRUE)


f0<-c(0, 1)
f1.new<-c(2, 1)
f2.new<-c(2, 1)
f3.new<-c(2, 1)

diff<-10
Loglikelihood.new<--10000

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

pii.old<-pii.new
A.old<-A.new
f1.old<-f1.new
f2.old<-f2.new
f3.old<-f3.new
Loglikelihood.old<-Loglikelihood.new
## updating the weights and probabilities of hidden states

bwfw.res<-bwfw.hmm3.Cartesian(x1, x2, x3, pii.old, A.old, f0, f1.old, f2.old, f3.old)

# the backward-forward variable

alpha<-bwfw.res$bw
beta<-bwfw.res$fw


f0x1<-dnorm(x1, f0[1], f0[2])
f0x2<-dnorm(x2, f0[1], f0[2])
f0x3<-dnorm(x3, f0[1], f0[2])
f1x1<-dnorm(x1, f1.old[1], f1.old[2])
f1x2<-dnorm(x2, f2.old[1], f2.old[2])
f1x3<-dnorm(x3, f3.old[1], f3.old[2])

pdtheta<-array(rep(0, 64*(NUM-1)),c((NUM-1), 8, 8))
b1<-rep(0, NUM-1)
for(k in 1:(NUM-1))
{

######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 1, 1]<-alpha[k, 1, 1, 1]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[1, 1]
      pdtheta[k, 1, 2]<-alpha[k, 1, 1, 1]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[1, 2]
      pdtheta[k, 1, 3]<-alpha[k, 1, 1, 1]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[1, 3]
      pdtheta[k, 1, 4]<-alpha[k, 1, 1, 1]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[1, 4]
      pdtheta[k, 1, 5]<-alpha[k, 1, 1, 1]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[1, 5]
      pdtheta[k, 1, 6]<-alpha[k, 1, 1, 1]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[1, 6]
      pdtheta[k, 1, 7]<-alpha[k, 1, 1, 1]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[1, 7]
      pdtheta[k, 1, 8]<-alpha[k, 1, 1, 1]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[1, 8]


######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 2, 1]<-alpha[k, 2, 1, 1]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[2, 1]
      pdtheta[k, 2, 2]<-alpha[k, 2, 1, 1]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[2, 2]
      pdtheta[k, 2, 3]<-alpha[k, 2, 1, 1]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[2, 3]
      pdtheta[k, 2, 4]<-alpha[k, 2, 1, 1]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[2, 4]
      pdtheta[k, 2, 5]<-alpha[k, 2, 1, 1]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[2, 5]
      pdtheta[k, 2, 6]<-alpha[k, 2, 1, 1]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[2, 6]
      pdtheta[k, 2, 7]<-alpha[k, 2, 1, 1]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[2, 7]
      pdtheta[k, 2, 8]<-alpha[k, 2, 1, 1]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[2, 8]
           
######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 3, 1]<-alpha[k, 1, 2, 1]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[3, 1]
      pdtheta[k, 3, 2]<-alpha[k, 1, 2, 1]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[3, 2]
      pdtheta[k, 3, 3]<-alpha[k, 1, 2, 1]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[3, 3]
      pdtheta[k, 3, 4]<-alpha[k, 1, 2, 1]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[3, 4]
      pdtheta[k, 3, 5]<-alpha[k, 1, 2, 1]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[3, 5]
      pdtheta[k, 3, 6]<-alpha[k, 1, 2, 1]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[3, 6]
      pdtheta[k, 3, 7]<-alpha[k, 1, 2, 1]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[3, 7]
      pdtheta[k, 3, 8]<-alpha[k, 1, 2, 1]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[3, 8]


######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 4, 1]<-alpha[k, 1, 1, 2]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[4, 1]
      pdtheta[k, 4, 2]<-alpha[k, 1, 1, 2]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[4, 2]
      pdtheta[k, 4, 3]<-alpha[k, 1, 1, 2]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[4, 3]
      pdtheta[k, 4, 4]<-alpha[k, 1, 1, 2]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[4, 4]
      pdtheta[k, 4, 5]<-alpha[k, 1, 1, 2]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[4, 5]
      pdtheta[k, 4, 6]<-alpha[k, 1, 1, 2]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[4, 6]
      pdtheta[k, 4, 7]<-alpha[k, 1, 1, 2]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[4, 7]
      pdtheta[k, 4, 8]<-alpha[k, 1, 1, 2]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[4, 8]

######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 5, 1]<-alpha[k, 2, 2, 1]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[5, 1]
      pdtheta[k, 5, 2]<-alpha[k, 2, 2, 1]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[5, 2]
      pdtheta[k, 5, 3]<-alpha[k, 2, 2, 1]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[5, 3]
      pdtheta[k, 5, 4]<-alpha[k, 2, 2, 1]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[5, 4]
      pdtheta[k, 5, 5]<-alpha[k, 2, 2, 1]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[5, 5]
      pdtheta[k, 5, 6]<-alpha[k, 2, 2, 1]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[5, 6]
      pdtheta[k, 5, 7]<-alpha[k, 2, 2, 1]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[5, 7]
      pdtheta[k, 5, 8]<-alpha[k, 2, 2, 1]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[5, 8]

######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 6, 1]<-alpha[k, 1, 2, 2]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[6, 1]
      pdtheta[k, 6, 2]<-alpha[k, 1, 2, 2]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[6, 2]
      pdtheta[k, 6, 3]<-alpha[k, 1, 2, 2]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[6, 3]
      pdtheta[k, 6, 4]<-alpha[k, 1, 2, 2]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[6, 4]
      pdtheta[k, 6, 5]<-alpha[k, 1, 2, 2]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[6, 5]
      pdtheta[k, 6, 6]<-alpha[k, 1, 2, 2]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[6, 6]
      pdtheta[k, 6, 7]<-alpha[k, 1, 2, 2]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[6, 7]
      pdtheta[k, 6, 8]<-alpha[k, 1, 2, 2]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[6, 8]

######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 7, 1]<-alpha[k, 2, 1, 2]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[7, 1]
      pdtheta[k, 7, 2]<-alpha[k, 2, 1, 2]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[7, 2]
      pdtheta[k, 7, 3]<-alpha[k, 2, 1, 2]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[7, 3]
      pdtheta[k, 7, 4]<-alpha[k, 2, 1, 2]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[7, 4]
      pdtheta[k, 7, 5]<-alpha[k, 2, 1, 2]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[7, 5]
      pdtheta[k, 7, 6]<-alpha[k, 2, 1, 2]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[7, 6]
      pdtheta[k, 7, 7]<-alpha[k, 2, 1, 2]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[7, 7]
      pdtheta[k, 7, 8]<-alpha[k, 2, 1, 2]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[7, 8]

######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)
      pdtheta[k, 8, 1]<-alpha[k, 2, 2, 2]*beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[8, 1]
      pdtheta[k, 8, 2]<-alpha[k, 2, 2, 2]*beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[8, 2]
      pdtheta[k, 8, 3]<-alpha[k, 2, 2, 2]*beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[8, 3]
      pdtheta[k, 8, 4]<-alpha[k, 2, 2, 2]*beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[8, 4]
      pdtheta[k, 8, 5]<-alpha[k, 2, 2, 2]*beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[8, 5]
      pdtheta[k, 8, 6]<-alpha[k, 2, 2, 2]*beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[8, 6]
      pdtheta[k, 8, 7]<-alpha[k, 2, 2, 2]*beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[8, 7]
      pdtheta[k, 8, 8]<-alpha[k, 2, 2, 2]*beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[8, 8]

      b1[k]<-1/(sum(pdtheta[k, 1, ])+sum(pdtheta[k, 2, ])+sum(pdtheta[k, 3, ])+sum(pdtheta[k, 4, ])+
                sum(pdtheta[k, 5, ])+sum(pdtheta[k, 6, ])+sum(pdtheta[k, 7, ])+sum(pdtheta[k, 8, ]))
      for(i in 1:8)
          for(j in 1:8)
                pdtheta[k, i, j]<-b1[k]*pdtheta[k, i, j]

}

# b. transition matrix 

for(i in 1:8)
{
     for(j in 1:8)
     {
          q1<-sum(pdtheta[, i, j])
          q2<-sum(pdtheta[, i, 1])+sum(pdtheta[, i, 2])+sum(pdtheta[, i, 3])+sum(pdtheta[, i, 4])+
              sum(pdtheta[, i, 5])+sum(pdtheta[, i, 6])+sum(pdtheta[, i, 7])+sum(pdtheta[, i, 8])
          A.new[i, j]<-q1/q2
     }
}

# c. non-null distribution 

ptheta1<-array(rep(0, 2*NUM),c(NUM, 2))
ptheta1[, 1]<-alpha[, 1, 1, 1]*beta[, 1, 1, 1]+
              alpha[, 1, 2, 1]*beta[, 1, 2, 1]+
              alpha[, 1, 1, 2]*beta[, 1, 1, 2]+
              alpha[, 1, 2, 2]*beta[, 1, 2, 2]
ptheta1[, 2]<-alpha[, 2, 1, 1]*beta[, 2, 1, 1]+
              alpha[, 2, 2, 1]*beta[, 2, 2, 1]+
              alpha[, 2, 1, 2]*beta[, 2, 1, 2]+
              alpha[, 2, 2, 2]*beta[, 2, 2, 2]
lf.s1<-ptheta1[, 2]/(ptheta1[, 1]+ptheta1[, 2])
q1<-sum(lf.s1)
q2<-sum(lf.s1*x1)
mu1<-q2/q1
q3<-sum(lf.s1*(x1-mu1)*(x1-mu1))
sd1<-sqrt(q3/q1)
f1.new<-c(mu1, sd1)


ptheta2<-array(rep(0, 2*NUM),c(NUM, 2))
ptheta2[, 1]<-alpha[, 1, 1, 1]*beta[, 1, 1, 1]+
              alpha[, 2, 1, 1]*beta[, 2, 1, 1]+
              alpha[, 1, 1, 2]*beta[, 1, 1, 2]+
              alpha[, 2, 1, 2]*beta[, 2, 1, 2]
ptheta2[, 2]<-alpha[, 1, 2, 1]*beta[, 1, 2, 1]+
              alpha[, 2, 2, 1]*beta[, 2, 2, 1]+
              alpha[, 1, 2, 2]*beta[, 1, 2, 2]+
              alpha[, 2, 2, 2]*beta[, 2, 2, 2]
lf.s2<-ptheta2[, 2]/(ptheta2[, 1]+ptheta2[, 2])
q1<-sum(lf.s2)
q2<-sum(lf.s2*x2)
mu2<-q2/q1
q3<-sum(lf.s2*(x2-mu2)*(x2-mu2))
sd2<-sqrt(q3/q1)
f2.new<-c(mu2, sd2)


ptheta3<-array(rep(0, 2*NUM),c(NUM, 2))
ptheta3[, 1]<-alpha[, 1, 1, 1]*beta[, 1, 1, 1]+
              alpha[, 2, 1, 1]*beta[, 2, 1, 1]+
              alpha[, 1, 2, 1]*beta[, 1, 2, 1]+
              alpha[, 2, 2, 1]*beta[, 2, 2, 1]
ptheta3[, 2]<-alpha[, 1, 1, 2]*beta[, 1, 1, 2]+
              alpha[, 2, 1, 2]*beta[, 2, 1, 2]+
              alpha[, 1, 2, 2]*beta[, 1, 2, 2]+
              alpha[, 2, 2, 2]*beta[, 2, 2, 2]
lf.s3<-ptheta3[, 2]/(ptheta3[, 1]+ptheta3[, 2])
q1<-sum(lf.s3)
q2<-sum(lf.s3*x3)
mu3<-q2/q1
q3<-sum(lf.s3*(x3-mu3)*(x3-mu3))
sd3<-sqrt(q3/q1)
f3.new<-c(mu3, sd3)


pii.new<-c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125)

c0<-bwfw.res$c0
Loglikelihood.new<--sum(log(c0))
df1<-abs(Loglikelihood.old-Loglikelihood.new)
diff<-df1

}


# g. return the results of the E-M algorithm

em.var<-list(pii=pii.new, A=A.new, f1=f1.new, f2=f2.new, f3=f3.new, ni=niter)
return (em.var)

}

