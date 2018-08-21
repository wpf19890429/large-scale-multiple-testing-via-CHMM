em.hmm.Cartesian<-function(x1, x2, maxiter=200)
{

#####################################################################################

NUM<-length(x1)
# precision tolerance level
ptol<-1e-4
niter<-0

### initializing model parameters

pii.new<-c(1, 0, 0, 0)

A.new<-matrix(c(c(0.25, 0.25, 0.25, 0.25), 
                c(0.25, 0.25, 0.25, 0.25),
                c(0.25, 0.25, 0.25, 0.25),
                c(0.25, 0.25, 0.25, 0.25)), 4, 4, byrow=TRUE)
f0<-c(0, 1)
f1.new<-c(2, 1)
f2.new<-c(2, 1)

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
Loglikelihood.old<-Loglikelihood.new
## updating the weights and probabilities of hidden states

bwfw.res<-bwfw.hmm.Cartesian(x1, x2, pii.old, A.old, f0, f1.old, f2.old)

# the backward-forward variable

alpha<-bwfw.res$bw
beta<-bwfw.res$fw

# the hidden states probabilities
gamma<-bwfw.res$pr
# the transition variables
dgamma<-bwfw.res$ts

## updating the parameter estimates

# a. initial state distribution



# transform variable of theta

pdtheta<-bwfw.res$pdthe


# b. transition matrix of Main effect

for (i in 1:2)
{
  for (j in 1:2)
  { 
      for(p in 1:2)
      {
            for(q in 1:2)
            {
                   q1<-sum(pdtheta[, i, j, p, q])
                   q2<-sum(pdtheta[, i, j, 1, 1])+sum(pdtheta[, i, j, 2, 1])+sum(pdtheta[, i, j, 1, 2])+sum(pdtheta[, i, j, 2, 2])     
                   if(j==1)
                          index_j<-j-1
                   else
                          index_j<-j
                   if(q==1)
                          index_q<-q-1
                   else
                          index_q<-q 
                   A.new[i+index_j, p+index_q]<-q1/q2 
            }
      } 
  }
}
# c. non-null distribution 

ptheta1<-array(rep(0, 2*NUM),c(NUM, 2))
ptheta1[, 1]<-alpha[, 1, 1]*beta[, 1, 1]+alpha[, 1, 2]*beta[, 1, 2]
ptheta1[, 2]<-alpha[, 2, 2]*beta[, 2, 2]+alpha[, 2, 1]*beta[, 2, 1]
lf.s1<-ptheta1[, 2]/(ptheta1[, 1]+ptheta1[, 2])
q1<-sum(lf.s1)
q2<-sum(lf.s1*x1)
mu1<-q2/q1
q3<-sum(lf.s1*(x1-mu1)*(x1-mu1))
sd1<-sqrt(q3/q1)
f1.new<-c(mu1, sd1)



ptheta2<-array(rep(0, 2*NUM),c(NUM, 2))
ptheta2[, 1]<-alpha[, 1, 1]*beta[, 1, 1]+alpha[, 2, 1]*beta[, 2, 1]
ptheta2[, 2]<-alpha[, 2, 2]*beta[, 2, 2]+alpha[, 1, 2]*beta[, 1, 2]
lf.s2<-ptheta2[, 2]/(ptheta2[, 1]+ptheta2[, 2])
q1<-sum(lf.s2)
q2<-sum(lf.s2*x2)
mu2<-q2/q1
q3<-sum(lf.s2*(x2-mu2)*(x2-mu2))
sd2<-sqrt(q3/q1)
f2.new<-c(mu2, sd2)

pii.new<-c(1, 0, 0, 0)

c0<-bwfw.res$c0
Loglikelihood.new<--sum(log(c0))
df1<-abs(Loglikelihood.old-Loglikelihood.new)
diff<-df1


}


# g. return the results of the E-M algorithm

em.var<-list(pii=pii.new, A=A.new, f1=f1.new, f2=f2.new, ni=niter)
return (em.var)

}


