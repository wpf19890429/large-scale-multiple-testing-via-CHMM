bwfw.hmm.Cartesian<-function(x1, x2, pii, A, f0, f1, f2)
{

###################################################################



## Initialize

NUM<-length(x1)

## Densities

f0x1<-dnorm(x1, f0[1], f0[2])
f0x2<-dnorm(x2, f0[1], f0[2])
f1x1<-dnorm(x1, f1[1], f1[2])
f1x2<-dnorm(x2, f2[1], f2[2])


## the backward-forward procedure

# a. the backward variables
# --rescaled 

alpha<-array(rep(0, NUM*4), c(NUM, 2, 2))
# scaling variable c_0
c0<-rep(0, NUM)

alpha[1, 1, 1]<-pii[1]*f0x1[1]*f0x2[1]
alpha[1, 2, 1]<-pii[2]*f1x1[1]*f0x2[1]
alpha[1, 1, 2]<-pii[3]*f0x1[1]*f1x2[1]
alpha[1, 2, 2]<-pii[4]*f1x1[1]*f1x2[1]

# rescaling alpha
c0[1]<-1/(alpha[1, 1, 1]+alpha[1, 2, 1]+alpha[1, 1, 2]+alpha[1, 2, 2])
for(i in 1:2)
    for(j in 1:2)
          alpha[1, i, j]<-c0[1]*alpha[1, i, j]

for (k in 1:(NUM-1))
{ 
  alpha[k+1, 1, 1]<-(alpha[k, 1, 1]*A[1, 1]+alpha[k, 2, 1]*A[2, 1]+
                     alpha[k, 1, 2]*A[3, 1]+alpha[k, 2, 2]*A[4, 1])*f0x1[k+1]*f0x2[k+1]
  alpha[k+1, 2, 1]<-(alpha[k, 1, 1]*A[1, 1]+alpha[k, 2, 1]*A[2, 1]+
                     alpha[k, 1, 2]*A[3, 1]+alpha[k, 2, 2]*A[4, 1])*f1x1[k+1]*f0x2[k+1]
  alpha[k+1, 1, 2]<-(alpha[k, 1, 1]*A[1, 1]+alpha[k, 2, 1]*A[2, 1]+
                     alpha[k, 1, 2]*A[3, 1]+alpha[k, 2, 2]*A[4, 1])*f0x1[k+1]*f1x2[k+1]
  alpha[k+1, 2, 2]<-(alpha[k, 1, 1]*A[1, 1]+alpha[k, 2, 1]*A[2, 1]+
                     alpha[k, 1, 2]*A[3, 1]+alpha[k, 2, 2]*A[4, 1])*f1x1[k+1]*f1x2[k+1]
  # rescaling alpha
  c0[k+1]<-1/(alpha[k+1, 1, 1]+alpha[k+1, 2, 1]+alpha[k+1, 1, 2]+alpha[k+1, 2, 2])
  for(i in 1:2)
    for(j in 1:2)
          alpha[k+1, i, j]<-c0[k+1]*alpha[k+1, i, j]
}

# b. the forward variables
# --rescaled

beta<-array(rep(0, NUM*4), c(NUM, 2, 2))

beta[NUM, 1, 1]<-1/4
beta[NUM, 2, 1]<-1/4
beta[NUM, 1, 2]<-1/4
beta[NUM, 2, 2]<-1/4
b0<-rep(0, NUM)


for (k in (NUM-1):1)
{ 
  beta[k, 1, 1]<-(beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[1, 1]+beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[1, 2]+
                  beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[1, 3]+beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[1, 4])
  beta[k, 2, 1]<-(beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[2, 1]+beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[2, 2]+
                  beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[2, 3]+beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[2, 4])
  beta[k, 1, 2]<-(beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[3, 1]+beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[3, 2]+
                  beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[3, 3]+beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[3, 4])
  beta[k, 2, 2]<-(beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[4, 1]+beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[4, 2]+
                  beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[4, 3]+beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[4, 4])
  # rescaling beta
  # using the same scaling factors as alpha 
  b0[k]<-1/(beta[k, 1, 1]+beta[k, 2, 1]+beta[k, 1, 2]+beta[k, 2, 2])
  for(i in 1:2)
    for(j in 1:2)
          beta[k, i, j]<-b0[k]*beta[k, i, j]
}

# c. lfdr variables
# --original
# --the same formulae hold for the rescaled alpha and beta

lfdr<-rep(0, NUM)


for (k in 1:NUM)
{ 
  q1<-alpha[k, 1, 1]*beta[k, 1, 1]+alpha[k, 2, 1]*beta[k, 2, 1]+alpha[k, 1, 2]*beta[k, 1, 2]
  q2<-alpha[k, 2, 2]*beta[k, 2, 2]
  lfdr[k]<-q1/(q1+q2)
}





lfdr1<-rep(0, NUM)
for (k in 1:NUM)
{ 
  q1<-alpha[k, 1, 1]*beta[k, 1, 1]+alpha[k, 1, 2]*beta[k, 1, 2]
  q2<-alpha[k, 2, 2]*beta[k, 2, 2]+alpha[k, 2, 1]*beta[k, 2, 1]
  lfdr1[k]<-q1/(q1+q2)
}


lfdr2<-rep(0, NUM)
for (k in 1:NUM)
{ 
  q1<-alpha[k, 2, 2]*beta[k, 2, 2]+alpha[k, 1, 2]*beta[k, 1, 2]
  q2<-alpha[k, 1, 1]*beta[k, 1, 1]+alpha[k, 2, 1]*beta[k, 2, 1]
  lfdr2[k]<-q1/(q1+q2)
}


# d. probabilities of hidden states
# -- and transition variables
# -- both are rescaled


pdtheta<-array(rep(0, 16*(NUM-1)),c((NUM-1), 2, 2, 2, 2))


b1<-rep(0, NUM-1)
for(k in 1:(NUM-1))
{
      pdtheta[k, 1, 1, 1, 1]<-alpha[k, 1, 1]*beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[1, 1]
      pdtheta[k, 2, 1, 1, 1]<-alpha[k, 2, 1]*beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[2, 1] 
      pdtheta[k, 1, 2, 1, 1]<-alpha[k, 1, 2]*beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[3, 1]
      pdtheta[k, 2, 2, 1, 1]<-alpha[k, 2, 2]*beta[k+1, 1, 1]*f0x1[k+1]*f0x2[k+1]*A[4, 1]

      pdtheta[k, 1, 1, 2, 1]<-alpha[k, 1, 1]*beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[1, 2]
      pdtheta[k, 2, 1, 2, 1]<-alpha[k, 2, 1]*beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[2, 2] 
      pdtheta[k, 1, 2, 2, 1]<-alpha[k, 1, 2]*beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[3, 2]
      pdtheta[k, 2, 2, 2, 1]<-alpha[k, 2, 2]*beta[k+1, 2, 1]*f1x1[k+1]*f0x2[k+1]*A[4, 2]

      pdtheta[k, 1, 1, 1, 2]<-alpha[k, 1, 1]*beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[1, 3]
      pdtheta[k, 2, 1, 1, 2]<-alpha[k, 2, 1]*beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[2, 3] 
      pdtheta[k, 1, 2, 1, 2]<-alpha[k, 1, 2]*beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[3, 3]
      pdtheta[k, 2, 2, 1, 2]<-alpha[k, 2, 2]*beta[k+1, 1, 2]*f0x1[k+1]*f1x2[k+1]*A[4, 3]

      pdtheta[k, 1, 1, 2, 2]<-alpha[k, 1, 1]*beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[1, 4]
      pdtheta[k, 2, 1, 2, 2]<-alpha[k, 2, 1]*beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[2, 4] 
      pdtheta[k, 1, 2, 2, 2]<-alpha[k, 1, 2]*beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[3, 4]
      pdtheta[k, 2, 2, 2, 2]<-alpha[k, 2, 2]*beta[k+1, 2, 2]*f1x1[k+1]*f1x2[k+1]*A[4, 4]
 
      b1[k]<-1/(pdtheta[k, 1, 1, 1, 1]+pdtheta[k, 1, 2, 1, 1]+pdtheta[k, 2, 1, 1, 1]+pdtheta[k, 2, 2, 1, 1]+
                pdtheta[k, 1, 1, 2, 1]+pdtheta[k, 1, 2, 2, 1]+pdtheta[k, 2, 1, 2, 1]+pdtheta[k, 2, 2, 2, 1]+
                pdtheta[k, 1, 1, 1, 2]+pdtheta[k, 1, 2, 1, 2]+pdtheta[k, 2, 1, 1, 2]+pdtheta[k, 2, 2, 1, 2]+
                pdtheta[k, 1, 1, 2, 2]+pdtheta[k, 1, 2, 2, 2]+pdtheta[k, 2, 1, 2, 2]+pdtheta[k, 2, 2, 2, 2])


      for(i in 1:2)
          for(j in 1:2)
               for(p in 1:2)
                   for(q in 1:2)
                             pdtheta[k, i, j, p, q]<-b1[k]*pdtheta[k, i, j, p, q]


    
}





# f. return the results of the bwfw proc.

bwfw.var<-list(bw=alpha, fw=beta, lsi=lfdr, lfdr1=lfdr1, lfdr2=lfdr2, pdthe=pdtheta, c0=c0)


return(bwfw.var)
  
}
