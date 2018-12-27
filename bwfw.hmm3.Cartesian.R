bwfw.hmm3.Cartesian<-function(x1, x2, x3, pii, A, f0, f1, f2, f3)
{

###################################################################


## Initialize

NUM<-length(x1)


## Densities

f0x1<-dnorm(x1, f0[1], f0[2])
f0x2<-dnorm(x2, f0[1], f0[2])
f0x3<-dnorm(x3, f0[1], f0[2])
f1x1<-dnorm(x1, f1[1], f1[2])
f1x2<-dnorm(x2, f2[1], f2[2])
f1x3<-dnorm(x3, f3[1], f3[2])


## the backward-forward procedure

# a. the backward variables
# --rescaled 

alpha<-array(rep(0, NUM*8), c(NUM, 2, 2, 2))
# scaling variable c_0
c0<-rep(0, NUM)

######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)

alpha[1, 1, 1, 1]<-pii[1]*f0x1[1]*f0x2[1]*f0x3[1]
alpha[1, 2, 1, 1]<-pii[2]*f1x1[1]*f0x2[1]*f0x3[1]
alpha[1, 1, 2, 1]<-pii[3]*f0x1[1]*f1x2[1]*f0x3[1]
alpha[1, 2, 2, 1]<-pii[5]*f1x1[1]*f1x2[1]*f0x3[1]
alpha[1, 1, 1, 2]<-pii[4]*f0x1[1]*f0x2[1]*f1x3[1]
alpha[1, 2, 1, 2]<-pii[7]*f1x1[1]*f0x2[1]*f1x3[1]
alpha[1, 1, 2, 2]<-pii[6]*f0x1[1]*f1x2[1]*f1x3[1]
alpha[1, 2, 2, 2]<-pii[8]*f1x1[1]*f1x2[1]*f1x3[1] 


# rescaling alpha
c0[1]<-1/(alpha[1, 1, 1, 1]+alpha[1, 2, 1, 1]+alpha[1, 1, 2, 1]+alpha[1, 2, 2, 1]+
          alpha[1, 1, 1, 2]+alpha[1, 2, 1, 2]+alpha[1, 1, 2, 2]+alpha[1, 2, 2, 2])
for(i in 1:2)
    for(j in 1:2)
        for(k in 1:2)
              alpha[1, i, j, k]<-c0[1]*alpha[1, i, j, k]
                                                                  
for (k in 1:(NUM-1))
{ 
      alpha[k+1, 1, 1, 1]<-(alpha[k, 1, 1, 1]*A[1, 1]+alpha[k, 2, 1, 1]*A[2, 1]+
                            alpha[k, 1, 2, 1]*A[3, 1]+alpha[k, 2, 2, 1]*A[5, 1]+
                            alpha[k, 1, 1, 2]*A[4, 1]+alpha[k, 2, 1, 2]*A[7, 1]+
                            alpha[k, 1, 2, 2]*A[6, 1]+alpha[k, 2, 2, 2]*A[8, 1])*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]
     
      alpha[k+1, 2, 1, 1]<-(alpha[k, 1, 1, 1]*A[1, 2]+alpha[k, 2, 1, 1]*A[2, 2]+
                            alpha[k, 1, 2, 1]*A[3, 2]+alpha[k, 2, 2, 1]*A[5, 2]+
                            alpha[k, 1, 1, 2]*A[4, 2]+alpha[k, 2, 1, 2]*A[7, 2]+
                            alpha[k, 1, 2, 2]*A[6, 2]+alpha[k, 2, 2, 2]*A[8, 2])*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]

      alpha[k+1, 1, 2, 1]<-(alpha[k, 1, 1, 1]*A[1, 3]+alpha[k, 2, 1, 1]*A[2, 3]+
                            alpha[k, 1, 2, 1]*A[3, 3]+alpha[k, 2, 2, 1]*A[5, 3]+
                            alpha[k, 1, 1, 2]*A[4, 3]+alpha[k, 2, 1, 2]*A[7, 3]+
                            alpha[k, 1, 2, 2]*A[6, 3]+alpha[k, 2, 2, 2]*A[8, 3])*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]

      alpha[k+1, 2, 2, 1]<-(alpha[k, 1, 1, 1]*A[1, 5]+alpha[k, 2, 1, 1]*A[2, 5]+
                            alpha[k, 1, 2, 1]*A[3, 5]+alpha[k, 2, 2, 1]*A[5, 5]+
                            alpha[k, 1, 1, 2]*A[4, 5]+alpha[k, 2, 1, 2]*A[7, 5]+
                            alpha[k, 1, 2, 2]*A[6, 5]+alpha[k, 2, 2, 2]*A[8, 5])*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]

      alpha[k+1, 1, 1, 2]<-(alpha[k, 1, 1, 1]*A[1, 4]+alpha[k, 2, 1, 1]*A[2, 4]+
                            alpha[k, 1, 2, 1]*A[3, 4]+alpha[k, 2, 2, 1]*A[5, 4]+
                            alpha[k, 1, 1, 2]*A[4, 4]+alpha[k, 2, 1, 2]*A[7, 4]+
                            alpha[k, 1, 2, 2]*A[6, 4]+alpha[k, 2, 2, 2]*A[8, 4])*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]

      alpha[k+1, 2, 1, 2]<-(alpha[k, 1, 1, 1]*A[1, 7]+alpha[k, 2, 1, 1]*A[2, 7]+
                            alpha[k, 1, 2, 1]*A[3, 7]+alpha[k, 2, 2, 1]*A[5, 7]+
                            alpha[k, 1, 1, 2]*A[4, 7]+alpha[k, 2, 1, 2]*A[7, 7]+
                            alpha[k, 1, 2, 2]*A[6, 7]+alpha[k, 2, 2, 2]*A[8, 7])*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]

      alpha[k+1, 1, 2, 2]<-(alpha[k, 1, 1, 1]*A[1, 6]+alpha[k, 2, 1, 1]*A[2, 6]+
                            alpha[k, 1, 2, 1]*A[3, 6]+alpha[k, 2, 2, 1]*A[5, 6]+
                            alpha[k, 1, 1, 2]*A[4, 6]+alpha[k, 2, 1, 2]*A[7, 6]+
                            alpha[k, 1, 2, 2]*A[6, 6]+alpha[k, 2, 2, 2]*A[8, 6])*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]

      alpha[k+1, 2, 2, 2]<-(alpha[k, 1, 1, 1]*A[1, 8]+alpha[k, 2, 1, 1]*A[2, 8]+
                            alpha[k, 1, 2, 1]*A[3, 8]+alpha[k, 2, 2, 1]*A[5, 8]+
                            alpha[k, 1, 1, 2]*A[4, 8]+alpha[k, 2, 1, 2]*A[7, 8]+
                            alpha[k, 1, 2, 2]*A[6, 8]+alpha[k, 2, 2, 2]*A[8, 8])*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]

      # rescaling alpha
      c0[k+1]<-1/(alpha[k+1, 1, 1, 1]+alpha[k+1, 2, 1, 1]+alpha[k+1, 1, 2, 1]+alpha[k+1, 2, 2, 1]+
                  alpha[k+1, 1, 1, 2]+alpha[k+1, 2, 1, 2]+alpha[k+1, 1, 2, 2]+alpha[k+1, 2, 2, 2])
      for(i in 1:2)
            for(j in 1:2)
                  for(p in 1:2)
                          alpha[k+1, i, j, p]<-c0[k+1]*alpha[k+1, i, j, p]
}

# b. the forward variables
# --rescaled

beta<-array(rep(0, NUM*8), c(NUM, 2, 2, 2))

beta[NUM, 1, 1, 1]<-1/8
beta[NUM, 2, 1, 1]<-1/8
beta[NUM, 1, 2, 1]<-1/8
beta[NUM, 2, 2, 1]<-1/8
beta[NUM, 1, 1, 2]<-1/8
beta[NUM, 2, 1, 2]<-1/8
beta[NUM, 1, 2, 2]<-1/8
beta[NUM, 2, 2, 2]<-1/8
b0<-rep(0, NUM)


######1:(0,0,0); 2:(1,0,0); 3:(0,1,0); 4:(0,0,1) 
######5:(1,1,0); 6:(0,1,1); 7:(1,0,1); 8:(1,1,1)

for (k in (NUM-1):1)
{ 
  beta[k, 1, 1, 1]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[1, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[1, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[1, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[1, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[1, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[1, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[1, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[1, 8])

  beta[k, 2, 1, 1]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[2, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[2, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[2, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[2, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[2, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[2, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[2, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[2, 8])

  beta[k, 1, 2, 1]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[3, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[3, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[3, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[3, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[3, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[3, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[3, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[3, 8])

  beta[k, 1, 1, 2]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[4, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[4, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[4, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[4, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[4, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[4, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[4, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[4, 8])

  beta[k, 2, 2, 1]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[5, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[5, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[5, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[5, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[5, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[5, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[5, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[5, 8])

  beta[k, 1, 2, 2]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[6, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[6, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[6, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[6, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[6, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[6, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[6, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[6, 8])

  beta[k, 2, 1, 2]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[7, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[7, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[7, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[7, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[7, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[7, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[7, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[7, 8])

  beta[k, 2, 2, 2]<-(beta[k+1, 1, 1, 1]*f0x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[8, 1]+
                     beta[k+1, 2, 1, 1]*f1x1[k+1]*f0x2[k+1]*f0x3[k+1]*A[8, 2]+
                     beta[k+1, 1, 2, 1]*f0x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[8, 3]+
                     beta[k+1, 1, 1, 2]*f0x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[8, 4]+
                     beta[k+1, 2, 2, 1]*f1x1[k+1]*f1x2[k+1]*f0x3[k+1]*A[8, 5]+
                     beta[k+1, 1, 2, 2]*f0x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[8, 6]+
                     beta[k+1, 2, 1, 2]*f1x1[k+1]*f0x2[k+1]*f1x3[k+1]*A[8, 7]+
                     beta[k+1, 2, 2, 2]*f1x1[k+1]*f1x2[k+1]*f1x3[k+1]*A[8, 8])


  # rescaling beta
  # using the same scaling factors as alpha 
  b0[k]<-1/(beta[k, 1, 1, 1]+beta[k, 2, 1, 1]+beta[k, 1, 2, 1]+beta[k, 2, 2, 1]+
            beta[k, 1, 1, 2]+beta[k, 2, 1, 2]+beta[k, 1, 2, 2]+beta[k, 2, 2, 2])
  for(i in 1:2)
       for(j in 1:2)
          for(p in 1:2)
               beta[k, i, j, p]<-b0[k]*beta[k, i, j, p]
}

# c. lfdr variables
# --original
# --the same formulae hold for the rescaled alpha and beta

lfdr<-rep(0, NUM)


for (k in 1:NUM)
{ 
  q1<-(alpha[k, 1, 1, 1]*beta[k, 1, 1, 1]+alpha[k, 2, 1, 1]*beta[k, 2, 1, 1]+
       alpha[k, 1, 2, 1]*beta[k, 1, 2, 1]+alpha[k, 1, 1, 2]*beta[k, 1, 1, 2])
  q2<-(alpha[k, 2, 2, 1]*beta[k, 2, 2, 1]+alpha[k, 1, 2, 2]*beta[k, 1, 2, 2]+
       alpha[k, 2, 1, 2]*beta[k, 2, 1, 2]+alpha[k, 2, 2, 2]*beta[k, 2, 2, 2])
  lfdr[k]<-q1/(q1+q2)
}



# f. return the results of the bwfw proc.

bwfw.var<-list(bw=alpha, fw=beta, repLIS=lfdr, c0=c0)


return(bwfw.var)
  
}
