##################################################################################
###MM algorithm for Plackett-Luce.

Bayesian.MM=function(init,tiny,alpha,beta,playersmatrix,n){
  iterationmatrix=matrix(0,1000,n)
  t=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  while(sum(abs(x0-x1)>tiny)>0){
    x0=x1
    for(i in 1:n){
      l=c()
      b1=c()
      for(y in 1:dim(playersmatrix)[1]){
        b=c()
        a=c()
        a1=c()
        players_at_yth_race=playersmatrix[y,]
        n0=sum(playersmatrix[y,]!=0)
        for(r in 1:(n0-1)){
          a[r]=sum(players_at_yth_race[r:n0]==i)
          #a[r]=1 if driver i receives rank no better than r in y-th race.
          b[r]=sum(exp(x0[players_at_yth_race[r:n0]]))
          a1[r]=a[r]/b[r]
        }
        b1[y]=sum(a1)
        l[y]=sum(players_at_yth_race[1:(n0-1)]==i)
        #l[y]=1 if driver i is not in the last ranking position in y-th race.
      }
      x1[i]=log(alpha-1+sum(l))-log(beta+sum(b1))
      #sum(l) is the total number that driver i is not in the last ranking position.
    }
    t=t+1
    iterationmatrix[t,]=x1
  }
  iterationmatrix=iterationmatrix[1:t,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

##################################################################################
###AccMM algorithm for Plackett-Luce.
Acc.Bayesian.MM=function(init,tiny,alpha,beta,playersmatrix,n){
  iterationmatrix=matrix(0,1000,n)
  t=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  while(sum(abs(x0-x1)>tiny)>0){
    x0=x1
    for(i in 1:n){
      l=c()
      b1=c()
      for(y in 1:dim(playersmatrix)[1]){
        b=c()
        a=c()
        a1=c()
        players_at_yth_race=playersmatrix[y,]
        n0=sum(playersmatrix[y,]!=0)
        for(r in 1:(n0-1)){
          a[r]=sum(players_at_yth_race[r:n0]==i)
          #a[r]=1 if driver i receives rank no better than r in y-th race.
          b[r]=sum(exp(x0[players_at_yth_race[r:n0]]))
          a1[r]=a[r]/b[r]
        }
        b1[y]=sum(a1)
        l[y]=sum(players_at_yth_race[1:(n0-1)]==i)
        #l[y]=1 if driver i is not in the last ranking position in y-th race.
      }
      x1[i]=log(alpha-1+sum(l))-log(beta+sum(b1))
      #sum(l) is the total number that driver i is not in the last ranking position.
    }
    c=log(n*(alpha-1)/beta)-log(sum(exp(x1)))
    x1=x1+c
    t=t+1
    iterationmatrix[t,]=x1
  }
  iterationmatrix=iterationmatrix[1:t,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}
##################################################################################
#value settings for experiments MM and AccMM:

init=rep(0,n)
tiny=10^(-4)

ImatrixMM0=Bayesian.MM(init,tiny,1,0,playersmatrix,n)$Iterationmatrix
ImatrixMM1.1=Bayesian.MM(init,tiny,1.01,0.01,playersmatrix,n)$Iterationmatrix
ImatrixMM1=Bayesian.MM(init,tiny,1.1,0.1,playersmatrix,n)$Iterationmatrix
ImatrixMM2=Bayesian.MM(init,tiny,2,1,playersmatrix,n)$Iterationmatrix
ImatrixMM3=Bayesian.MM(init,tiny,11,10,playersmatrix,n)$Iterationmatrix

dim(ImatrixMM0)
dim(ImatrixMM1.1)
dim(ImatrixMM1)
dim(ImatrixMM2)
dim(ImatrixMM3)

AccImatrixMM1.1=Acc.Bayesian.MM(init,tiny,1.01,0.01,playersmatrix,n)$Iterationmatrix
AccImatrixMM1=Acc.Bayesian.MM(init,tiny,1.1,0.1,playersmatrix,n)$Iterationmatrix
AccImatrixMM2=Acc.Bayesian.MM(init,tiny,2,1,playersmatrix,n)$Iterationmatrix
AccImatrixMM3=Acc.Bayesian.MM(init,tiny,11,10,playersmatrix,n)$Iterationmatrix

dim(AccImatrixMM1.1)
dim(AccImatrixMM1)
dim(AccImatrixMM2)
dim(AccImatrixMM3)
