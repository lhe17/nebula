pmg_ll = function(para,X,offset,Y,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = para[nb+1]

  temp = pmg_ll_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nind,k,beta,sigma)
  return(temp)
}


pmg_der = function(para,X,offset,Y,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = para[nb+1]

  exps = exp(sigma[1])
  alpha = 1/(exps-1)
  alpha_pr = -exps/((exps-1)^2)

  temp = pmg_der_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  temp[nb+1] = temp[nb+1] - alpha_pr*(sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha))
  return(temp)
}

pmg_hes = function(para,X,offset,Y,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = para[nb+1]
  
  exps = exp(sigma)
  alpha = 1/(exps-1)
  exps_m = (exps-1)^2
  alpha_pr = -exps/exps_m
  alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;
  
  grdig = sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha)
  
  hes = pmg_hes_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  hes[nb+1,nb+1] = hes[nb+1,nb+1] + alpha_dpr*grdig + alpha_pr*alpha_pr*(sum(Trigamma(cumsumy[posind]+alpha))-length(posind)*Trigamma(alpha))
  return(hes)
}


