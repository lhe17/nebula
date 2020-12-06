## penalized likelihood

pml_ll_der = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k,sigmas)
{
  betas = para[1:nb]
  logw = para[(nb+1):(nb+k)]
  temp = pml_ll_der_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,betas,logw,sigmas)
  return(list("objective"=temp$fn,"gradient"=temp$gr))
}

pml_ll_der2 = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k,sigmas)
{
  betas = para[1:nb]
  logw = para[(nb+1):(nb+k)]
  temp = pml_ll_der_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,betas,logw,sigmas)
  # return(list("objective"=temp$fn,"gradient"=temp$gr))
  objective = temp$fn
  attr(objective, "gradient") = temp$gr
  objective
}


pql_ll = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k,betas,reml,eps,ord)
{
  exps = exp(para[1])
  alpha = 1/(exps-1)
  gamma = para[2]
  lambda = 1/(sqrt(exps)*(exps-1))

  repml = opt_pml(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,betas,para,reml,eps,ord)
  if(repml$second < (-1))
  {stop()}

  loglik = repml$loglik + nind*gamma*log(gamma) + k*alpha*log(lambda) - k*Lgamma(alpha)
  loglik = loglik + (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  loglik = loglik - 0.5*repml$logdet + log(1+repml$second)
  -loglik
}

pql_gamma_ll = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k,betas,reml,eps,gamma,ord)
{
  exps = exp(para[1])
  alpha = 1/(exps-1)
  lambda = 1/(sqrt(exps)*(exps-1))

  repml = opt_pml(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,betas,c(para[1],gamma),reml,eps,ord)
  if(repml$second < (-1))
  {stop()}

  loglik = repml$loglik + nind*gamma*log(gamma) + k*alpha*log(lambda) - k*Lgamma(alpha)
  loglik = loglik + (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  loglik = loglik - 0.5*repml$logdet + log(1+repml$second)
  -loglik
}

pql_nbm_ll = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k,betas,reml,eps,ord)
{
  alpha = para[1]
  gamma = para[2]

  repml = opt_pml_nbm(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,betas,para,reml,eps,1)
  loglik = repml$loglik + nind*gamma*log(gamma) + sum(Lgamma(Y+gamma)) -length(Y)*Lgamma(gamma)
  loglik = loglik - 0.5*repml$logdet
  -loglik

}

pql_nbm_gamma_ll = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k,betas,reml,eps,gamma,ord)
{
  alpha = para[1]
  
  repml = opt_pml_nbm(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,betas,c(para[1],gamma),reml,eps,1)
  loglik = repml$loglik + nind*gamma*log(gamma) + sum(Lgamma(Y+gamma)) -length(Y)*Lgamma(gamma)
  loglik = loglik - 0.5*repml$logdet
  -loglik
}
