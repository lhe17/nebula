ptmg_ll = function(para,X,offset,Y,n_one,ytwo,fam,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = para[(nb+1):(nb+2)]

  exps = exp(sigma[1])
  alpha = 1/(exps-1)
  gamma = sigma[2]

  fn = ptmg_ll_eigen(X,offset,Y,fam-1,fid-1,as.double(cumsumy),posind-1,posindy,nind,k,beta,sigma)
  fn = fn + (sum(Lgamma(cumsumy[posind]+alpha))-length(posind)*Lgamma(alpha))
  fn = fn + (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  return(-fn)
}

ptmg_der = function(para,X,offset,Y,n_one,ytwo,fam,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = para[(nb+1):(nb+2)]

  exps = exp(sigma[1])
  alpha = 1/(exps-1)
  gamma = sigma[2]
  alpha_pr = -exps/((exps-1)^2)

  temp = ptmg_der_eigen(X,offset,Y,fam-1,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  temp[nb+1] = temp[nb+1] - alpha_pr*(sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha))
  temp[nb+2] = temp[nb+2] - (sum(Digamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Digamma(gamma) + sum(n_one)/gamma + n_one[2]/(gamma+1))
  return(temp)
}

ptmg_ll_der = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = para[(nb+1):(nb+2)]

  exps = exp(sigma[1])
  alpha = 1/(exps-1)
  gamma = sigma[2]
  alpha_pr = -exps/((exps-1)^2)

  temp = ptmg_ll_der_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  temp$fn = temp$fn - (sum(Lgamma(cumsumy[posind]+alpha))-length(posind)*Lgamma(alpha))
  temp$fn = temp$fn - (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  gr = temp$gr
  gr[nb+1] = gr[nb+1] - alpha_pr*(sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha))
  gr[nb+2] = gr[nb+2] - (sum(Digamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Digamma(gamma) + sum(n_one)/gamma + n_one[2]/(gamma+1))
  return(list("objective"=temp$fn,"gradient"=gr))
}

ptmg_ll_der2 = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = exp(para[(nb+1):(nb+2)])

  exps = exp(sigma[1])
  alpha = 1/(exps-1)
  gamma = sigma[2]
  alpha_pr = -exps/((exps-1)^2)

  temp = ptmg_ll_der_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  temp$fn = temp$fn - (sum(Lgamma(cumsumy[posind]+alpha))-length(posind)*Lgamma(alpha))
  temp$fn = temp$fn - (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  gr = temp$gr
  gr[nb+1] = gr[nb+1] - alpha_pr*(sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha))
  gr[nb+2] = gr[nb+2] - (sum(Digamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Digamma(gamma) + sum(n_one)/gamma + n_one[2]/(gamma+1))
  gr[nb+1] = gr[nb+1]*sigma[1]
  gr[nb+2] = gr[nb+2]*sigma[2]
  # return(list("objective"=temp$fn,"gradient"=gr))
  objective = temp$fn
  attr(objective, "gradient") = gr
  objective
}


ptmg_ll_der_hes2 = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = exp(para[(nb+1):(nb+2)])

  exps = exp(sigma[1])
  exps_m = (exps-1)^2
  alpha = 1/(exps-1)
  gamma = sigma[2]
  alpha_pr = -exps/((exps-1)^2)
  alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;

  temp = ptmg_ll_der_hes_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  temp$fn = temp$fn - (sum(Lgamma(cumsumy[posind]+alpha))-length(posind)*Lgamma(alpha))
  temp$fn = temp$fn - (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  gr = temp$gr
  grdig = sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha)
  gr[nb+1] = gr[nb+1] - alpha_pr*grdig
  gr[nb+2] = gr[nb+2] - (sum(Digamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Digamma(gamma) + sum(n_one)/gamma + n_one[2]/(gamma+1))
  gr[nb+1] = gr[nb+1]*sigma[1]
  gr[nb+2] = gr[nb+2]*sigma[2]
  hes = temp$hes
  hes[nb+2,nb+2] = hes[nb+2,nb+2] - (sum(Trigamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Trigamma(gamma) - sum(n_one)/gamma/gamma - n_one[2]/(gamma+1)/(gamma+1))
  hes[nb+1,nb+1] = hes[nb+1,nb+1] - alpha_dpr*grdig - alpha_pr*alpha_pr*(sum(Trigamma(cumsumy[posind]+alpha))-length(posind)*Trigamma(alpha))

  hes[nb+2,nb+2] = hes[nb+2,nb+2]*sigma[2]*sigma[2] + gr[nb+2]
  hes[nb+1,nb+1] = hes[nb+1,nb+1]*sigma[1]*sigma[1] + gr[nb+1]
  hes[nb+1,nb+2] = hes[nb+2,nb+1] = hes[nb+1,nb+2]*sigma[2]*sigma[1]
  hes[nb+1,1:nb] = hes[1:nb,nb+1] = hes[nb+1,1:nb]*sigma[1]
  hes[nb+2,1:nb] = hes[1:nb,nb+2] = hes[nb+2,1:nb]*sigma[2]
  objective = temp$fn
  attr(objective, "gradient") = gr
  attr(objective, "hessian") = hes
  objective
  # return(list("objective"=-objective,"gradient"=-gr, "hessian"=-hes))
}

ptmg_ll_der_hes3 = function(para,X,offset,Y,n_one,ytwo,fid,cumsumy,posind,posindy,nb,nind,k)
{
  beta = para[1:nb]
  sigma = exp(para[(nb+1):(nb+2)])

  exps = exp(sigma[1])
  exps_m = (exps-1)^2
  alpha = 1/(exps-1)
  gamma = sigma[2]
  alpha_pr = -exps/((exps-1)^2)
  alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;

  temp = ptmg_ll_der_hes_eigen(X,offset,Y,fid-1,as.double(cumsumy),posind-1,posindy,nb,nind,k,beta,sigma)
  temp$fn = temp$fn - (sum(Lgamma(cumsumy[posind]+alpha))-length(posind)*Lgamma(alpha))
  temp$fn = temp$fn - (sum(Lgamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Lgamma(gamma) + sum(n_one)*log(gamma)+ n_one[2]*log(gamma+1))
  gr = temp$gr
  grdig = sum(Digamma(cumsumy[posind]+alpha))-length(posind)*Digamma(alpha)
  gr[nb+1] = gr[nb+1] - alpha_pr*grdig
  gr[nb+2] = gr[nb+2] - (sum(Digamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Digamma(gamma) + sum(n_one)/gamma + n_one[2]/(gamma+1))
  gr[nb+1] = gr[nb+1]*sigma[1]
  gr[nb+2] = gr[nb+2]*sigma[2]
  hes = temp$hes
  hes[nb+2,nb+2] = hes[nb+2,nb+2] - (sum(Trigamma(ytwo+gamma))-(length(posindy)-sum(n_one))*Trigamma(gamma) - sum(n_one)/gamma/gamma - n_one[2]/(gamma+1)/(gamma+1))
  hes[nb+1,nb+1] = hes[nb+1,nb+1] - alpha_dpr*grdig - alpha_pr*alpha_pr*(sum(Trigamma(cumsumy[posind]+alpha))-length(posind)*Trigamma(alpha))

  hes[nb+2,nb+2] = hes[nb+2,nb+2]*sigma[2]*sigma[2] + gr[nb+2]
  hes[nb+1,nb+1] = hes[nb+1,nb+1]*sigma[1]*sigma[1] + gr[nb+1]
  hes[nb+1,nb+2] = hes[nb+2,nb+1] = hes[nb+1,nb+2]*sigma[2]*sigma[1]
  hes[nb+1,1:nb] = hes[1:nb,nb+1] = hes[nb+1,1:nb]*sigma[1]
  hes[nb+2,1:nb] = hes[1:nb,nb+2] = hes[nb+2,1:nb]*sigma[2]
  objective = temp$fn
  #attr(objective, "gradient") = gr
  #attr(objective, "hessian") = hes
  #objective
  return(list("value"=objective,"gradient"=gr, "hessian"=hes))
}
