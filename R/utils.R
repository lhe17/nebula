## other useful functions

check_reml = function(reml, model)
{
  if((reml!=0) & (reml!=1))
  {stop("reml must be either 0 or 1.")}
  
  if((reml==1) & (model!='NBLMM'))
  {
    reml <- 0
    warning("The value of reml is changed to zero because reml=1 is only supported for NBLMM in the current version.")
  }
  reml
}

check_conv = function(repml, conv)
{
  if(is.nan(repml$loglik))
  {conv = -30}else{
    if(repml$iter==50)
    {
      conv = -20
    }else{
      if(repml$damp==11)
      {conv = -10}else{
        if(repml$damp==12)
        {conv = -40}
      }
    }
  }
  
  # if(sum(is.na(diag(repml$var)))>0 | sum(diag(repml$var)<0,na.rm=T)>0)
  if(RSpectra::eigs_sym(repml$var,1,which='SA')$values[1] < 1e-8)
  {
    conv = -25
  }
  conv
}