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

check_conv = function(repml, conv, nb, vare, min, max, cutoff = 1e-8)
{
  if(conv==1)
  {
    if(vare[1]==max[1] | vare[2]==min[2])
    {
      conv = -60
    }else{
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
    }
  }
  
  if(nb>1)
  {
    if(min(eigen(repml$var)$values) < cutoff)
    {conv = -25}
  }
  
  conv
}