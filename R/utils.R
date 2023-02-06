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

