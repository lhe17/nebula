// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// #include <unsupported/Eigen/SpecialFunctions>

// [[Rcpp::export]]
Rcpp::List opt_pml(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                   const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                   const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                   const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                   const Eigen::VectorXd & beta_c,
                   const Eigen::VectorXd & sigma_c, const int reml, const double eps, const int ord)
{
  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);
  
  double loglik = 0;
  
  Eigen::VectorXd logw = Eigen::VectorXd::Constant(k_c,1,0);
  // Eigen::VectorXd logw = logw_c;
  Eigen::VectorXd beta = beta_c;
  
  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}
  
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
  }
  
  Eigen::VectorXd yx = Eigen::VectorXd::Zero(nb_c);
  for(unsigned int i=0;i<nelem;i++)
  {yx += X_c.row(posindy_c(i)).transpose()*Y_c(i);}
  
  Eigen::VectorXd extb = offset_c + X_c*beta;
  
  Eigen::VectorXd w = logw.array().exp();
  
  for(unsigned int i=0;i<nelem;i++)
  {loglik += extb(posindy_c(i))*Y_c(i);}
  
  loglik += logw.dot(cumsumy_c);
  
  for(int i=0;i<k_c;i++)
  {
    extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
  }
  extb = extb.array().exp();
  
  Eigen::ArrayXd extbphil = (extb.array()+gamma).log();
  loglik -= gamma*extbphil.sum();
  for(unsigned int i=0;i<nelem;i++)
  {loglik -= Y_c(i)*extbphil(posindy_c(i));}
  
  loglik += alpha*logw.sum() - lambda*w.sum();
  
  //double eps = 1e-6;
  double loglikp = 0;
  double likdif = 0;
  int step = 0;
  Eigen::MatrixXd vb(nb_c,nb_c);
  Eigen::MatrixXd vb2(nb_c,nb_c);
  Eigen::VectorXd vw(k_c);
  Eigen::MatrixXd vwb(k_c,nb_c);
  Eigen::VectorXd gstar_extb_phi(nind_c);
  int stepd = 0;
  int maxstd = 10;
  int maxstep = 50;
  double convd = 0.01;
  
  //while((step==0)||(likdif>abs(eps*loglik)))
  while((step==0)||((likdif>eps)&&(step<maxstep)))
  {
    step++;
    
    // double damp = 1;
    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);
    Eigen::VectorXd damp = Eigen::VectorXd::Constant(nb_c,1);
    
    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
    Eigen::VectorXd dw(k_c);
    for(int i=0;i<k_c;i++)
    {
      dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    dw = cumsumy_c - dw - lambda*w;
    dw = dw.array() + alpha;
    
    gstar_extb_phi = gstar_extb_phi.array()/(extb.array()+gamma);
    
    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw += lambda*w;
    
    Eigen::MatrixXd xgsetbp = X_c.array().colwise()*gstar_extb_phi.array();
    
    for(int i=0;i<k_c;i++)
    {
      vwb.row(i) = xgsetbp.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }
    vwb = gamma*vwb;
    for(int i=0;i<nb_c;i++)
    {
      for(int j=i;j<nb_c;j++)
      {
        vb(i,j) = X_c.col(i).dot(xgsetbp.col(j));
        if(i!=j)
        {
          vb(j,i) = vb(i,j);
        }
      }
    }
    // vb = X_c.transpose()*vb;
    vb = gamma*vb;
    //Eigen::MatrixXd temp = vwb.array().colwise()/vw.array();
    //vb2 = vb - vwb.transpose()*temp;
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array().sqrt();
    vb2 = vb - temp.transpose()*temp;
    
    Eigen::VectorXd dwvw = dw.array()/vw.array();
    // Eigen::VectorXd dbvwbvbdw = vb.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    //beta = beta + stepbeta;
    Eigen::VectorXd new_b = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    //logw = logw + steplogw;
    Eigen::VectorXd new_w = logw + steplogw;
    
    loglikp = loglik;
    loglik = 0;
    
    extb = offset_c + X_c*new_b;
    w = new_w.array().exp();
    
    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += new_w.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
    }
    extb = extb.array().exp();
    
    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik += alpha*new_w.sum() - lambda*w.sum();
    
    likdif = loglik - loglikp;
    stepd = 0;
    double minstep = 40;
    
    while((likdif<0)||(std::isinf(loglik)))
    {
      stepd++;
      minstep = minstep/2;
      
      if(stepd>maxstd)
      {
        likdif = 0;
        loglik = loglikp;
        double mabsdb = db.cwiseAbs().maxCoeff();
        double mabsdw = dw.cwiseAbs().maxCoeff();
        if((mabsdb>convd) || (mabsdw>convd))
        {stepd++;}
        break;
      }
      
      for(int i=0;i<nb_c;i++)
      {
        if((stepbeta(i)<40)&&(stepbeta(i)>-40))
        {
          damp(i) = damp(i)/2;
          new_b(i) = beta(i) + stepbeta(i)*damp(i);
        }else{
          if(stepbeta(i)>0)
          {
            new_b(i) = beta(i) + minstep;
          }else{
            new_b(i) = beta(i) - minstep;
          }
        }
      }
      
      for(int i=0;i<k_c;i++)
      {
        if((steplogw(i)<40)&&(steplogw(i)>-40))
        {
          damp_w(i)=damp_w(i)/2;
          new_w(i) = logw(i) + steplogw(i)*damp_w(i);
        }else{
          if(steplogw(i)>0)
          {
            new_w(i) = logw(i) + minstep;
          }else{
            new_w(i) = logw(i) - minstep;
          }
        }
      }
      
      loglik = 0;
      
      extb = offset_c + X_c*new_b;
      w = new_w.array().exp();
      
      for(unsigned int i=0;i<nelem;i++)
      {loglik += extb(posindy_c(i))*Y_c(i);}
      loglik += new_w.dot(cumsumy_c);
      for(int i=0;i<k_c;i++)
      {
        extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
      }
      extb = extb.array().exp();
      
      extbphil = (extb.array()+gamma).log();
      loglik -= gamma*extbphil.sum();
      for(unsigned int i=0;i<nelem;i++)
      {loglik -= Y_c(i)*extbphil(posindy_c(i));}
      loglik += alpha*new_w.sum() - lambda*w.sum();
      
      likdif = loglik - loglikp;
      
    }
    beta = new_b;
    logw = new_w;
  }
  
  double logdet;
  if(reml==1)
  {
    logdet = vw.array().abs().log().sum();
    double absdet = vb2.determinant();
    if(absdet<0)
      absdet *= -1;
    logdet += log(absdet);
  }else{
    logdet = vw.array().abs().log().sum();
  }
  
  double sec_ord = 0;
  if(ord>1)
  {
    Eigen::VectorXd third_der = Eigen::VectorXd::Zero(k_c);
    Eigen::VectorXd four_der = Eigen::VectorXd::Zero(k_c);
    Eigen::VectorXd temp_der = Eigen::VectorXd::Zero(k_c);
    
    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::ArrayXd extbg = extb.array()+gamma;
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw += lambda*w;
    Eigen::ArrayXd vws = vw.array()*vw.array();
    
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    
    gstar = gstar_extb_phi.array()*(gamma-extb.array());
    for(int i=0;i<k_c;i++)
    {
      third_der(i) = gstar.segment(fid_c(i),len(i)).sum();
    }
    third_der = gamma*third_der;
    third_der += lambda*w;
    temp_der = third_der.array()*third_der.array()/(vws*vw.array());
    sec_ord += 5*temp_der.sum()/24;
    
    if(ord>2)
    {
      
      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      Eigen::ArrayXd extbp = extb.array()*extb.array();
      gstar = gstar_extb_phi.array()*(gamma*gamma+extbp-(4*gamma)*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()/vws;
      sec_ord -= temp_der.sum()/8;
      
      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      gstar = gstar_extb_phi.array()*(gamma*gamma*gamma-(11*gamma*gamma)*extb.array()+(11*gamma)*extbp-extbp*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()*third_der.array()/(vws*vws);
      sec_ord += 7*temp_der.sum()/48;
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            // Rcpp::Named("var") = vb2.inverse(),
                            Rcpp::Named("var") = vb2,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
                            Rcpp::Named("damp") = stepd,
                            Rcpp::Named("second") = sec_ord);
  
}

// [[Rcpp::export]]
Rcpp::List opt_pml_nbm(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                       const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                       const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                       const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                       const Eigen::VectorXd & beta_c,
                       const Eigen::VectorXd & sigma_c, const int reml, const double eps, const int ord)
{
  double alpha = sigma_c(0);
  double gamma = sigma_c(1);
  
  double loglik = 0;
  
  Eigen::VectorXd logw = Eigen::VectorXd::Constant(k_c,1,0);
  // Eigen::VectorXd logw = logw_c;
  Eigen::VectorXd beta = beta_c;
  
  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}
  
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
  }
  
  Eigen::VectorXd yx = Eigen::VectorXd::Zero(nb_c);
  for(unsigned int i=0;i<nelem;i++)
  {yx += X_c.row(posindy_c(i)).transpose()*Y_c(i);}
  
  Eigen::VectorXd extb = offset_c + X_c*beta;
  
  for(unsigned int i=0;i<nelem;i++)
  {loglik += extb(posindy_c(i))*Y_c(i);}
  
  loglik += logw.dot(cumsumy_c);
  
  for(int i=0;i<k_c;i++)
  {
    extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
  }
  extb = extb.array().exp();
  
  Eigen::ArrayXd extbphil = (extb.array()+gamma).log();
  loglik -= gamma*extbphil.sum();
  for(unsigned int i=0;i<nelem;i++)
  {loglik -= Y_c(i)*extbphil(posindy_c(i));}
  
  loglik = loglik - logw.dot(logw)/alpha/2 - k_c/2.0*log(alpha);
  
  
  // double eps = 1e-6;
  double loglikp = 0;
  double likdif = 0;
  int step = 0;
  Eigen::MatrixXd vb(nb_c,nb_c);
  Eigen::VectorXd vw(k_c);
  Eigen::MatrixXd vb2(nb_c,nb_c);
  Eigen::MatrixXd vwb(k_c,nb_c);
  Eigen::VectorXd gstar_extb_phi(nind_c);
  int stepd = 0;
  int maxstd = 10;
  int maxstep = 50;
  double convd = 0.01;
  
  //while((step==0)||(likdif>abs(eps*loglik)))
  while((step==0)||((likdif>eps)&&(step<maxstep)))
  {
    step++;
    
    // double damp = 1;
    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);
    Eigen::VectorXd damp = Eigen::VectorXd::Constant(nb_c,1);
    
    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
    Eigen::VectorXd dw(k_c);
    for(int i=0;i<k_c;i++)
    {
      dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    dw = cumsumy_c-dw-logw/alpha;
    
    gstar_extb_phi = gstar_extb_phi.array()/(extb.array()+gamma);
    
    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw = vw.array() + 1/alpha;
    
    Eigen::MatrixXd xgsetbp = X_c.array().colwise()*gstar_extb_phi.array();
    
    for(int i=0;i<k_c;i++)
    {
      vwb.row(i) = xgsetbp.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }
    vwb = gamma*vwb;
    for(int i=0;i<nb_c;i++)
    {
      for(int j=i;j<nb_c;j++)
      {
        vb(i,j) = X_c.col(i).dot(xgsetbp.col(j));
        if(i!=j)
        {
          vb(j,i) = vb(i,j);
        }
      }
    }
    // vb = X_c.transpose()*vb;
    vb = gamma*vb;
    // Eigen::MatrixXd temp = vwb.array().colwise()/vw.array();
    // vb2 = vb - vwb.transpose()*temp;
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array().sqrt();
    vb2 = vb - temp.transpose()*temp;
    
    Eigen::VectorXd dwvw = dw.array()/vw.array();
    // Eigen::VectorXd dbvwbvbdw = vb.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd new_b = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    Eigen::VectorXd new_w = logw + steplogw;
    
    loglikp = loglik;
    loglik = 0;
    
    extb = offset_c + X_c*new_b;
    
    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += new_w.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
    }
    extb = extb.array().exp();
    
    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik = loglik - new_w.dot(new_w)/alpha/2 - k_c/2.0*log(alpha);
    
    likdif = loglik - loglikp;
    stepd = 0;
    double minstep = 40;
    
    while((likdif<0)||(std::isinf(loglik)))
    {
      stepd++;
      minstep = minstep/2;
      
      if(stepd>maxstd)
      {
        likdif = 0;
        loglik = loglikp;
        double mabsdb = db.cwiseAbs().maxCoeff();
        double mabsdw = dw.cwiseAbs().maxCoeff();
        if((mabsdb>convd) || (mabsdw>convd))
        {stepd++;}
        break;
      }
      
      for(int i=0;i<nb_c;i++)
      {
        if((stepbeta(i)<40)&&(stepbeta(i)>-40))
        {
          damp(i)=damp(i)/2;
          new_b(i) = beta(i) + stepbeta(i)*damp(i);
        }else{
          if(stepbeta(i)>0)
          {
            new_b(i) = beta(i) + minstep;
          }else{
            new_b(i) = beta(i) - minstep;
          }
        }
      }
      
      for(int i=0;i<k_c;i++)
      {
        if((steplogw(i)<40)&&(steplogw(i)>-40))
        {
          damp_w(i)=damp_w(i)/2;
          new_w(i) = logw(i) + steplogw(i)*damp_w(i);
        }else{
          if(steplogw(i)>0)
          {
            new_w(i) = logw(i) + minstep;
          }else{
            new_w(i) = logw(i) - minstep;
          }
        }
      }
      
      loglik = 0;
      
      extb = offset_c + X_c*new_b;
      
      for(unsigned int i=0;i<nelem;i++)
      {loglik += extb(posindy_c(i))*Y_c(i);}
      loglik += new_w.dot(cumsumy_c);
      for(int i=0;i<k_c;i++)
      {
        extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
      }
      extb = extb.array().exp();
      
      extbphil = (extb.array()+gamma).log();
      loglik -= gamma*extbphil.sum();
      for(unsigned int i=0;i<nelem;i++)
      {loglik -= Y_c(i)*extbphil(posindy_c(i));}
      loglik = loglik - new_w.dot(new_w)/alpha/2 - k_c/2.0*log(alpha);
      
      likdif = loglik - loglikp;
    }
    
    beta = new_b;
    logw = new_w;
  }
  
  double logdet;
  if(reml==1)
  {
    logdet = vw.array().abs().log().sum();
    double absdet = vb2.determinant();
    if(absdet<0)
      absdet *= -1;
    logdet += log(absdet);
  }else{
    logdet = vw.array().abs().log().sum();
  }
  // double logdet = vw.array().log().sum();
  
  double sec_ord = 0;
  
  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            //Rcpp::Named("var") = vb2.inverse(),
                            Rcpp::Named("var") = vb2,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
                            Rcpp::Named("damp") = stepd,
                            Rcpp::Named("second") = sec_ord);
  
}

