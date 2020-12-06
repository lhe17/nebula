// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// #include <unsupported/Eigen/SpecialFunctions>


// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::COLAMDOrdering<int> CO;

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double pmg_ll_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::VectorXd & offset_c,
                    const Eigen::VectorXd & Y_c, const Eigen::VectorXd & fid_c,
                    const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXd & posind_c,
                    const Eigen::VectorXd & posindy_c, const int nind_c, const int k_c,
                    const Eigen::VectorXd & beta_c,const double sigma_c)
{

    double term1 = 0;

    Eigen::VectorXd xtb = offset_c + X_c*beta_c;

    Eigen::VectorXd extb = xtb.array().exp();

    Eigen::VectorXd cumsumxtb(k_c);
    for(int i=0;i<k_c;i++)
    {cumsumxtb(i)= extb.segment(fid_c(i),fid_c(i+1)-fid_c(i)).sum();}

    double exps = exp(sigma_c);
    double alpha = 1/(exps-1);
    double lambda = 1/(sqrt(exps)*(exps-1));

    unsigned int nelem = posind_c.size();
    Eigen::VectorXd cumsumy_ca = cumsumy_c.array() + alpha;
    for(unsigned int i=0;i<nelem;i++)
    {term1 += lgamma(cumsumy_ca(posind_c(i)));}
    term1 -= nelem*lgamma(alpha);

    nelem = posindy_c.size();
    for(unsigned int i=0;i<nelem;i++)
    {term1 += xtb(posindy_c(i))*Y_c(i);}

    term1 += k_c*alpha*log(lambda);

    Eigen::ArrayXd cumsumxtb_t = cumsumxtb.array()+lambda;
    term1 -= cumsumy_ca.dot(cumsumxtb_t.log().matrix());

    term1 = term1*(-1);
    return term1;

}


// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
Eigen::VectorXd pmg_der_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::VectorXd & offset_c,
                              const Eigen::VectorXd & Y_c, const Eigen::VectorXd & fid_c,
                              const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXd & posind_c,
                              const Eigen::VectorXd & posindy_c, const int nb_c, const int nind_c, const int k_c,
                              const Eigen::VectorXd & beta_c,const double sigma_c)
{
    Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+1);

    Eigen::VectorXd xtb = offset_c + X_c*beta_c;
    Eigen::VectorXd extb = xtb.array().exp();
    Eigen::ArrayXd cumsumxtb(k_c);
    Eigen::ArrayXd len(k_c);
    for(int i=0;i<k_c;i++)
    {
        len(i) = fid_c[i+1]-fid_c[i];
        cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
    }

    double exps = exp(sigma_c);
    double alpha = 1/(exps-1);
    double lambda = 1/(sqrt(exps)*(exps-1));

    double alpha_pr = -exps/pow(exps-1,2);
    double lambda_pr = (1-3*exps)/(2*sqrt(exps)*pow(exps-1,2));

    Eigen::ArrayXd ystar = cumsumy_c.array() + alpha;
    Eigen::ArrayXd mustar = cumsumxtb + lambda;
    Eigen::ArrayXd ymustar = ystar/mustar;

    Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
    Eigen::MatrixXd xexb_f(k_c,nb_c);
    for(int i=0;i<k_c;i++)
    {
        xexb_f.row(i) = xexb.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }

    Eigen::VectorXd db = Eigen::VectorXd::Zero(nb_c);
    db -= xexb_f.transpose()*ymustar.matrix();
    unsigned int nelem = posindy_c.size();
    for(unsigned int i=0;i<nelem;i++)
    {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

    double dtau = 0;
    double ldm = log(lambda)*k_c - mustar.log().sum();
    double adlmy = alpha/lambda*k_c - ymustar.sum();
    dtau += alpha_pr*ldm + lambda_pr*adlmy;

    der.segment(0,nb_c) = db;
    der[nb_c] = dtau;

    return (-1)*der;
}

// [[Rcpp::export]]
Eigen::MatrixXd pmg_hes_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::VectorXd & offset_c,
                              const Eigen::VectorXd & Y_c, const Eigen::VectorXd & fid_c,
                              const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXd & posind_c,
                              const Eigen::VectorXd & posindy_c, const int nb_c, const int nind_c, const int k_c,
                              const Eigen::VectorXd & beta_c,const double sigma_c)
{
  
  Eigen::VectorXd xtb = offset_c + X_c*beta_c;
  Eigen::VectorXd extb = xtb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }
  
  double exps = exp(sigma_c);
  double exps_m = exps-1;
  double exps_s = sqrt(exps);
  double alpha = 1/exps_m;
  double lambda = alpha/exps_s;
  double log_lambda = log(lambda);
  
  exps_m = pow(exps_m,2);
  double alpha_pr = -exps/exps_m;
  double lambda_pr = (1-3*exps)/(2*exps_s*exps_m);
  double alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;
  double lambda_dpr = -3*exps/(2*exps_s*exps_m) + (3*exps-1)/(4*exps_s*exps_m) + (3*exps-1)*exps_s/(exps_m*(exps-1));
  
  Eigen::ArrayXd ystar = cumsumy_c.array() + alpha;
  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::ArrayXd ymustar = ystar/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;
  Eigen::ArrayXd imustar = 1/mustar;
  
  Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  for(int i=0;i<k_c;i++)
  {
    xexb_f.row(i) = xexb.block(fid_c[i],0,len(i),nb_c).colwise().sum();
  }
  
  Eigen::MatrixXd hes(nb_c+1,nb_c+1); 
  Eigen::VectorXd tempb;
  for(int i=0;i<nb_c;i++)
  {
    for(int j=i;j<nb_c;j++)
    {
        tempb = xexb.col(i).cwiseProduct(X_c.col(j));
        double tempc = 0;
        for(int k=0;k<k_c;k++)
        {
          tempc += ymumustar(k)*xexb_f(k,i)*xexb_f(k,j) - ymustar(k)*tempb.segment(fid_c(k),len(k)).sum();
        }
        hes(i,j) = tempc;
        if(i!=j)
        {hes(j,i) = tempc;}
    }
  }
  tempb = lambda_pr*ymumustar-alpha_pr*imustar;
  hes.block(0,nb_c,nb_c,1) = xexb_f.transpose()*tempb;
  hes.block(nb_c,0,1,nb_c) = hes.block(0,nb_c,nb_c,1).transpose();
   
  double hes_sigma = alpha_dpr*log_lambda + 2*alpha_pr*lambda_pr/lambda - alpha*pow(lambda_pr/lambda,2) + alpha/lambda*lambda_dpr;
  hes_sigma = hes_sigma*k_c;
  hes_sigma = hes_sigma - (2*alpha_pr*lambda_pr*imustar.sum() + alpha_dpr*mustar.log().sum() + lambda_dpr*ymustar.sum() - lambda_pr*lambda_pr*ymumustar.sum());
  hes(nb_c,nb_c) = hes_sigma;
  
  return hes;
}


// and the inner product returns a scalar
//
// [[Rcpp::export]]
double ptmg_ll_eigen(const Eigen::MatrixXd & X_c, const Eigen::VectorXd & offset_c,
                     const Eigen::VectorXd & Y_c,
                     const Eigen::VectorXd & fam_c, const Eigen::VectorXi & fid_c,
                     const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                     const Eigen::VectorXi & posindy_c, const int nind_c, const int k_c,
                     const Eigen::VectorXd & beta_c, const Eigen::VectorXd & sigma_c)
{
  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);

  double term1 = 0;

  Eigen::VectorXd extb = offset_c + X_c*beta_c;
  unsigned int nelem = posindy_c.size();

  for(unsigned int i=0;i<nelem;i++)
  {term1 += extb(posindy_c(i))*Y_c(i);}

  extb = extb.array().exp();

  Eigen::VectorXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }

  // unsigned int nelem_1 = posind_c.size();
  Eigen::VectorXd cumsumy_ca = cumsumy_c.array() + alpha;

  Eigen::ArrayXd cumsumxtb_t = cumsumxtb.array() + lambda;
  term1 -= cumsumy_ca.dot(cumsumxtb_t.log().matrix());

  term1 += k_c*alpha*log(lambda);
  term1 += nind_c*gamma*log(gamma);

  Eigen::VectorXd sum_elgcp(nind_c);
  cumsumy_ca = (cumsumy_ca.array()/cumsumxtb_t).matrix();


  for(int i=0;i<k_c;i++)
  {
    sum_elgcp.segment(fid_c(i),len(i)) = cumsumy_ca(i)*extb.segment(fid_c(i),len(i));
  }

  term1 += cumsumy_ca.dot(cumsumxtb.matrix());
  // term1 += gamma/phi*sum_elgcp.sum();

  // sum_elgcp = (phi+sum_elgcp.array()).log();
  sum_elgcp = (sum_elgcp.array()+gamma).log();
  term1 -= gamma*sum_elgcp.sum();

  for(unsigned int i=0;i<nelem;i++)
  {term1 -= Y_c(i)*sum_elgcp(posindy_c(i));}

  // term1 = term1*(-1);

  return term1;

}

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Eigen::VectorXd ptmg_der_eigen(const Eigen::MatrixXd & X_c, const Eigen::VectorXd & offset_c,
                               const Eigen::VectorXd & Y_c, const Eigen::VectorXd & fam_c, const Eigen::VectorXi & fid_c,
                               const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                               const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                               const Eigen::VectorXd & beta_c,const Eigen::VectorXd & sigma_c) {

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+2);

  Eigen::VectorXd xtb = offset_c + X_c*beta_c;
  Eigen::VectorXd extb = xtb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  for(int i=0;i<k_c;i++)
  {cumsumxtb(i)= extb.segment(fid_c(i),fid_c(i+1)-fid_c(i)).sum();}

  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);

  double alpha_pr = -exps/pow(exps-1,2);
  double lambda_pr = (1-3*exps)/(2*sqrt(exps)*pow(exps-1,2));

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;
  Eigen::ArrayXd mustar = lambda + cumsumxtb;
  Eigen::ArrayXd ymustar = ystar.array()/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;

  // Eigen::VectorXd gstar = Y_c.array() + gamma;
  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}
  Eigen::VectorXd gstar_phiymustar(nind_c);
  Eigen::ArrayXd len(k_c);
  double slpey = 0;

  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    gstar_phiymustar.segment(fid_c(i),len(i)) = ymustar(i)*extb.segment(fid_c(i),len(i));
  }
  gstar_phiymustar = gamma + gstar_phiymustar.array();
  slpey = gstar_phiymustar.array().log().sum();
  gstar_phiymustar = gstar.array()/gstar_phiymustar.array();

  Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  for(int i=0;i<k_c;i++)
  {
    // len(i) = fid_c[i+1]-fid_c[i];
    xexb_f.row(i) = xexb.block(fid_c[i],0,len(i),nb_c).colwise().sum();
  }
  xexb_f.transposeInPlace();

  Eigen::MatrixXd dbeta_41(k_c,nb_c);
  Eigen::ArrayXd dbeta_42(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start,lent).dot(extb.segment(start,lent));
  }
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd ymumustar_csxtb = ymumustar*cumsumxtb;
  Eigen::ArrayXd ymumustar_dbeta = ymumustar*dbeta_42;
  // Eigen::VectorXd db = - pow(sexp2,2)*(xexb_f*ymumustar_csxtb.matrix()) - dbeta_41*ymustar.matrix() + sexp2*(xexb_f*(ymumustar*dbeta_42).matrix());
  Eigen::VectorXd db = xexb_f*(ymumustar_dbeta - ymumustar_csxtb).matrix() - dbeta_41*ymustar.matrix();
  for(unsigned int i=0;i<nelem;i++)
  {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

  double dtau = 0;
  double ldm = log(lambda)*k_c - mustar.log().sum();
  double adlmy = alpha/lambda*k_c - ymustar.sum();
  dtau += alpha_pr*(cumsumxtb/mustar).sum() - lambda_pr*ymumustar_csxtb.sum();
  dtau -= alpha_pr*(dbeta_42/mustar).sum() - lambda_pr*ymumustar_dbeta.sum();
  dtau += alpha_pr*ldm + lambda_pr*adlmy;

  double dtau2 = log(gamma)*nind_c + nind_c - slpey - gstar_phiymustar.sum();

  der.segment(0,nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c+1] = dtau2;

  return (-1)*der;
}


// [[Rcpp::export]]
Eigen::MatrixXd call_cumsumy(const Eigen::MappedSparseMatrix<double> count, const Eigen::VectorXi & fid, const int k, const int ng)
{

  Eigen::MatrixXd cumsumy(ng,k);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(ng);
  int temp_k = 0;

  for (int i=0; i<count.outerSize(); ++i)
  {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(count,i); it; ++it)
    {
      temp[it.row()] += it.value();
    }
    if(i==(fid[temp_k+1]-1))
    {
      cumsumy.col(temp_k) = temp;
      temp = Eigen::VectorXd::Zero(ng);
      temp_k++;
    }
  }
  return cumsumy;

}

// [[Rcpp::export]]
Rcpp::List call_posindy(const Eigen::MappedSparseMatrix<double> count, const int k, const int nc)
{

  // Eigen::SparseVector<int> cck = count.col(k);
  int nnz = count.col(k).nonZeros();
  Eigen::VectorXi temp(nnz);
  Eigen::VectorXi value(nnz);
  Eigen::VectorXi ytwo(nnz);
  double mct = 0;
  int n_one = 0;
  int n_two = 0;
  int temp_k = 0;
  int nth = 0;
  Eigen::VectorXi n_onetwo(2);


  for (Eigen::MappedSparseMatrix<double>::InnerIterator it(count,k); it; ++it)
  //for (Eigen::SparseVector<int>::InnerIterator it(cck); it; ++it)
  {
    //if(it.value()>0)
    //{
      temp[temp_k] = it.row();
      // temp[temp_k] = it.index();
      int vt = it.value();
      value[temp_k] = vt;
      mct += vt;
      if(vt==1)
      {
        n_one++;
      }else{
        if(vt==2)
        {
          n_two++;
        }else{
          ytwo[nth++] = vt;
        }
      }
      temp_k++;
    //}
  }
  // mct = mct/count.rows();
  mct = mct/nc;
  n_onetwo(0) = n_one;
  n_onetwo(1) = n_two;

  return Rcpp::List::create(Rcpp::Named("posindy") = temp,
                            Rcpp::Named("Y") = value,
                            Rcpp::Named("mct") = mct,
                            Rcpp::Named("n_onetwo") = n_onetwo,
                            Rcpp::Named("ytwo") = ytwo.head(nth));

}


// [[Rcpp::export]]
Rcpp::List center_m(const Eigen::Map<Eigen::MatrixXd> & X_c)
{
  Eigen::MatrixXd cm = X_c.rowwise() - X_c.colwise().mean();
  Eigen::VectorXd sds = (cm.cwiseProduct(cm)).colwise().mean();
  sds = sds.array().sqrt();
  int nc = X_c.cols();
  for(int i=0;i<nc;i++)
  {
    if(sds(i)>0)
    {
      cm.col(i) = cm.col(i)/sds(i);
    }else{
      if(X_c(0,i)!=0)
      {
        if(X_c(0,i)!=1)
        {
          cm.col(i) = X_c.col(i)/X_c(0,i);
        }else{
          cm.col(i) = X_c.col(i);
        }
      }else{
        sds(i) = -1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("pred") = cm,
                            Rcpp::Named("sds") = sds);

}

// [[Rcpp::export]]
Rcpp::List cv_offset(const Eigen::Map<Eigen::VectorXd> & offset_c, int input, const int nc)
{
  Eigen::VectorXd offset(nc);
  double cv = 0;
  double moffset = 1;
  if(input==1)
  {
    offset = offset_c;
    moffset = offset.mean();
  }else{
    offset = Eigen::VectorXd::Constant(nc,1);
  }
  
  if(moffset>0)
  {
    cv = (offset.array()-moffset).square().sum();
    cv = sqrt(cv/nc)/moffset;
  }
  offset = offset.array().log();
  if(input==1)
  {
    moffset = offset.mean();
  }else{
    moffset = 0;
  }
  return Rcpp::List::create(Rcpp::Named("offset") = offset,
                            Rcpp::Named("moffset") = moffset,
                            Rcpp::Named("cv") = cv);

}

// [[Rcpp::export]]
Rcpp::List ptmg_ll_der_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                     const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                     const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                     const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                     const Eigen::VectorXd & beta_c, const Eigen::VectorXd & sigma_c)
{
  double exps = exp(sigma_c(0));
  double exps_m = exps-1;
  double exps_s = sqrt(exps);
  double alpha = 1/exps_m;
  double lambda = alpha/exps_s;
  double gamma = sigma_c(1);
  double log_lambda = log(lambda);
  double log_gamma = log(gamma);

  exps_m = pow(exps_m,2);
  double alpha_pr = -exps/exps_m;
  double lambda_pr = (1-3*exps)/(2*exps_s*exps_m);

  double term1 = 0;

  Eigen::VectorXd extb = offset_c + X_c*beta_c;
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {term1 += extb(posindy_c(i))*Y_c(i);}


  extb = extb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;

  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::VectorXd mustar_log = mustar.log().matrix();

  term1 -= ystar.dot(mustar_log);

  term1 += k_c*alpha*log_lambda;
  term1 += nind_c*gamma*log_gamma;

  Eigen::VectorXd sum_elgcp(nind_c);
  Eigen::ArrayXd ymustar = ystar.array()/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;
  Eigen::VectorXd gstar_phiymustar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  for(unsigned int i=0;i<nelem;i++)
  {gstar_phiymustar(posindy_c(i)) += Y_c(i);}

  for(int i=0;i<k_c;i++)
  {
    sum_elgcp.segment(fid_c(i),len(i)) = ymustar(i)*extb.segment(fid_c(i),len(i));
  }

  term1 += sum_elgcp.sum();

  sum_elgcp = sum_elgcp.array() + gamma;
  gstar_phiymustar = gstar_phiymustar.array()/sum_elgcp.array();

  sum_elgcp = sum_elgcp.array().log();
  double slpey = sum_elgcp.sum();
  term1 -= gamma*slpey;
  for(unsigned int i=0;i<nelem;i++)
  {term1 -= Y_c(i)*sum_elgcp(posindy_c(i));}
  term1 = term1*(-1);

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+2);

  Eigen::MatrixXd xexb;
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  Eigen::MatrixXd dbeta_41(k_c,nb_c);

  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    xexb = X_c.block(start,0,lent,nb_c).array().colwise() * extb.segment(start,lent).array();
    xexb_f.row(i) = xexb.colwise().sum();
    dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb;
  }
  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd dbeta_42(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start,lent).dot(extb.segment(start,lent));
    // dbeta_42[i] = gstar_phiymustar.segment(start,lent).sum();
  }

  Eigen::ArrayXd ymumustar_dbeta_csxtb = ymumustar*(dbeta_42-cumsumxtb);
  Eigen::VectorXd db = xexb_f*ymumustar_dbeta_csxtb.matrix() - dbeta_41*ymustar.matrix();
  for(unsigned int i=0;i<nelem;i++)
  {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}
  double dtau = 0;
  double ldm = log_lambda*k_c - mustar_log.sum();
  double adlmy = exps_s*k_c - ymustar.sum();

  dtau = dtau + alpha_pr*((cumsumxtb-dbeta_42)/mustar).sum();
  dtau = dtau + lambda_pr*ymumustar_dbeta_csxtb.sum();
  dtau += alpha_pr*ldm + lambda_pr*adlmy;

  double dtau2 = log_gamma*nind_c + nind_c - slpey - gstar_phiymustar.sum();

  der.segment(0,nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c+1] = dtau2;

  return Rcpp::List::create(Rcpp::Named("fn") = term1,
                            Rcpp::Named("gr") = (-1)*der);

}


// [[Rcpp::export]]
Rcpp::List ptmg_ll_der_hes_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                             const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                             const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                             const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                             const Eigen::VectorXd & beta_c, const Eigen::VectorXd & sigma_c)
{
  double exps = exp(sigma_c(0));
  double exps_m = exps-1;
  double exps_s = sqrt(exps);
  double alpha = 1/exps_m;
  double lambda = alpha/exps_s;
  double gamma = sigma_c(1);
  double log_lambda = log(lambda);
  double log_gamma = log(gamma);

  exps_m = pow(exps_m,2);
  double alpha_pr = -exps/exps_m;
  double lambda_pr = (1-3*exps)/(2*exps_s*exps_m);
  double alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;
  double lambda_dpr = -3*exps/(2*exps_s*exps_m) + (3*exps-1)/(4*exps_s*exps_m) + (3*exps-1)*exps_s/(exps_m*(exps-1));

  double term1 = 0;
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(nb_c+2,nb_c+2);
  double hes_sigma = alpha_dpr*log_lambda + 2*alpha_pr*lambda_pr/lambda - alpha*pow(lambda_pr/lambda,2) + alpha/lambda*lambda_dpr;
  hes_sigma = hes_sigma*k_c;

  Eigen::VectorXd extb = offset_c + X_c*beta_c;
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {term1 += extb(posindy_c(i))*Y_c(i);}

  extb = extb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;

  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::VectorXd mustar_log = mustar.log().matrix();

  term1 -= ystar.dot(mustar_log);

  term1 += k_c*alpha*log_lambda;
  term1 += nind_c*gamma*log_gamma;

  Eigen::VectorXd sum_elgcp(nind_c);
  Eigen::ArrayXd ymustar = ystar.array()/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;
  Eigen::VectorXd gstar_phiymustar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  for(unsigned int i=0;i<nelem;i++)
  {gstar_phiymustar(posindy_c(i)) += Y_c(i);}

  for(int i=0;i<k_c;i++)
  {
    sum_elgcp.segment(fid_c(i),len(i)) = ymustar(i)*extb.segment(fid_c(i),len(i));
  }

  term1 += sum_elgcp.sum();

  sum_elgcp = sum_elgcp.array() + gamma;
  Eigen::ArrayXd tempa = 1/sum_elgcp.array();
  hessian(nb_c+1,nb_c+1) = nind_c/gamma - 2*tempa.sum();
  gstar_phiymustar = gstar_phiymustar.array()*tempa;

  sum_elgcp = sum_elgcp.array().log();
  double slpey = sum_elgcp.sum();
  term1 -= gamma*slpey;
  for(unsigned int i=0;i<nelem;i++)
  {term1 -= Y_c(i)*sum_elgcp(posindy_c(i));}
  term1 = term1*(-1);

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+2);

  Eigen::MatrixXd xexb = X_c.array().colwise()*extb.array();
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  Eigen::MatrixXd dbeta_41(k_c,nb_c);

  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // xexb = X_c.block(start,0,lent,nb_c).array().colwise() * extb.segment(start,lent).array();
    xexb_f.row(i) = xexb.block(start,0,lent,nb_c).colwise().sum();
    dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
  }
  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd dbeta_42(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start,lent).dot(extb.segment(start,lent));
    // dbeta_42[i] = gstar_phiymustar.segment(start,lent).sum();
  }

  dbeta_42 = dbeta_42-cumsumxtb;

  Eigen::ArrayXd ymumustar_dbeta_csxtb = ymumustar*dbeta_42;
  Eigen::VectorXd db = xexb_f*ymumustar_dbeta_csxtb.matrix() - dbeta_41*ymustar.matrix();
  for(unsigned int i=0;i<nelem;i++)
  {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}
  double dtau = 0;
  double ldm = log_lambda*k_c - mustar_log.sum();
  double adlmy = exps_s*k_c - ymustar.sum();

  Eigen::ArrayXd imustar = 1/mustar;
  Eigen::VectorXd hes_b_l = dbeta_42*imustar;
  double dbim = hes_b_l.sum();
  dtau = dtau - alpha_pr*dbim;
  dtau = dtau + lambda_pr*ymumustar_dbeta_csxtb.sum();
  dtau += alpha_pr*ldm + lambda_pr*adlmy;

  // not necessary
  hes_sigma -= dbim*alpha_dpr - 2*alpha_pr*lambda_pr*hes_b_l.dot(imustar.matrix()) - lambda_dpr*ymumustar_dbeta_csxtb.sum() + 2*lambda_pr*lambda_pr*(ymumustar_dbeta_csxtb*imustar).sum();

  double dtau2 = log_gamma*nind_c + nind_c - slpey - gstar_phiymustar.sum();

  Eigen::VectorXd hes_sigma_beta = -alpha_pr*(xexb_f*imustar.matrix()) + lambda_pr*(xexb_f*ymumustar.matrix());

  hes_sigma = hes_sigma - (2*alpha_pr*lambda_pr*imustar.sum() + alpha_dpr*mustar_log.sum() + lambda_dpr*ymustar.sum() - lambda_pr*lambda_pr*ymumustar.sum());
  Eigen::VectorXd apmmlpymm = alpha_pr*imustar - lambda_pr*ymumustar;


  hes_b_l = gstar_phiymustar.array()*tempa;
  hessian(nb_c+1,nb_c+1) = hessian(nb_c+1,nb_c+1) + hes_b_l.sum();

  Eigen::VectorXd hes_b_l_f(k_c);
  Eigen::MatrixXd hes_b_l_xexb = xexb.array().colwise()*(hes_b_l.array()-tempa);
  Eigen::MatrixXd hes_b_l_xexb_f(k_c,nb_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_xexb_f.row(i) = hes_b_l_xexb.block(start,0,lent,nb_c).colwise().sum();
    // hes_b_l_f(i) = hes_b_l.segment(start,lent).sum();
  }
  Eigen::VectorXd hes_gamma_beta = hes_b_l_xexb_f.transpose()*ymustar.matrix();

  hes_b_l = hes_b_l.array()*extb.array();
  Eigen::VectorXd hes_b2_l = hes_b_l - extb.cwiseProduct(tempa.matrix());
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_f(i) = hes_b2_l.segment(start,lent).sum();
  }
  hes_gamma_beta = hes_gamma_beta - xexb_f*hes_b_l_f.cwiseProduct(ymumustar.matrix());
  hessian.row(nb_c+1).head(nb_c) = hes_gamma_beta;
  hessian.col(nb_c+1).head(nb_c) = hes_gamma_beta;
  hessian(nb_c+1,nb_c) = apmmlpymm.dot(hes_b_l_f);
  hessian(nb_c,nb_c+1) = hessian(nb_c+1,nb_c);

  hes_b2_l = hes_b_l.array()*extb.array();

  Eigen::VectorXd tempb;
  Eigen::VectorXd tempc;
  hes_b_l_xexb = xexb.array().colwise()*hes_b_l.array();
  Eigen::VectorXd hes_b2_l_f(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_xexb_f.row(i) = hes_b_l_xexb.block(start,0,lent,nb_c).colwise().sum();
    hes_b2_l_f(i) = hes_b2_l.segment(start,lent).sum();
  }
  // hes_b_l_xexb_f.transposeInPlace();

  hes_sigma_beta += hes_b_l_xexb_f.transpose()*(apmmlpymm.cwiseProduct(ymustar.matrix())) - xexb_f*(apmmlpymm.cwiseProduct(hes_b2_l_f.cwiseProduct(ymumustar.matrix())));
  hessian.row(nb_c).head(nb_c) = hes_sigma_beta;
  hessian.col(nb_c).head(nb_c) = hes_sigma_beta;

  apmmlpymm = apmmlpymm.cwiseProduct(apmmlpymm);
  hes_sigma += hes_b2_l_f.dot(apmmlpymm);
  hessian(nb_c,nb_c) = hes_sigma;

  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();
  apmmlpymm = ymustar.cwiseProduct(ymustar);
  Eigen::VectorXd ymuymuimu = apmmlpymm.cwiseProduct(imustar.matrix());
  hes_b2_l_f = hes_b2_l_f.array()*ymuymuimu.array()*imustar;

  for(int i=0;i<nb_c;i++)
  {
    for(int j=i;j<nb_c;j++)
    {
      tempb = xexb.col(i).cwiseProduct(X_c.col(j));
      tempc = tempb.cwiseProduct(hes_b_l);

      for(int k=0;k<k_c;k++)
      {
        int start = fid_c[k];
        int lent = len(k);

        hessian(i,j) -= tempb.segment(start,lent).dot(gstar_phiymustar.segment(start,lent))*ymustar(k);
        hessian(i,j) += tempc.segment(start,lent).sum()*apmmlpymm(k);
        //not necessary
        hessian(i,j) += ymumustar_dbeta_csxtb(k)*tempb.segment(start,lent).sum();
      }
      hes_b_l_f = dbeta_41.col(i).cwiseProduct(xexb_f.col(j)) + dbeta_41.col(j).cwiseProduct(xexb_f.col(i));
      hessian(i,j) += hes_b_l_f.dot(ymumustar.matrix());
      hessian(i,j) -= ymuymuimu.dot(hes_b_l_xexb_f.col(i).cwiseProduct(xexb_f.col(j))+hes_b_l_xexb_f.col(j).cwiseProduct(xexb_f.col(i)));
      tempb = xexb_f.col(i).cwiseProduct(xexb_f.col(j));
      hessian(i,j) -= tempb.dot(ymumustar.matrix());
      hessian(i,j) += hes_b2_l_f.dot(tempb);
      //not necessary
      hessian(i,j) += (-2)*tempb.dot((ymumustar_dbeta_csxtb*imustar).matrix());

      if(i!=j)
      {hessian(j,i) = hessian(i,j);}
    }
  }

  der.segment(0,nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c+1] = dtau2;

  return Rcpp::List::create(Rcpp::Named("fn") = term1,
                            Rcpp::Named("gr") = (-1)*der,
                            Rcpp::Named("hes") = (-1)*hessian);

}


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

  //while((step==0)||(likdif>abs(eps*loglik)))
  while((step==0)||(likdif>eps))
  {
    step++;

    double damp = 1;
    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);

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
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array();
    vb2 = vb - vwb.transpose()*temp;

    Eigen::VectorXd dwvw = dw.array()/vw.array();
    // Eigen::VectorXd dbvwbvbdw = vb.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    beta = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    logw = logw + steplogw;

    loglikp = loglik;
    loglik = 0;

    extb = offset_c + X_c*beta;
    w = logw.array().exp();

    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += logw.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
    }
    extb = extb.array().exp();

    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik += alpha*logw.sum() - lambda*w.sum();

    likdif = loglik - loglikp;
    while((likdif<0)||(std::isinf(loglik)))
    {
      damp = damp/2;

      if(damp<0.01)
      {
        likdif = 0;
        loglik = loglikp;
        break;
      }

      for(int i=0;i<k_c;i++)
      {
        if(steplogw(i)<10)
          {damp_w(i)=damp_w(i)/2;}else{
            if((damp<0.5))
            {
              damp_w(i) = 0;
            }else{
              damp_w(i) = 1-0.1/steplogw(i);
            }
          }
      }

      Eigen::VectorXd new_b = damp*stepbeta;
      // Eigen::VectorXd new_w = damp*steplogw;
      Eigen::VectorXd new_w = steplogw.cwiseProduct(damp_w);
      beta = beta - new_b;
      logw = logw - new_w;

      loglik = 0;

      extb = offset_c + X_c*beta;
      w = logw.array().exp();

      for(unsigned int i=0;i<nelem;i++)
      {loglik += extb(posindy_c(i))*Y_c(i);}
      loglik += logw.dot(cumsumy_c);
      for(int i=0;i<k_c;i++)
      {
        extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
      }
      extb = extb.array().exp();

      extbphil = (extb.array()+gamma).log();
      loglik -= gamma*extbphil.sum();
      for(unsigned int i=0;i<nelem;i++)
      {loglik -= Y_c(i)*extbphil(posindy_c(i));}
      loglik += alpha*logw.sum() - lambda*w.sum();

      likdif = loglik - loglikp;

    }

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
    
      four_der = four_der.array()*four_der.array()/(vws*vws);
      sec_ord += 35*four_der.sum()/384;
      
      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      gstar = gstar_extb_phi.array()*(gamma*gamma*gamma-(11*gamma*gamma)*extb.array()+(11*gamma)*extbp-extbp*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()*third_der.array()/(vws*vws);
      sec_ord += 7*temp_der.sum()/96;
    }
     
  }

  // double logdet = vw.array().log().sum();

  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            Rcpp::Named("var") = vb2.inverse(),
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
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

  loglik = loglik - logw.dot(logw)/alpha/2 - k_c/2*log(alpha);


  // double eps = 1e-6;
  double loglikp = 0;
  double likdif = 0;
  int step = 0;
  Eigen::MatrixXd vb(nb_c,nb_c);
  Eigen::VectorXd vw(k_c);
  Eigen::MatrixXd vb2(nb_c,nb_c);
  Eigen::MatrixXd vwb(k_c,nb_c);
  Eigen::VectorXd gstar_extb_phi(nind_c);

  //while((step==0)||(likdif>abs(eps*loglik)))
  while((step==0)||(likdif>eps))
  {
    step++;

    double damp = 1;
    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);

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
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array();
    vb2 = vb - vwb.transpose()*temp;

    Eigen::VectorXd dwvw = dw.array()/vw.array();
    // Eigen::VectorXd dbvwbvbdw = vb.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    beta = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    logw = logw + steplogw;

    loglikp = loglik;
    loglik = 0;

    extb = offset_c + X_c*beta;

    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += logw.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
    }
    extb = extb.array().exp();

    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik = loglik - logw.dot(logw)/alpha/2 - k_c/2*log(alpha);

    likdif = loglik - loglikp;
    while((likdif<0)||(std::isinf(loglik)))
    {
      damp = damp/2;
      if(damp<1e-2)
      {
        likdif = 0;
        loglik = loglikp;
        break;
      }
      
      for(int i=0;i<k_c;i++)
      {
        if(steplogw(i)<10)
          {damp_w(i)=damp_w(i)/2;}else{
            if((damp<0.5))
            {
              damp_w(i) = 0;
            }else{
              damp_w(i) = 1-0.1/steplogw(i);
            }
          }
      }
      
      Eigen::VectorXd new_b = damp*stepbeta;
      // Eigen::VectorXd new_w = damp*steplogw;
      Eigen::VectorXd new_w = steplogw.cwiseProduct(damp_w);
      beta = beta - new_b;
      logw = logw - new_w;

      loglik = 0;

      extb = offset_c + X_c*beta;

      for(unsigned int i=0;i<nelem;i++)
      {loglik += extb(posindy_c(i))*Y_c(i);}
      loglik += logw.dot(cumsumy_c);
      for(int i=0;i<k_c;i++)
      {
        extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
      }
      extb = extb.array().exp();

      extbphil = (extb.array()+gamma).log();
      loglik -= gamma*extbphil.sum();
      for(unsigned int i=0;i<nelem;i++)
      {loglik -= Y_c(i)*extbphil(posindy_c(i));}
      loglik = loglik - logw.dot(logw)/alpha/2 - k_c/2*log(alpha);

      likdif = loglik - loglikp;

    }

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
    vw = vw.array() + 1/alpha;
    Eigen::ArrayXd vws = vw.array()*vw.array();
    
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    
    gstar = gstar_extb_phi.array()*(gamma-extb.array());
    for(int i=0;i<k_c;i++)
    {
      third_der(i) = gstar.segment(fid_c(i),len(i)).sum();
    }
    third_der = gamma*third_der;
    temp_der = third_der.array()*third_der.array()/(vws*vw.array());
    sec_ord += 5*temp_der.sum()/24;
    
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    Eigen::ArrayXd extbp = extb.array()*extb.array();
    gstar = gstar_extb_phi.array()*(gamma*gamma+extbp-4*gamma*extb.array());
    for(int i=0;i<k_c;i++)
    {
      four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
    }
    four_der = gamma*four_der;
    temp_der = four_der.array()/vws;
    sec_ord -= temp_der.sum()/8;
    
    if(ord>2)
    {
      four_der = four_der.array()*four_der.array()/(vws*vws);
      sec_ord += 35*four_der.sum()/384;
      
      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      gstar = gstar_extb_phi.array()*(gamma*gamma*gamma-(11*gamma*gamma)*extb.array()+(11*gamma)*extbp-extbp*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      temp_der = four_der.array()*third_der.array()/(vws*vws);
      sec_ord += 7*temp_der.sum()/96;
    }
    
  }
  

  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            Rcpp::Named("var") = vb2.inverse(),
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
                            Rcpp::Named("second") = sec_ord);

}


// [[Rcpp::export]]
Rcpp::List pml_ll_der_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                            const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                            const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                            const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                            const Eigen::VectorXd & beta_c, const Eigen::VectorXd & logw_c, const Eigen::VectorXd & sigma_c)
{
  double alpha = sigma_c(0);
  double gamma = sigma_c(1);

  double loglik = 0;

  Eigen::VectorXd logw = logw_c;
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

  loglik += alpha*logw.sum() - alpha*w.sum();

  Eigen::VectorXd gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
  Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
  Eigen::VectorXd dw(k_c);
  for(int i=0;i<k_c;i++)
  {
    dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
  }
  dw = cumsumy_c - dw - alpha*w;
  dw = dw.array() + alpha;

  Eigen::VectorXd gr(k_c+nb_c);
  gr.head(nb_c) = db;
  gr.tail(k_c) = dw;

  return Rcpp::List::create(Rcpp::Named("fn") = -loglik,
                            Rcpp::Named("gr") = -gr);

}

