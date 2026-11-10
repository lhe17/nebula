// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// #include <unsupported/Eigen/SpecialFunctions>

// [[Rcpp::export]]
Eigen::MatrixXd call_cumsumy(const Eigen::Map<Eigen::SparseMatrix<double>> count, const Eigen::VectorXi & fid, const int k, const int ng)
{
  
  Eigen::MatrixXd cumsumy(ng,k);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(ng);
  int temp_k = 0;
  
  for (int i=0; i<count.outerSize(); ++i)
  {
    for (Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator it(count,i); it; ++it)
    {
      if(it.value()>0)
      {
        temp[it.row()] += it.value();
      }
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
Rcpp::List call_posindy(const Eigen::Map<Eigen::SparseMatrix<double>> count, const int k, const int nc)
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
  
  
  for (Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator it(count,k); it; ++it)
    //for (Eigen::SparseVector<int>::InnerIterator it(cck); it; ++it)
  {
    if(it.value()>0)
    {
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
    }
  }
  // mct = mct/count.rows();
  mct = mct/nc;
  n_onetwo(0) = n_one;
  n_onetwo(1) = n_two;
  
  return Rcpp::List::create(Rcpp::Named("posindy") = temp.head(temp_k),
                            Rcpp::Named("Y") = value.head(temp_k),
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
  int nr = X_c.rows();
  for(int i=0;i<nc;i++)
  {
    if(sds(i)>0)
    {
      cm.col(i) = cm.col(i)/sds(i);
    }else{
      if(X_c(0,i)!=0)
      {
        /*
         if(X_c(0,i)!=1)
         {
         cm.col(i) = X_c.col(i)/X_c(0,i);
         }else{
         cm.col(i) = X_c.col(i);
         }
         */
        cm.col(i) = Eigen::VectorXd::Constant(nr,1);
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
  double mexpoffset = moffset;
  
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
                            Rcpp::Named("mexpoffset") = mexpoffset,
                            Rcpp::Named("cv") = cv);
  
}

// [[Rcpp::export]]
double get_cv(const Eigen::Map<Eigen::VectorXd> offset_c, const Eigen::Map<Eigen::MatrixXd> X_c, 
              const Eigen::VectorXd & beta_c, const Eigen::VectorXi & cell_ind, const int ncell, const int nc)
{
  Eigen::ArrayXd extb = offset_c;
  for(int i=0; i<ncell; i++)
  {
    int ind = cell_ind(i)-1;
    extb = extb + X_c.col(ind).array()*beta_c(ind);
  }
  extb = extb.array().exp();
  double cv = 0;
  double moffset = extb.mean();
  if(moffset>0)
  {
    cv = (extb-moffset).square().sum();
    cv = cv/nc/(moffset*moffset);
  }
  return cv;
}

// [[Rcpp::export]]
Eigen::VectorXd get_cell(const Eigen::Map<Eigen::MatrixXd> X_c, const Eigen::VectorXi & fid_c, const int nb_c, const int k_c)
{
  Eigen::VectorXd iscell = Eigen::VectorXd::Zero(nb_c);
  for(int i=0; i<nb_c; i++)
  {
    for(int j=0;j<k_c;j++)
    {
      int flag = 0;
      int end = fid_c[j+1];
      int start = fid_c(j);
      double temp = X_c(start,i);
      for(int k=start;k<end;k++)
      {
        if(X_c(k,i)!=temp)
        {
          flag = 1;
          break;
        }
      }
      if(flag==1)
      {
        iscell(i) = 1;
        break;
      }
    }
  }
  return iscell;
}

