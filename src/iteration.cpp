#include <RcppArmadillo.h>   
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List iter (int max_iter,
                 float tol, 
                 arma::mat VD,
                 arma::mat VD2,
                 arma::mat V,
                 arma::mat B,
                 double alpha,
                 double beta,
                 double nu,
                 double kappa) {
  
  //initiate number of iteration and improvement for while loop
  int    noi = 0;                
  double improvement = 5000;
  
  //initiate objects for svd and matrix manipulation
  arma::vec obj;
  arma::mat Z;
  arma::mat U;
  arma::vec d;
  arma::mat VV;
  
  while (noi <= max_iter && improvement > tol){
    Z = VD2 * V.t() * B;        
    svd(U, d, VV, Z);
    int n = U.n_rows;
    int p = VV.n_cols;
    arma::mat UU = U.submat(0, 0, n-1, p-1);       //subset U for matrix conformity
    arma::mat A = UU * VV.t();                     //update A: X'XB = UDV'

    //Proximal Gradient Descent to Update B
    arma::mat grad = VD2 * ( V.t() * (A - B) ) - beta * B;
    arma::mat B_temp = B + nu * grad;
    
    //l1 soft-threshold
    int ncol = B_temp.n_cols;
    int nrow = B_temp.n_rows;
    for (int i = 0; i< nrow; i++ ){       //correct B based on kappa value
      for (int j = 0; j < ncol; j++){
        if ( B_temp(i,j) > kappa ){         
          B(i,j) = B_temp(i,j) - kappa;
        }else if ( B_temp(i,j) <= -kappa ){
          B(i,j) = B_temp(i,j) + kappa;
        }
      }
    }

    //compute residual
    arma::mat R = VD.t() - ( VD.t() * B ) * A.t();
    
    //compute objective function, append to the objective vector
    int sz = obj.size();
    obj.resize(sz + 1);
    double R2 = 0.5 * accu( pow(R,2) ) + alpha * accu( abs(B) ) + 0.5 * beta * accu( pow(B,2) );
    obj(sz) = R2;

    //break if not improving
    if (noi > 0){
      improvement = ( obj(noi - 1) - obj(noi) ) / obj(noi);
    }
    
    //next iter
    noi += 1;
  }

  //return a list of svd$d, mat B, and mat obj
  Rcpp::List ret;
  ret["d"] = d;
  ret["B"] = B;

  return ret;
}

