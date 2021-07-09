#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace Eigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::LLT<MatrixXd> llt;

inline MatrixXd AtA(MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint());
}
inline MatrixXd mappedAtA(Map<MatrixXd>& A) {
  int n(A.cols());
  return MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint());
}

inline arma::umat get_locations(arma::sp_mat& B)
{
  
  // Make const iterator
  arma::sp_mat::const_iterator start = B.begin();
  arma::sp_mat::const_iterator end   = B.end();
  
  // Calculate number of points
  int n = std::distance(start, end);
  
  // Kill process if no values are found (very sparse matrix)
  if (n <= 0) { Rcpp::stop("No values found!"); }
  
  // Build a location storage matrix
  arma::umat locs(2, n);
  
  // Create a vector to store each row information in. (Row, Col)
  arma::uvec temp(2);
  
  // Start collecting locations
  arma::sp_mat::const_iterator it = start;
  for(int i = 0; i < n; ++i)
  {
    temp(0) = it.row();
    temp(1) = it.col();
    locs.col(i) = temp;
    ++it; // increment
  }
  
  return locs;
}

/*** R
lambda_trans = function(LambdaS, LambdaSt, len, ind, j) {
  LambdaS@x <- LambdaSt@x <- len
  LambdaS@x[which(ind == j)] <- LambdaSt@x[which(ind == j)] <- 1
  diagonal <- diag(LambdaS)
  diag(LambdaS) <- diag(LambdaSt) <- 0
  Dj <- LambdaS + LambdaSt
  diag(Dj) <- diagonal
  # Dj <- as.matrix(Dj)
  return(Dj)
}
*/

inline Map<MatrixXd> reshape (MatrixXd &b, const size_t n, const size_t m) {
  return Map<MatrixXd>(b.data(), n, m);
}

inline VectorXd dnormLog(VectorXd x, VectorXd means, double sds) {
  int n = x.size();
  VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x[i], means[i], sds, TRUE);
  }
  return res;
}


// [[Rcpp::export]]
List cAICeigen(
    Map<VectorXd> y,
    Map<VectorXd> y_hat,
    Map<MatrixXd> X,
    double& sigma2,
    Map<MatrixXd> G,
    Map<MatrixXd> Z,
    int& q,
    int& m
  ) {
  const int n = y.size();
  
  VectorXd e = y - y_hat;
  
  MatrixXd D = G / sigma2;
  
  llt Lambda(D);
  // SpMat Lambda = llt.matrixL();
  MatrixXd Lambdat = Lambda.matrixU();
  
  arma::mat armaG = arma::mat(G.data(), G.rows(), G.cols(), false, false);
  arma::sp_mat spLambdat = arma::sp_mat(chol(armaG / sigma2));
  
  SparseMatrix<double> splambdat = Lambdat.sparseView();
  S4 lambdat(wrap(splambdat));
  S4 lambda(wrap(splambdat.adjoint()));
  
  
  //----------------------------------------------
  // V0inv
  MatrixXd V0 = MatrixXd::Identity(n, n) + Z * AtA(Lambdat) * Z.adjoint();
  MatrixXd LLt = Lambdat * mappedAtA(Z) * Lambdat.adjoint() + MatrixXd::Identity(Lambdat.cols(), Lambdat.cols());
  llt cholLLt(LLt);
  
  MatrixXd Linv = cholLLt.matrixL().solve(MatrixXd::Identity(LLt.cols(), LLt.cols()));
  MatrixXd LLZt = Linv * (Lambdat * Z.adjoint());
  MatrixXd V0inv = MatrixXd::Identity(n, n) - AtA(LLZt);
  
  //----------------------------------------------
  // A
  MatrixXd XtV0invX = X.adjoint() * (V0inv * X);
  llt Rx(XtV0invX);
  
  MatrixXd XRxinv = X * Rx.matrixU().solve(MatrixXd::Identity(X.cols(), X.cols()));
  MatrixXd XRxinvtV0inv = XRxinv.adjoint() * V0inv;
  MatrixXd A = V0inv - AtA(XRxinvtV0inv);
  
  //----------------------------------------------
  // W
  int nc = q + 1;
  IntegerVector ncvec1(seq_len(nc));
  IntegerVector ncvec2(seq_len(nc));
  IntegerVector rowIndices(sum(ncvec2));
  int indp = 0;
  for (int i = 0; i < nc; ++i) {
    int ncp = ncvec2[i];
    std::fill(rowIndices.begin() + indp, rowIndices.begin() + indp + ncp, ncvec1[i]);
    indp += ncp;
  }
  
  IntegerVector colIndices(sum(ncvec2));
  indp = 0;
  for (int i = 0; i < nc; ++i) {
    int ncp = ncvec2[i];
    for (int j = 0; j < ncp; ++j) {
      colIndices[indp + j] = ncvec2[j];
    }
    indp += ncp;
  }
  
  IntegerVector LindTemplate;
  IntegerVector col_ind(colIndices.length());
  for (int i; i < colIndices.length(); ++i) {
    col_ind[i] = Rf_choose(colIndices[i], 2);
  }
  LindTemplate = rowIndices + nc * (colIndices - 1) - col_ind;
  
  IntegerVector Lind;
  Lind = rep(LindTemplate, m);
  
  IntegerVector ind = Lind;
  arma::umat gl = get_locations(spLambdat);
  IntegerVector len = rep(0, gl.n_cols);
  
  int n_LindT = LindTemplate.length();
  std::vector<MatrixXd> Wlist(n_LindT);
  std::vector<MatrixXd> eWelist(n_LindT);
  
  
  Environment Gol = Environment::global_env();
  Function lambda_trans = Gol["lambda_trans"];
  for (int j = 0; j < n_LindT; ++j) {
    S4 Djr;
    Djr = lambda_trans(lambda, lambdat, len, ind, j + 1);
    SparseMatrix<double> sDj(as<Map<SparseMatrix<double>> >(Djr));
    MatrixXd Dj = MatrixXd(sDj);
    
    Wlist[j] = Z * Dj * Z.adjoint();
    
    eWelist[j] = e.adjoint() * Wlist[j] * e;
  }
  
  //-------------------------------------------------------
  // B
  double tye = y.adjoint() * e;
  
  std::vector<MatrixXd> WAlist(n_LindT);
  std::vector<MatrixXd> tv_WAlist(n_LindT);
  std::vector<MatrixXd> v_WAlist(n_LindT);

  for (int i = 0; i < n_LindT; ++i) {
    WAlist[i] = Wlist[i] * A;
    MatrixXd wat = WAlist[i].adjoint();
    int n_wa = pow(WAlist[i].cols(), 2);
    v_WAlist[i] = reshape(WAlist[i], n_wa, 1);
    tv_WAlist[i] = reshape(wat, n_wa, 1);

    // Map<RowVectorXf> rvxf(WAlist[i].data(), WAlist[i].size());
    // tv_WAlist[i] = WAlist[i].reshaped().transpose();
  }

  MatrixXd B(n_LindT, n_LindT);
  MatrixXd C(n_LindT, n);

  for (int j = 0; j < n_LindT; ++j) {
    // MatrixXd Wj = Wlist[j];
    // MatrixXd eWje = eWelist[j];

    C.row(j) = e.adjoint() * (Wlist[j] * A) - eWelist[j] * (e.adjoint()) / (2 * tye);

    for (int k = 0; k < n_LindT; ++k) {
      double WkAWjA = (tv_WAlist[j].adjoint() * v_WAlist[k]).sum();
      
      VectorXd tmp = eWelist[j] * eWelist[k] / (2 * tye);
      VectorXd tmp1 = e.adjoint() * Wlist[k] * (A * (Wlist[j] * e));
      
      B(j, k) = - tye * WkAWjA / (2 * (n - X.cols())) - tmp[0] + tmp1[0];
      B(k, j) = B(j, k);
    }
  }
  
  //-------------------------------------------------------
  // df
  MatrixXd Lambday = B.inverse() * C;
  
  double df = n - A.diagonal().sum();
  for (int j = 0; j < n_LindT; ++j) {
    VectorXd tmp = Lambday.row(j) * (A * (Wlist[j] * e));
    df = df + tmp[0];
  }
  df = df + 1;
  
  //================================================================
  // cAIC
  double bc = df;
  
  VectorXd cllv = dnormLog(y, y_hat, sqrt(sigma2));
  double cll = cllv.sum();
  
  double cAIC = -2 * cll + 2 * bc;
  
  return List::create(Named("cAIC") = cAIC,
                      Named("bc") = bc,
                      Named("cll") = cll);
}