#include <RcppEigen.h>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;               	      // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers


Eigen::MatrixXd appendOne(MatrixXd x) {
  //Returns column matrix of original with 1.0 appended
  MatrixXd temp_x = x;

  temp_x.resize(temp_x.rows() + 1, 1);

  for (int i = 0; i < temp_x.rows(); i++) {
    temp_x(i, 0) = x(i, 0);
  }
  temp_x(temp_x.rows() - 1, 0) = 1.0;
  return temp_x;
}

Eigen::MatrixXd sorted_matrix_by_feat(Eigen::MatrixXd A, int feat) {

  std::vector<Eigen::VectorXd> vec;
  for (size_t i = 0; i < A.rows(); ++i)
    vec.push_back(A.row(i));

  std::sort(vec.begin(), vec.end(),
            [feat](const Eigen::VectorXd& lhs, const Eigen::VectorXd& rhs) {
            return lhs(feat) < rhs(feat);
  });

  for (size_t i = 0; i < A.rows(); ++i)
    A.row(i) = vec[i];

  return A;
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd update_A_inv_C(Eigen::MatrixXd a, Eigen::MatrixXd new_x, bool leftNode) {
  MatrixXd temp_x = appendOne(new_x);

   //Initilize z_K
   MatrixXd z_K = a * temp_x;

   //Update A using Sherman–Morrison formula corresponding to right or left side
   if (leftNode) {
     MatrixXd g_K = (z_K * z_K.transpose()) / (1 + (temp_x.transpose() * z_K)(0,0));
     return a - g_K;
   } else {
     MatrixXd g_K = (z_K * z_K.transpose()) / (1 - (temp_x.transpose() * z_K)(0,0));
     return a + g_K;
   }
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd update_S_k(Eigen::MatrixXd prev, Eigen::MatrixXd next_obs, double next_y, bool left) {
  if (left) {
    return prev + (next_y * (appendOne(next_obs)));
  } else {
    return prev - (next_y * (appendOne(next_obs)));
  }
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd update_G_k(Eigen::MatrixXd prev, Eigen::MatrixXd next_obs, bool left) {
  MatrixXd t = appendOne(next_obs);
  if (left) {
    return prev + (t * t.transpose());
  } else {
    return prev - (t * t.transpose());
  }
}



double compute_RSS(Eigen::MatrixXd A_r, Eigen::MatrixXd A_l,
                   Eigen::MatrixXd S_r, Eigen::MatrixXd S_l,
                   Eigen::MatrixXd G_r, Eigen::MatrixXd G_l) {

  return ((S_l.transpose() * A_l * G_l * A_l * S_l)(0,0) +
          (S_r.transpose() * A_r * G_r * A_r * S_r)(0,0) -
          (2 * S_l.transpose() * A_l * S_l)(0,0) -
          (2 * S_r.transpose() * A_r * S_r)(0,0));
}

//' @export
// [[Rcpp::export]]
Eigen::VectorXd compute_Best_RSS(Eigen::Map<Eigen::MatrixXd> training_data,
                                 double lambda,
                                 int splitting_feat,
                                 Eigen::Map<Eigen::VectorXd> linear_features) {


    //Take only given features when doing the linear fit
    MatrixXd temp = MatrixXd(training_data.rows(), linear_features.size() + 1);
    VectorXd v = VectorXd(1);
    v << (training_data.rows() - 1);

    VectorXd used_columns(linear_features.rows() + 1);
    used_columns << linear_features, v;

    for (size_t i = 0; i < used_columns.size(); i++) {
        temp.col(i) = training_data.col(linear_features(i));
    }

    training_data = temp;

    //Sort training_data by splitting_feat
    training_data = sorted_matrix_by_feat(training_data, splitting_feat);

    int best_obs;
    int current_obs = 1;
    double best_RSS = std::numeric_limits<double>::infinity();
    double current_RSS;


    int n = training_data.rows();
    int d = training_data.cols() - 1;
    MatrixXd first_x = training_data.block(0,0,1,d);
    MatrixXd appended_first = appendOne(first_x.transpose());

    //Initialize left A inverse
    MatrixXd i = MatrixXd::Identity(d + 1, d + 1);
    i(d,d) = 0;
    MatrixXd A_l = ((appended_first*(appended_first.transpose())) + lambda * i).inverse();


    //Initialize right A inverse
    MatrixXd A_r = lambda * i;
    for (int k = 1; k < n; k++) {
      MatrixXd temp = training_data.block(k,0,1,d);
      MatrixXd next_x = appendOne(temp);
      A_r += (next_x*(next_x.transpose()));
    }

    A_r = A_r.inverse();


    //Initialize S_l and S_r
    MatrixXd S_l = training_data(0,d) * (appended_first.transpose());

    MatrixXd S_r = training_data(1,d) *
                    ((appendOne(training_data.block(1,0,1,d))).transpose());

    for (int k = 2; k < n; k++) {
        S_r += training_data(k,d) *
                    ((appendOne(training_data.block(k,0,1,d))).transpose());
    }

    //Initialize G_l and G_r

    MatrixXd G_l = appended_first*(appended_first.transpose());

    MatrixXd G_r = appendOne(training_data.block(1,0,1,d)) *
                  (appendOne(training_data.block(1,0,1,d)).transpose());

    for (int k = 0; k < n; k++) {
      G_r += appendOne(training_data.block(k,0,1,d)) *
        (appendOne(training_data.block(k,0,1,d)).transpose());
    }

    //Evaluate first RSS
    current_RSS = compute_RSS(A_r, A_l, S_r, S_l, G_r, G_l);

    //Evaluate RSS at each point

    for (int k = 1; k < n; k++) {
      MatrixXd current_x = training_data.block(k,0,1,d).transpose();
      double current_y = training_data(k,d);
      current_obs++;

      A_l = update_A_inv_C(A_l, current_x, true);
      A_r = update_A_inv_C(A_r, current_x, false);

      S_l = update_S_k(S_l, current_x, current_y, true);
      S_r = update_S_k(S_r, current_x, current_y, false);

      G_l = update_G_k(G_l, current_x, true);
      G_r = update_G_k(G_r, current_x, false);

      current_RSS = compute_RSS(A_r, A_l, S_r, S_l, G_r, G_l);

      if (current_RSS < best_RSS) {
        best_obs = current_obs;
      }
    }

    VectorXd results;
    results << best_obs, best_RSS;
    return results;
}


/*** R
set.seed(309814)
# feat <- matrix(rnorm(80), ncol = 4)
# lambda <- .2
# feat_current <- feat[-20, ]
# feat_current_c <- cbind(feat_current, 1)
# A_inv_current <- solve(t(feat_current_c) %*% feat_current_c +
#                          lambda * diag(rep(1, ncol(feat_current_c))))
# new_x <- feat[20, ]
# A_inv_algo <- update_A_inv_R(A_inv_current = A_inv_current, new_x = new_x,
#                          leftnode = TRUE)
# feat_c <- cbind(feat, 1)
#
# A_inv_truth <- solve(t(feat_c) %*% feat_c +
#                        lambda * diag(rep(1, ncol(feat_c)))) # truth
#
# l <- update_A_inv_C(A_inv_current, new_x, TRUE)
# l
# A_inv_algo
# A_inv_truth

a <- matrix(rnorm(20), ncol = 5)



a


solve(a)

  */
