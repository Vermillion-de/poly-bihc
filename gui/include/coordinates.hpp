#include "types.hpp"

tuple<MatrixXd, MatrixXd> method1__G1n__G1_( 
    const MatrixXd& g1n, const MatrixXd& g1_, 
    const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig,
    const MatrixXd& G1n, const MatrixXd& G1_, 
    const MatrixXd& G2n, const MatrixXd& G2_)
{
    const MatrixXd A = inv(Ig - g1n) * g1_;
    const MatrixXd B = g2n * A + g2_;
    const MatrixXd C = G2n * A + G2_;
    const MatrixXd D = C * inv(B);
    return {D * (Ig - g1n), -D * g1_};
}

tuple<MatrixXd, MatrixXd> method2__G1n__G1_( 
    const MatrixXd& g1n, const MatrixXd& g1_, 
    const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig,
    const MatrixXd& G1n, const MatrixXd& G1_, 
    const MatrixXd& G2n, const MatrixXd& G2_)
{
    const MatrixXd A = inv(g1_) * (Ig - g1n); // d = Ac
    const MatrixXd B = g2_ * A + g2n;
    const MatrixXd C = G2n + G2_ * A;
    const MatrixXd D = C * inv(B);
    return {D * (Ig - g1n), -D * g1_};
}

tuple<MatrixXd, MatrixXd> bih_reg_02(
    const MatrixXd& g1n, const MatrixXd& g1_, 
    const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig,  
    const MatrixXd& G1n, const MatrixXd& G1_, 
    const MatrixXd& G2n, const MatrixXd& G2_, const vector<real> lamb)
{
    const int c_v = g1n.cols(), c_n = g1_.cols();
    const int c_l = g2n.cols(), c_m = g2_.cols();
    const MatrixXd mv_0 = MatrixXd::Zero(g1n.rows(), g1n.cols()); 
    const MatrixXd mn_0 = MatrixXd::Zero(g1_.rows(), g1_.cols()); 
    const MatrixXd m_vn = optimize({g1n-Ig,  g2n}, {g1_,   g2_}, {c_l, c_m}, 
                                   {mv_0, Ig-g1n}, {mn_0, -g1_}, {c_v, c_n}, lamb);
    return vnlm_to_vn({G1n, G1_, G2n, G2_}, m_vn);
}
