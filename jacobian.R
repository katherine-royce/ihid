#' Return eigenvalues of Jacobian of IH mutation model 
#' 
#' Computes the eigenvalues of the Jacobian matrix of the IH mutation model given original parameters and populations.
#' 
#' @return e An ordered list of the eigenvalues of the Jacobian
#' @export 

jacobian <- function(m_w, beta_w, gamma_w, m_d, beta_dw, gamma_dw, beta_dm, gamma_dm, m_h, beta_h, gamma_h, p_d, p_h, mu, S_w,I_w,    S_d,I_d,T_d,     S_h,I_h){

row1 <- c(beta_w*I_w-m_w, -beta_w*S_w,               0,           0,                                    0,                             0,            0,    0,                       0,                            0)
row2 <- c(beta_w*I_w,     beta_w*S_w-(gamma_w+m_w),  0,           0,                                    0,                             0,            0,    0,                       0,                            0)
row3 <- c(0,              gamma_w,                   -m_w,        0,                                    0,                             0,            0,    0,                       0,                            0)
row4 <- c(0,              -p_d*S_d,                  0,           -beta_dw*I_d-beta_dm*T_d-p_d*I_w-m_d, -beta_dw*S_d,                  -beta_dm*S_d, 0,    0,                       0,                            0)
row5 <- c(0,              p_d*S_d,                   0,           beta_dw*I_d+p_d*I_w,                  beta_dw*S_d-(mu+gamma_dw+m_d), 0,            0,    0,                       0,                            0)
row6 <- c(0,              0,                         beta_dm*T_d, mu,                                   beta_dm*S_d-(gamma_dm+m_d),    0,            0,    0,                       0,                            0)
row7 <- c(0,              0,                         0,           0,                                    gamma_dw,                      gamma_dm,     -m_d, 0,                       0,                            0)
row8 <- c(0,              0,                         0,           0,                                    0,                             -p_h*S_h,     0,    -beta_h*I_h-p_h*T_d-m_h, -beta_h*S_h,                  0)
row9 <- c(0,              0,                         0,           0,                                    0,                             p_h*S_h,      0,    beta_h*I_h+p_h*T_d,      beta_h*S_h-(gamma_h+m_h),     0)
row10 <- c(0,              0,                         0,           0,                                    0,                             0,            0,    0,                       gamma_h,                  -m_h)
A <- rbind(row1, row2, row3, row4, row5,row6, row7, row8, row9, row10)

e <- sort(eigen(A)$values);
return(e)
}
