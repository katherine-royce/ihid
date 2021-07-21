%% Return eigenvalues of Jacobian of mutation model given initial params and equilibrium
% Kate Royce
% October 2020

function e = jacobian(m_w, beta_w, gamma_w, m_d, beta_dw, gamma_dw, beta_dm, gamma_dm, m_h, beta_h, gamma_h, p_d, p_h, mu, S_w,I_w,    S_d,I_d,T_d,     S_h,I_h)
    
    A = [beta_w*I_w-m_w, -beta_w*S_w,               0,           0,                                    0,                             0,            0,    0,                       0,                            0;
         beta_w*I_w,     beta_w*S_w-(gamma_w+m_w),  0,           0,                                    0,                             0,            0,    0,                       0,                            0;
         0,              gamma_w,                   -m_w,        0,                                    0,                             0,            0,    0,                       0,                            0;
         0,              -p_d*S_d,                  0,           -beta_dw*I_d-beta_dm*T_d-p_d*I_w-m_d, -beta_dw*S_d,                  -beta_dm*S_d, 0,    0,                       0,                            0;
         0,              p_d*S_d,                   0,           beta_dw*I_d+p_d*I_w,                  beta_dw*S_d-(mu+gamma_dw+m_d), 0,            0,    0,                       0,                            0;
         0,              0,                         beta_dm*T_d, mu,                                   beta_dm*S_d-(gamma_dm+m_d),    0,            0,    0,                       0,                            0;
         0,              0,                         0,           0,                                    gamma_dw,                      gamma_dm,     -m_d, 0,                       0,                            0;
         0,              0,                         0,           0,                                    0,                             -p_h*S_h,     0,    -beta_h*I_h-p_h*T_d-m_h, -beta_h*S_h,                  0;
         0,              0,                         0,           0,                                    0,                             p_h*S_h,      0,    beta_h*I_h+p_h*T_d,      beta_h*S_h-(gamma_h+m_h),     0;
         0,              0,                         0,           0,                                    0,                             0,            0,    0,                       gamma_h,                  -m_h];
    e = sort(eig(A)); 
end
