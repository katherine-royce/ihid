%% Compute next-generation matrix for IH mutation model
% Kate Royce
% October 2020

function A = nextgen(beta_w, gamma_w, m_w, beta_dw, gamma_dw, beta_dm, gamma_dm, m_d, mu, beta_h, gamma_h, m_h, p_d, p_h, S_w, S_d, S_h)
    F = [beta_w*S_w, 0, 0, 0;
        p_d*S_d, beta_dw*S_d, 0, 0;
        0, 0, beta_dm*S_d, 0;
        0, 0, p_h*S_h, beta_h*S_h];
    V = [gamma_w+m_w, 0, 0, 0;
        0, mu+gamma_dw+m_d, 0, 0;
        0, -mu, gamma_dm+m_d, 0;
        0, 0, 0, gamma_h+m_h];
    A = F*inv(V);
end