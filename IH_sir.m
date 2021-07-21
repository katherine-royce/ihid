%% Simulates the progression of a zoonosis with an intermediate host with the goal of identifying the said host
% Kate Royce
% September 2020

function y = IH_sir(t, x, b_w, m_w, beta_w, gamma_w, b_d, m_d, beta_dw, gamma_dw, beta_dm, gamma_dm, p_d, b_h, m_h, beta_h, gamma_h, p_h, mu)
    y=zeros(10,1);
    
    % traditional model in wild compartment
    y(1)=b_w-beta_w*x(1)*x(2)-m_w*x(1);                                %S_w
    y(2)=beta_w*x(1)*x(2)-gamma_w*x(2)-m_w*x(2);                       %I_w
    y(3)=gamma_w*x(2)-m_w*x(3);                                        %R_w
    
    % 2 types of infections in domestic animals: caught from wild animals
    % OR transmissible to humans
    y(4)=b_d-beta_dw*x(4)*x(5)-p_d*x(2)*x(4)-beta_dm*x(4)*x(6)-m_d*x(4); %S_d
    y(5)=beta_dw*x(4)*x(5)+p_d*x(2)*x(4)-mu*x(5)-gamma_dw*x(5)-m_d*x(5); %I_d 
    y(6)=mu*x(5)+beta_dm*x(4)*x(6)-gamma_dm*x(6)-m_d*x(6);               %T_d
    y(7)=gamma_dw*x(5)+gamma_dm*x(6)-m_d*x(7);                           %R_d
    
    % traditional model in human compartment seeded by type T_d
    y(8)=b_h-beta_h*x(8)*x(9)-p_h*x(6)*x(8)-m_h*x(8);                  %S_h
    y(9)=beta_h*x(8)*x(9)+p_h*x(6)*x(8)-gamma_h*x(9)-m_h*x(9);         %I_h
    y(10)=gamma_h*x(9)-m_h*x(10);                                      %R_h
end
