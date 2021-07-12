#' Compute next-generation matrix for IH mutation model
#' 
#' This function computes the next-generation matrix for the IH SIR model given the original parameters and susceptibles.
#' 
#' @return A The next-generation matrix
#' @export 

nextgen <- function(beta_w, gamma_w, m_w, beta_dw, gamma_dw, beta_dm, gamma_dm, m_d, mu, beta_h, gamma_h, m_h, wild_transmission, human_transmission, S_w, S_d, S_h){
F <- rbind(c(beta_w*S_w, 0, 0, 0),
     c(wild_transmission*S_d, beta_dw*S_d, 0, 0),
     c(0, 0, beta_dm*S_d, 0),
     c(0, 0, human_transmission*S_h, beta_h*S_h))
V <- rbind(c(gamma_w+m_w, 0, 0, 0),
     c(0, mu+gamma_dw+m_d, 0, 0),
     c(0, -mu, gamma_dm+m_d, 0),
     c(0, 0, 0, gamma_h+m_h))
A = F*solve(V) # FV^(-1)
return(A)
}