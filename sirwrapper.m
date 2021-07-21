%% Parameterize and graph a SIR model of a zoonosis with an intermediate host
%  Kate Royce
%  September 2020

clc
%% Set parameters for other compartments

% all parameters are per day unless otherwise noted
b_r = 0.0625;             % wild birth rate 
m_r = 0.0625;             % wild background death rate                                 
b_h = 0.0118;             % human birth rate
m_h = 0.0118;             % human background death rate

%SARS and COVID
% beta_r = 0.25;            % wild transmission rate
% gamma_r = 0.14;           % wild recovery rate; SARS R0_w = 1.78
% gamma_h = 0.035;          % human recovery rate

%MERS               
beta_r = 0.08; % MERS is NOT highly transmissible in reservoir, but mild illness 
gamma_r = 0.14; 
gamma_h = 0.035;          % 1/day (28.4 day recovery period)



%% Set parameters for intermediate host compartment

T = readtable('MERS_species_params.xlsx');
h = height(T);
T.Beta_h = zeros(h, 1); %8
T.Final_Ih = zeros(h, 1); %9
T.Time_To_Infection = zeros(h, 1); %10
T.Max_Humans_Infected = zeros(h,1); %11
T.Time_To_Max = zeros(h, 1); %12
T.R0 = zeros(h, 1); %13
T.Rank = zeros(h,1); %14

C = table2cell(T(:,1), zeros(h,1), zeros(h,1), zeros(h,1)); %to store data

init_Sr = 0.9; % seed initial infection in reservoir at equilibrium Sw
init_Ir = 1-init_Sr; % (same for all IH simulations)

for i = 1:h
    % set IH species parameters
    
    %biological--measures biological factors that make IH similar to other species
    reservoir_similarity = T.res_sim(i); % unitless
    human_similarity = T.hum_sim(i); % unitless

    %contact--measures how frequently IH comes into contact with other species
    reservoir_contact = T.res_con(i); %contact/day
    human_contact = T.hum_con(i); %contact/day
    
    %transmission risk--measures if average contact between IH and other species is likely to transmit virus
    reservoir_risk = T.res_risk(i); % infections/contact
    human_risk = T.hum_risk(i);   % infections/contact

    % Model parameters
    b_i = 1./17; %birth and mortality rates are at replacement 
    m_i = 1./17; %individuals/day
    beta_iw = reservoir_similarity.*beta_r;         % transmission rate of wildtype in IH
    gamma_iw = gamma_r;                             % recovery rate from wildtype in IH, 1/day
    beta_im = beta_r;                               % transmission rate of mutant strain in IH
    gamma_im = gamma_r;                             % recovery rate from mutant in IH
    wild_transmission = reservoir_contact.*reservoir_risk;  % rate of spillover wild -> intermediate (infections/day)
    human_transmission = human_contact.*human_risk;         % rate of spillover intermediate -> humans
    mu = human_similarity.*human_contact.*human_risk;       % mutation rate (infections/day)
    beta_h = human_similarity.*beta_iw;                     % transmission rate in humans
    %the more intuitive possibility: 
    %beta_h = human_similarity.*beta_im;
    
    %% Run simulation

    % length of simulation
    time = [0,3000];

    % initial population values (proportions in each compartment) 
    population_vector = [init_Sr,init_Ir,0,    1,0,0,0,     1,0,0]; 

    % wrapper function for ODE solver
    sirWrapper = @(t,x) IH_sir(t,x, b_r, m_r, beta_r, gamma_r, b_i, m_i, beta_iw, gamma_iw, beta_im, gamma_im, wild_transmission, b_h, m_h, beta_h, gamma_h, human_transmission, mu);

    % solves system of ODEs
    [t,x] = ode45(sirWrapper, time, population_vector); 
    
    % stores data
    C{i,2} = t;
    C{i,3} = x;

    %% calculate stability of equilibria
    S_r = x(end, 1);
    I_r = x(end, 2);
    S_i = x(end, 4);
    I_i = x(end, 5);
    T_i = x(end, 6);
    S_h = x(end, 8);
    I_h = x(end, 9);

    T{i,8} = beta_h;
    T{i,9} = I_h;
    
    % index where %of humans susceptible drops below 0.99 for first time
    % (assuming spillover into population of at least 100 individuals, this
    % assures one 'real' individual has it)
    init_time = find(x(:,8)<0.99, 1, 'first'); 
    
    [max_Ih,max_time] = max(x(:,9));
    T{i,11} = max_Ih; % max proportion of humans infected
    
    if max_Ih > 10^(-2) % the epidemic establishes itself in humans
        T{i,10} = init_time;
        T{i,12} = max_time; % time to max proportion of humans infected
    else % no significant foothold in humans, so indicate with non-integer timestep
        T{i,10} = 0.5;
        init_time = 0.5;
        max_time = 0.5;
        T{i,12} = 0.5;
    end
    
    % calculate stability of equilibria
    % disease-free equilibria
    dfe = jacobian(m_r, beta_r, gamma_r, m_i, beta_iw, gamma_iw, beta_im, gamma_im, m_h, beta_h, gamma_h, wild_transmission, human_transmission, mu, 1,0,  1,0,0,  1,0);
    dfeTest = real(dfe)<0; %stable if all eigenvalues have real part < 0 ie this array is all 1s

    % endemic equilibria 
    ee = jacobian(m_r, beta_r, gamma_r, m_i, beta_iw, gamma_iw, beta_im, gamma_im, m_h, beta_h, gamma_h, wild_transmission, human_transmission, mu, S_r,I_r,    S_i,I_i,T_i,     S_h,I_h);
    eeTest = real(ee)<0; %stable if all eigenvalues have real part < 0 ie this array is all 1s

    % calculate R_0
    S_r_naive = b_r/m_r;
    S_i_naive = b_i/m_i;
    S_h_naive = b_h/m_h;
    A = nextgen(beta_r, gamma_r, m_r, beta_iw, gamma_iw, beta_im, gamma_im, m_i, mu, beta_h, gamma_h, m_h, wild_transmission, human_transmission, S_r_naive, S_i_naive, S_h_naive);     

    C{i,4} = A;
    global_r0 = max(eig(A));
    T{i,13} = global_r0;
    T{i,14} = (global_r0+max_Ih.*100)./2; %rank is average of global r0 and max percent infected humans
    
end 

[maximum,idx] = max(T.Rank);
species_name = string(T.Species(idx))
rows = (T.Max_Humans_Infected > .2 | T.R0 > 1.5); %change test based on pathogen!
cols = [1,9:14];
sort = T(rows,cols);
output = sortrows(sort, 7, 'descend')

produce_graph(C{idx, 1},C{idx,2},C{idx,3});
