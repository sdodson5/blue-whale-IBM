function s_new = stpm_state_selection_social_personal(krill_loc, sst_loc, par,s_old, forage_eff, sum_calls,south_flag,mig_strat)

% gamma is a numWhale x 4 matrix with transition probabilities.
% Dependent upon the previous state (s)
% Probability of transition to state 3 (south migration) foraging efficiency and social calls - MULTIPLICATION OF
% THE TWO!

nWhales = length(krill_loc);

stateMatrix34 = repmat(ceil(s_old./2) - 1,1,4); % 1 if in state 3,4 0 if states 1,2
stateMatrix12 = ones(size(stateMatrix34)) - stateMatrix34;

% Probability of transitioning to S3: strategy dependent
sigma_FE = (1 + exp(par.state3_pref(1).*(forage_eff - par.state3_pref(2)))).^(-1); % Based on foraging efficiency
sigma_S = (1 + exp(par.state3_pref(5).*(sum_calls - par.state3_pref(6)))).^(-1);  % Based on social calls. sum_calls gives sum of total signals received/number of signals received = average of signals received

switch mig_strat

    case 'personal'
    
        s3_trans = south_flag.*sigma_FE; % Probability of transition to S3: only FE
        w1 = (1 - (s3_trans))./(1 + par.w(4));
        w2 = par.w(4).*w1;


    case 'personal_social'

        s3_trans = south_flag.*sigma_FE.*sigma_S; % Probability of transition to S3: FE + social
        w1 = (1 - (s3_trans))./(1 + par.w(4));
        w2 = par.w(4).*w1;


    case 'social'

        s3_trans = south_flag.*sigma_S; % Probability of transition to S3: social
        w1 = (1 - (s3_trans))./(1 + par.w(4));
        w2 = par.w(4).*w1;

end

alpha = par.envir_pref(1) + par.envir_pref(2).*(sst_loc + par.envir_pref(3));
gamma_E = exp(alpha)./(exp(alpha) + 1);  % Environmental influence
gamma_K1 = (1 + exp(par.krill_pref(1).*(-krill_loc + par.krill_pref(2)))).^(-1); % Krill influence for states 1,2

gamma_K2 = (1 + exp(par.krill_pref(1).*(-krill_loc + par.krill_pref(3)))).^(-1);  % Krill influence for states 3,4

% Transition rates for whales currently in state 1 or 2
P2_12 = w1.*gamma_E + w2.*gamma_K1; % Transition to state 2 (forage)
P1_12 = w1 + w2 - P2_12;            % Transition to state 1 (transit)
P3_12 = s3_trans;                    % Transition to state 3 (south transit)
P4_12 = zeros(nWhales,1);% Transition to state 4 (south forage)

% Transition rates for whales currently in state 3 or 4
P2_34 = zeros(nWhales,1);  % Transition to state 2 (north forage)
P1_34 = zeros(nWhales,1);  % Transition to state 1 (north transit)
P3_34 = 1 - gamma_K2;                  % Transition to state 3 (south transit) 
P4_34 = gamma_K2;                      % Transition to state 4 (south forage) 

% Form full transition matrices
TP1 = [P1_12, P2_12, P3_12, P4_12]; % Transition probabilities for whales currently in states 1 and 2
TP2 = [P1_34, P2_34, P3_34, P4_34]; % Transition probabilities for whales currently in states 3 and 4

transProb = stateMatrix12 .* TP1 + stateMatrix34 .* TP2;  % Selects which transition matrix is relevant for each whale. 

% Select new states for each agent.
s_new = zeros(nWhales,1);

for j = 1:nWhales
    s_new(j) = randsample(4,1,'true',transProb(j,:));  % Select the state for each whale
            
end
   


