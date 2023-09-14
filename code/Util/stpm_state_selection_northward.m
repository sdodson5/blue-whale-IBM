function new_s = stpm_state_selection_northward(krill_loc, sst_loc,good_whales, par)

% gamma_K: trans prob based on krill
% gamma_E: trans prob based on environment
% gamma  : trans prob based on combo of krill and environ
% Size of gamma matrices depends on krill_loc and sst_loc inputs 
% Foraging probability specifically for the northward states S_1 and S_2


% Rows are whales, columns are locations 
alpha = par.envir_pref(1) + par.envir_pref(2).*(sst_loc + par.envir_pref(3));
gamma_E = exp(alpha)./(exp(alpha) + 1); 
gamma_K = (1 + exp(par.krill_pref(1).*(-krill_loc + par.krill_pref(2)))).^(-1);  % For states 1,2

gamma   = (par.w(1).* gamma_E + par.w(2).*gamma_K)./par.w(3);

new_s = 1 + double( rand(length(good_whales),1) < gamma); % Select state based on transition rate
