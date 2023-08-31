function [calling_whales, deltaCI] = find_calling_whales(sex_idx, good_whale_vec, S,par,call_prcnt)

% Return indices of calling whales and their song message (deltaCI)
% S = state vector for last par.rate time steps


calling_whales = zeros(par.numWhales,1);
deltaCI = zeros(par.numWhales,1);

% Probablistically set which male whales are calling
tmp_call = rand( length(sex_idx),1 );
tmp_call( tmp_call > call_prcnt ) = 0;
tmp_call( tmp_call > 0 ) = 1;
calling_whales(sex_idx) = good_whale_vec(sex_idx).*tmp_call; % 0 = no call, 1 = call. Only good_whales can call

% Determine DeltaCI: slideing window approach on the past par.rate
% behavioral states

x = sum( mod(S(calling_whales>0,:)-1,2), 2 )./par.rate;
deltaCI(calling_whales>0) =  par.call_singal_sent(1).*tanh(par.call_singal_sent(2).*(x-par.call_singal_sent(3)));

