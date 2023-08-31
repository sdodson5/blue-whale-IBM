function [x_star, y_star, ars_on] = ARS_4state(idx2, sst_tmp, krill_tmp, par, X_domain, Y_domain, ars_whales)
% Area restricted search to find the highest probability of foraging
% Take in index of current grid cell of whale
% (x_star, y_star) gives location (in meters) of highest prob of foraging

nRow = size(sst_tmp,1); 
ars_on = zeros(par.numWhales,1);  % Include all whales

x_star = zeros(size(idx2,1),1); % Only whales that do full ARS
y_star = zeros(size(idx2,1),1);

%% Regular ARS whales
r_x = zeros(size(idx2,1),2*par.search_rad+1);   % Rows are whales, each column is a neighboring cell
c_y = zeros(size(idx2,1),2*par.search_rad-1); 

tr_y = repmat(idx2(:,2) + par.search_rad, 1,2*par.search_rad+1); br_y = repmat(idx2(:,2) - par.search_rad,1, 2*par.search_rad+1);
lc_x = repmat(idx2(:,1) - par.search_rad, 1,2*par.search_rad-1); rc_x = repmat(idx2(:,1) + par.search_rad,1, 2*par.search_rad-1);

r_x(:,1) = idx2(:,1) - par.search_rad;
c_y(:,1) = idx2(:,2) - par.search_rad+1;

for j = 2:2*par.search_rad+1
    r_x(:,j) = r_x(:,j-1) + 1; 
       
end

for k = 2:2*par.search_rad-1
    c_y(:,k) = c_y(:,k-1) + 1;
        
end
    
nearby_idx_x = [r_x, r_x, lc_x, rc_x];  
nearby_idx_y = [tr_y, br_y, c_y, c_y];

nearby_idx = (nearby_idx_x -1).*nRow + nearby_idx_y;  % Cell indicies of search cells

nearby_sst   = sst_tmp(nearby_idx);     % Get the sst and krill on the search radius
nearby_krill = krill_tmp(nearby_idx); 

% Compute probability of foraging at neighboring grid cells
alpha = par.envir_pref(1) + par.envir_pref(2).*(nearby_sst + par.envir_pref(3));
gamma_E = exp(alpha)./(exp(alpha) + 1); 
gamma_K = (1 + exp(par.krill_pref(1).*(-nearby_krill + par.krill_pref(2)))).^(-1);

gammas_nearby   = (par.w(1).* gamma_E + par.w(2).*gamma_K)./par.w(3); % Foraging probability at nearby grid cells
gammas_nearby(isnan(gammas_nearby)) = 0;  % Set 0 probability of wanting to move onto land (cells with NaN)

ars_idx = find( max(gammas_nearby,[],2) - min(gammas_nearby,[],2) > par.search_thres);  % Whales that do ARS. Don't ARS if all directions are equal
ars_on(ars_whales(ars_idx)) = 1;

[~,mx_gamma] = max(gammas_nearby(ars_idx,:),[],2);  % Get indices of max gammas

if ~isempty(mx_gamma)
    
    max_idx = sub2ind(size(nearby_idx),ars_idx,mx_gamma);

    % Location of highest foraging prob
    x_star(ars_idx) = X_domain(nearby_idx(max_idx));
    y_star(ars_idx) = Y_domain(nearby_idx(max_idx)); 
    
end
   











