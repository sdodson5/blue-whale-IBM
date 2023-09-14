function [X, Y, sst, krill, s, forage_eff] = blueWhale_IBM_nullModels(par,sst_data,krill_data,grid_pars,mig_strat)
% Simulates NWHALES whales with ARS search
% Once whales exit, they cannot reenter the domain
% Monitor & record foraging efficiency
% Null model 1 mig_strat = 'no_mig': No migration
% Null model 2 mig_strat = 'doy': Southward migration date pre-selected from a normal distribution

N = par.rate*par.numDays;
numWhales = par.numWhales; % Rename so I don't have to type par each time

% (x,y) values - Values of cells are the X,Y-coordinates
X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); 
X_domain = repmat(X_domain(1:end-1),grid_pars.numX,1);
Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); 
Y_domain = repmat(Y_domain(1:end-1)',1,grid_pars.numY);

% Preallocate matrices: Rows are whales, columns are time steps
X       = zeros(numWhales,N);   % X coordinate (in meters from south-west corner of domain)
Y       = zeros(numWhales,N);   % Y coordinate (in meters from south-west corner of domain)
sst     = zeros(numWhales,N);   % SST
krill   = zeros(numWhales,N);   % Krill
s       = zeros(numWhales,N);   % State (1 or 2)
direc   = zeros(numWhales,1);   % Direction whale moved to get to current location
steps   = zeros(numWhales,1);   % Step lengths
turns   = zeros(numWhales,1);   % Turning angles
good_whale_vec = zeros(numWhales,N);    % Which whales are in the ocean
forage_eff = zeros(numWhales,N);

doy =  par.doy_start:1/par.rate:(par.doy_start+par.numDays); % Day of year vector

s2_matrix = zeros(numWhales,N);

x_new = zeros(numWhales,1); 
y_new = zeros(numWhales,1);

state2_correction_angle = pi; % Mean turning angle for S2 is pi. If whale conducts ARS in S2, need to add pi to angle selected to counteract this

% Initialize varaibles
s(:,1)     = 1;               % Start all whales in transit
direc(:)   = pi/2;            % State all whale pointing due north


% Select migration dates (if necessary)
switch mig_strat
    case 'doy' % Normally distributed based on day of year

        migration_dates = normrnd(par.state3_pref(1),par.state3_pref(2),[par.numWhales,1]); % Pre-selected migration dates for each whale
end

% case 'no_mig' has no migration dates

%% Sample initial condition from bounding box
sst_on = ones(numWhales,1);     % 1 if sst has not been selected 

%  initial condition sampled from areas of climatologically high krill:
%  previously identified based on ROMS data
X(:,1) = par.start_box(1) + (par.start_box(2) - par.start_box(1)).* rand(numWhales,1); % X coord
ymin = -par.start_box(3).*X(:,1) + par.start_box(4);
Y(:,1) = ymin + (par.start_box(5) - ymin).* rand(numWhales,1); % Y coord: depends on the X-coordinate

  
idx       = coordinateToGridCell([X(:,1),Y(:,1)],grid_pars);            % Maps (X,Y) coording in meters to matrix (column, row). Order of idx is (column, row) = (X,Y)!
idx_cell  = sub2ind([grid_pars.numX,grid_pars.numY],idx(:,2),idx(:,1)); % Maps matrix (row,column) to cell number

sst_tmp     = sst_data(:,:,par.doy_start);
krill_tmp   = krill_data(:,:,par.doy_start);

sst_loc   = sst_tmp( idx_cell );    % SST at whale locations
sst_on    = sst_on.*isnan(sst_loc); % 1 if whale at a location with no sst (i.e. on land)

while sum(sst_on)  % Resample for whales with no SST (on land)
    
    bad_idx = find(sst_on > 0);  % Whales that are on land
    
    % Resample initial conditions
    X(bad_idx,1) = par.start_box(1) + (par.start_box(2) - par.start_box(1)).* rand(length(bad_idx),1);
    ymin = -par.start_box(3).*X(bad_idx,1) + par.start_box(4);
    Y(bad_idx,1) = ymin + (par.start_box(5) - ymin).* rand(length(bad_idx),1); % Y coord: depends on the X-coordinate
   
    idx(bad_idx,:)     = coordinateToGridCell([X(bad_idx,1),Y(bad_idx,1)],grid_pars); % Order of idx is (column, row) = (X,Y)!
    idx_cell(bad_idx)  = sub2ind([grid_pars.numX,grid_pars.numY],idx(bad_idx,2),idx(bad_idx,1));

    % Determine if new locatio is on land
    sst_loc(bad_idx)   = sst_tmp( idx_cell(bad_idx) );
    sst_on    = sst_on.*isnan(sst_loc);
    
end

% Record sst and krill
sst(:,1)   = sst_loc;
krill(:,1) = krill_tmp( idx_cell );

% Sample steps and turns
steps(:) = gamrnd(par.mu(s(:,1)).^2./par.sigma(s(:,1)).^2, par.sigma(s(:,1)).^2./par.mu(s(:,1)));
turns(:) = mod(vmrand(par.kappa_mu(s(:,1)),par.kappa(s(:,1))),2*pi);

%% Generate rest of sequence

run_thresh = 2*N;  % To make sure an infinite while loop isn't created
good_whales = find(sst_on == 0); % Whales in the ocean 

for k = 2:N  % Time step

    sst_on = ones(numWhales,1);  % 1 if have not yet selected sst/krill for new whale location
    
    %% Generate new location for whales that were in ocean on last time step
    x_new(good_whales) = cos(direc(good_whales) + turns(good_whales)); 
    y_new(good_whales) = sin(direc(good_whales) + turns(good_whales));
    
    X(good_whales,k) = X(good_whales,k-1) + steps(good_whales).*x_new(good_whales);
    Y(good_whales,k) = Y(good_whales,k-1) + steps(good_whales).*y_new(good_whales);
    
    % Determine if whales have stepped out of ROMS domain (different than being beached on land)
    [whale_in_bounds, good_whales] = find_whales_in_bounds(X(good_whales,k), Y(good_whales,k),grid_pars,par,good_whales);  % Function to find which whales are still in bounds after most recent time step. 
    
    % Whale_in bounds is vector of {0,1}'s of length numWhales, good whales gives indices of whales in bounds
    good_whale_vec(good_whales,k) = 1; % Mostly for debugging purposes
    
    % Check to make sure there are still some whales in the domain
    if isempty(good_whales)
       disp(['All Whales Out of Bounds on Time Step ' num2str(k)]); 
       break;
        
    end
    
    % Extract SST/Krill at new location
    idx = zeros(numWhales,2); idx_cell = zeros(numWhales,1);
    sst_loc = zeros(numWhales,1);
    
    idx(good_whales,:)      = coordinateToGridCell([X(good_whales,k),Y(good_whales,k)],grid_pars); 
    idx_cell(good_whales)   = sub2ind([grid_pars.numX,grid_pars.numY],idx(good_whales,2),idx(good_whales,1));
    sst_tmp   = sst_data(  :,:,floor(k/par.rate) + par.doy_start );
    krill_tmp = krill_data(:,:,floor(k/par.rate) + par.doy_start );
    
    sst_loc(good_whales)   = sst_tmp( idx_cell(good_whales) );
    sst_on    = whale_in_bounds.*sst_on.*isnan(sst_loc);  % 1 corresponds to whale on land - only consider whales that are still in bounds
    
    % Resample any bad locations (beached whales)
    run_times = 0;
    while (sum(sst_on) && run_times < run_thresh)  % as long as sst hasn't been selected, and run times < threshold value (to ensure finite loop)
        run_times = run_times + 1;
        
        bad_idx = find(sst_on > 0); % whales on land
        idx_bad_whales_s = sub2ind(size(par.mu),bad_idx,s(bad_idx,k-1) ); % Find the cell index
    
        % Select new steps and turns
        steps(bad_idx) = gamrnd(par.mu(idx_bad_whales_s).^2./par.sigma(idx_bad_whales_s).^2, par.sigma(idx_bad_whales_s).^2./par.mu(idx_bad_whales_s));
         if k > ceil(N/4) % Add northward pull in first quarter of the simulation/year
            turns(bad_idx) = sign(-1 + 2*rand(length(bad_idx),1) )*pi/2 + turns(bad_idx); % Add or subtract 90 degrees with prob 1/2 each 
        else
            turns(bad_idx) = pi/2 + turns(bad_idx); % Add pi/2 (90 degree turn northward) during first quarter of simulation. 
        end
        
        % Calculate new positions
        x_new(bad_idx) = cos(direc(bad_idx) + turns(bad_idx));
        y_new(bad_idx) = sin(direc(bad_idx) + turns(bad_idx));
        X(bad_idx,k) = X(bad_idx,k-1) + steps(bad_idx).*x_new(bad_idx);
        Y(bad_idx,k) = Y(bad_idx,k-1) + steps(bad_idx).*y_new(bad_idx);
        
        idx(bad_idx,:) = coordinateToGridCell([X(bad_idx,k),Y(bad_idx,k)],grid_pars); 
        idx_cell(bad_idx)  = sub2ind([grid_pars.numX,grid_pars.numY],idx(bad_idx,2),idx(bad_idx,1));

        sst_loc(bad_idx)   = sst_tmp( idx_cell(bad_idx) );
        sst_on    = whale_in_bounds.*sst_on.*isnan(sst_loc); 
    end 
    
    % record sst and krill
    sst(good_whales,k)   = sst_loc(good_whales);
    krill(good_whales,k) = krill_tmp( idx_cell(good_whales) );
    
    
    
    
   %% Select new behavioral states

    switch mig_strat
        case 'doy'
            % List whales that are in the southward migration states
            south_mig_vec = find(migration_dates <= ceil(doy(k)) ); % Find indices of whales whose assigned migration dates are <= the current date
            
            south_whales = intersect(good_whales,south_mig_vec); % Current southward migrating whales
            north_whales = setdiff(good_whales,south_mig_vec);   % Current northward migrating whales
            
            % State selection for northward whales: states = 1,2 (was easier to separate the states)
            s(north_whales,k) =  stpm_state_selection_northward(krill(north_whales,k) , sst(north_whales,k), north_whales, par);
            
             % State selection for southward whales: states = 3,4
            s(south_whales,k) =  stpm_state_selection_southward(krill(south_whales,k) , sst(south_whales,k), south_whales, par);
            
        case 'no_mig'
            % State selection for states = 1,2 (no mig only has 1 & 2)
            s(good_whales,k) =  stpm_state_selection_northward(krill(good_whales,k) , sst(good_whales,k), good_whales, par);
    end
   
    s2_matrix(s(:,k)==2,k) = 1; % 1 if whale is in state 2
    
    % Find foraging efficiency history
    if k > par.state3_pref(4)*par.rate % Need history in order to calculate quantity
        forage_eff(:,k) = sum(krill(:,k-par.state3_pref(3)*par.rate:k).*s2_matrix(:,k-par.state3_pref(3)*par.rate:k),2)./(par.state3_pref(3)*par.rate);     
    end
    
    %% Determine next movement update

    % Conduct ARS
    x_star = zeros(numWhales,1); y_star = zeros(numWhales,1);
    
    s12_whales = union(find(s(:,k) == 1),find(s(:,k) ==2));     % whales in states 1 and 2 
    ars_whales = intersect(good_whales( idx(good_whales,1) > par.search_rad + 1), good_whales( idx(good_whales,2) > par.search_rad + 1)  );  % Remove whales whos search radius extends out of the domain.
    s12_whales = intersect(ars_whales,s12_whales);
    [x_star(s12_whales), y_star(s12_whales), ars_on] = ARS_4state(idx(s12_whales,:), sst_tmp,krill_tmp, par, X_domain, Y_domain, s12_whales);  % Location of highest prob of foraging
      
    x_new(ars_on==1) = x_star(ars_on==1) - X(ars_on==1,k); 
    y_new(ars_on==1) = y_star(ars_on==1) - Y(ars_on==1,k);
    x_new(ars_on==1) = x_new(ars_on==1)./(sqrt(x_new(ars_on==1).^2 + y_new(ars_on==1).^2));
    
    % Point state 3 whales toward breeding grounds
    x_new(s(:,k)==3) = par.x_star3 - X(s(:,k)==3,k);
    y_new(s(:,k)==3) = par.y_star3 - Y(s(:,k)==3,k);
    x_new(s(:,k)==3) = x_new(s(:,k)==3)./(sqrt(x_new(s(:,k)==3).^2 + y_new(s(:,k)==3).^2));
    
    % If the whales are in state 2, the angle is centered around pi - need
    % to update to 0 if the ARS is on. Do this by adding pi if whale in
    % state 2 and ARS on.
    state2_correction_on = zeros(par.numWhales,1);
    state2_correction_on(s12_whales) = (s(s12_whales,k)-1).*ars_on(s12_whales);
    state2_correction_on = state2_correction_on(good_whales);
        
    direc(good_whales) = sign(y_new(good_whales)).*acos(x_new(good_whales)); % Select the direction the whales moved to get from k-1 to kth location. If the whales conduct ARS, this is updated to point toward area of highest krill
    
    % Select new step lengths and turning angles    
    idx_good_whales_s    = sub2ind(size(par.mu),good_whales,s(good_whales,k) );  % Tells cell index of par.mu, sigma, kappa, kappa_mu to select for each whale. 
    steps(good_whales) = gamrnd(par.mu(idx_good_whales_s).^2./par.sigma(idx_good_whales_s).^2, par.sigma(idx_good_whales_s).^2./par.mu(idx_good_whales_s));
    turns(good_whales) = mod(vmrand(par.kappa_mu(idx_good_whales_s),par.kappa(idx_good_whales_s)),2*pi) + state2_correction_on.*state2_correction_angle;
    
     
end








