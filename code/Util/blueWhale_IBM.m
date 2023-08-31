function [X, Y, sst, krill, s,meanSumCalls,forage_eff] = blueWhale_IBM(par,sst_data,krill_data,grid_pars,mig_strat)

% Simulates NWHALES whales 
% Once whales exit, they cannot reenter the domain
% Call signal based on the average over M days
% Call signal based on Inverse Square Law formula
% Increasing percentage of whales calling throughout the season: linear increase
% Southward migration based on foraging efficiency and social
% communications

N = par.rate*par.numDays;
numWhales = par.numWhales; % Rename so I don't have to type par each time

% (x,y) values - Values of cells are the X,Y-coordinates
X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); 
X_domain = repmat(X_domain(1:end-1),grid_pars.numX,1);
Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); 
Y_domain = repmat(Y_domain(1:end-1)',1,grid_pars.numY);

% Preallocate matrices: Rows are agents, columns are time steps
X       = zeros(numWhales,N);   % X coordinate (in meters from south-west corner of domain)
Y       = zeros(numWhales,N);   % Y coordinate (in meters from south-west corner of domain)
sst     = zeros(numWhales,N);   % SST
krill   = zeros(numWhales,N);   % Krill
s       = zeros(numWhales,N);   % State (1, 2, 3, or 4)
direc   = zeros(numWhales,1);   % Direction whale moved to get to current location
steps   = zeros(numWhales,1);   % Step lengths
turns   = zeros(numWhales,1);   % Turning angles
good_whale_vec = zeros(numWhales,N);    % Which whales are in the ocean
deltaCI  = zeros(numWhales, N);
calling_whales = zeros(numWhales,N);
sum_calls = zeros(numWhales,N);
meanSumCalls = zeros(numWhales,N); % mean sum calls over number of days
numWhalesHeard = zeros(numWhales,N);
forage_eff = zeros(numWhales,N);

s2_matrix = zeros(numWhales,N);

x_new = zeros(numWhales,1); 
y_new = zeros(numWhales,1);

state2_correction_angle = pi; % Mean turning angle for S2 is pi. If whale conducts ARS in S2, need to add pi to angle selected to counteract this

% Initialize varaibles
s(:,1)     = 1;               % Start all whales in transit
direc(:)   = pi/2;            % Start all whale pointing due north


sex = randi([0,1],numWhales,1);
sex_idx = find(sex);

south_flag = 0;

MB_hyrdophone = 1e5.*[5.654648091270412, 4.250000000000000]; % Location of MB hydrophone coordinates in terms of the ROMS data
%% Select initial condition 
sst_on = ones(numWhales,1);     % Flag to determine if initial conditions are valid: 1 if sst has not been selected 

%  initial condition sampled from areas of climatologically high krill:
%  previously identified based on ROMS data
X(:,1) = par.start_box(1) + (par.start_box(2) - par.start_box(1)).* rand(numWhales,1); % X coord
ymin = -par.start_box(3).*X(:,1) + par.start_box(4);
Y(:,1) = ymin + (par.start_box(5) - ymin).* rand(numWhales,1); % Y coord: depends on the X-coordinate
    
idx       = coordinateToGridCell([X(:,1),Y(:,1)],grid_pars);            % Maps (X,Y) coording in meters to matrix (column, row). Order of idx is (column, row) = (X,Y)!
idx_cell  = sub2ind([grid_pars.numX,grid_pars.numY],idx(:,2),idx(:,1)); % Maps matrix (row,column) to cell number

% SST & krill at current time step
sst_tmp     = sst_data(:,:,1);
krill_tmp   = krill_data(:,:,1);

sst_loc   = sst_tmp( idx_cell );    % SST at whale locations
sst_on    = sst_on.*isnan(sst_loc); % 1 if whale at a location with no sst (i.e. on land, land = NaN values)

while sum(sst_on)  % Resample for whales with no SST (on land)
    
    bad_idx = find(sst_on > 0);  % Whales that are on land
    
    % resample initial condition
    X(bad_idx,1) = par.start_box(1) + (par.start_box(2) - par.start_box(1)).* rand(length(bad_idx),1);
    ymin = -par.start_box(3).*X(bad_idx,1) + par.start_box(4);
    Y(bad_idx,1) = ymin + (par.start_box(5) - ymin).* rand(length(bad_idx),1); % Y coord: depends on the X-coordinate

    idx(bad_idx,:)     = coordinateToGridCell([X(bad_idx,1),Y(bad_idx,1)],grid_pars); % Order of idx is (column, row) = (X,Y)!
    idx_cell(bad_idx)  = sub2ind([grid_pars.numX,grid_pars.numY],idx(bad_idx,2),idx(bad_idx,1));

    % Determine if new location is on land
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
    
    % Determine if whales have stepped out of ROMS domain (different than being on land)
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
    sst_tmp   = sst_data(  :,:,1+floor(k/par.rate) );
    krill_tmp = krill_data(:,:,1+floor(k/par.rate) );
    
    sst_loc(good_whales)   = sst_tmp( idx_cell(good_whales) );
    sst_on    = whale_in_bounds.*sst_on.*isnan(sst_loc);  % 1 corresponds to whale on land - only consider whales that are still in bounds
    
    % Resample any bad locations (whales on land)
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
            turns(bad_idx) = pi/2 + turns(bad_idx); % Add pi/2 (90 degree turn northward ) during first quarter of simulation. 
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
    
    % record sst and krill at location of each agent
    sst(good_whales,k)   = sst_loc(good_whales);
    krill(good_whales,k) = krill_tmp( idx_cell(good_whales) );
    
%% Social behavior    
    if k > par.rate % Need some minimal history
        
        %% Call senders (calling whales)
        % Update percentage of agents calling based on date
        call_prcnt = ((par.call_prcnt(2) - par.call_prcnt(1))./par.call_prcnt(4)).*k + par.call_prcnt(1);
        
        % Select which male whales are singing and find their delta CI
        [calling_whales(:,k), deltaCI(:,k)] = find_calling_whales(sex_idx, good_whale_vec(:,k), s(:,k-par.rate:k-1),par,call_prcnt);        
        
        calling_idx = find(calling_whales(good_whales,k)); % Index of calling whales in the good_whales vec
        calling_idx2 = find(calling_whales(:,k)); % Index of calling whales among all whales 
        
        %% Call receivers (all whales)
        % Find distances between all good whales     
        distMat = squareform(pdist([X(good_whales,k),Y(good_whales,k); MB_hyrdophone],'euclidean')); % Distances in meters

        hydroDist = distMat(calling_idx,end);
        distMat = distMat(calling_idx,1:end-1); %  Restrict rows to calling whales, columns are all whales. 
        distMat(distMat > par.maxCallDist) = 0; % Set distances  = 0 if outside of maxCallDistance
        
        call_signal = (par.A0 - abs(20.*log10(1./distMat))).*repmat(deltaCI(calling_idx2,k),1,length(good_whales)); % INVERSE SQUARE LAW: A0 represents source level in dB
        call_signal(distMat==0) = 0; % Set the call signals of far away whales (or themselves) = 0

        % Determine average call signal of current time step
        tmp_call_sum = sum(call_signal~=0,1); % How many calling whales are heard - turns into logical array and then sums the number of non-zero entries
        if isempty(tmp_call_sum)
            tmp_call_sum = 0;
        end
        numWhalesHeard(good_whales,k) = tmp_call_sum;
        sum_calls(good_whales,k) = sum(call_signal); % Sum of all non-zero calls heard
            
        if k > par.state3_pref(7)*par.rate % Average over previous time steps
           
            tot_whales_heard = sum(numWhalesHeard(good_whales,k-par.state3_pref(7)*par.rate:k),2);
            tot_whales_heard(tot_whales_heard ==0) = 1; %  If no calling whales heard, set to 1 so that not dividing by zero below.
            meanSumCalls(good_whales,k) = sum(sum_calls(good_whales,k-par.state3_pref(7)*par.rate:k),2)./tot_whales_heard; % Find the average of the sum calls over the last par.state3_pref(7) days.
            
        end
                    

    end
    
  
    
    %% Select new behavioral states based on current conditions and migration strategy
    % migration strategy built into next function
    s(good_whales,k) = stpm_state_selection_social_personal(krill(good_whales,k), sst(good_whales,k), par,s(good_whales,k-1), forage_eff(good_whales,k-1), meanSumCalls(good_whales,k),south_flag,mig_strat);
    s2_matrix(s(:,k)==2,k) = 1; % 1 if whale is in state 2
    
    % Find foraging efficiency history
    if k > par.state3_pref(4)*par.rate % Hold off southward transition until the southward transition probabilities can be calculated (need enough history)
        
        forage_eff(:,k) = sum(krill(:,k-par.state3_pref(3)*par.rate:k).*s2_matrix(:,k-par.state3_pref(3)*par.rate:k),2)./(par.state3_pref(3)*par.rate);
        south_flag = 1;
        
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
    idx_good_whales_s  = sub2ind(size(par.mu),good_whales,s(good_whales,k) );  % Tells cell index of par.mu, sigma, kappa, kappa_mu to select for each whale. 
    steps(good_whales) = gamrnd(par.mu(idx_good_whales_s).^2./par.sigma(idx_good_whales_s).^2, par.sigma(idx_good_whales_s).^2./par.mu(idx_good_whales_s));
    turns(good_whales) = mod(vmrand(par.kappa_mu(idx_good_whales_s),par.kappa(idx_good_whales_s)),2*pi) + state2_correction_on.*state2_correction_angle;
    
     
end





