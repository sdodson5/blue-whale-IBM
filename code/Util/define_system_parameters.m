
% Define default system parameters and store in Matlab structure 'par'

% Basic parameters
par.rate        =  4;           % Samples/day
par.numWhales   = 2000;         % Number of whales
par.numDays     = 184;          % Number of days
par.search_rad  = 2;            % search radius for ARS (numbers of grid cells)
par.search_thres = 0.05;        % Conduct ARS if difference between min and max on search radius are greater than this threshold
par.w           = [1;2;3;2];      % Weights: Environmental, krill, normalization factor 
par.start_box   = [500000,700000,5/6,800000,800000]; % Information for location of start box

start_month_val = '07'; % Starting month/day: for the whales entering the domain 5 = May
final_month_val = '07'; % Ending month/day for whales entering the domain. Start dates are uniformly distributed between these two date


% Calling parameters
max_call_rad = 125; % Max call radius (km)
par.call_prcnt = [0.05,0.3,0,720]; %[min, max, min timstep, max time step]
par.A0 = 180;
par.call_singal_sent = [1;8;0.38]; 
par.maxCallDist = max_call_rad*1000; % Maximum call distance (meters)

% Southward migration parameters

% Parameter order: c_1, c_2, number of days to average over, number of days
% delay before southward migration allowed, c_3, c_4
par.state3_pref = [15, 0.2, 10, 20, 0.22, -10, 10]; 
par.x_star3 = 800000; % x-coordinate: state 3 heading
par.y_star3 = 3000; % y-coordinate: state 3 heading

% Movement parameters (same as in Dodson et al (2020), Ecological Modelling)
mu = [22200, 6300, 22200, 6300 ];  % Rate 4 mean step lengths: s1, s2, s3, s4
sigma = [12900, 5820, 12900, 5820]; % Rate 4 std dev step lengths: s1, s2, s3, s4

kappa       = [10, 3, 10, 3];          % Std dev turning angle: s1, s2, s3, s4
kappa_mu    = [0, pi, 0, pi];          % Mean turning angle: s1, s2, s3, s4

par.krill_pref  = [5, 0.3, 0.6];         % Krill preferences: [steepness of transtion, krill threshold]
par.envir_pref  = [-1,-0.2, -16];   % Environmental (SST) preferences for the SDM 

% Turn into vectors: helpful for selection of behavioral states 
par.kappa      = repmat(kappa,par.numWhales,1)     ;
par.kappa_mu   = repmat(kappa_mu,par.numWhales,1)  ;
par.mu         = repmat(mu,par.numWhales,1)        ;
par.sigma      = repmat(sigma,par.numWhales,1)     ;

