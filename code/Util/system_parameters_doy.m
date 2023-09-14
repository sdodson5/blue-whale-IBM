
% Define default parameters for the day of year null model

% Define system parameters and store in Matlab structure 'par'
par.rate        =  4;           % Samples/day
par.numWhales   = 2000;         % Number of whales
par.numDays     = 184;          % Number of days
par.search_rad  = 2;            % search radius for ARS (numbers of grid cells)
par.search_thres = 0.05;        % Conduct ARS if difference between min and max on search radius are greater than this threshold
par.w           = [1;2;3;2];      % Weights: Environmental, krill, normalization factor  <--- !!! Set to [1;0;1] for SST only, [0;1;1] for krill only
par.start_box   = [500000,700000,5/6,800000,800000];%[400000 700000 300000 600000];%[600000 900000 100000 200000];%[400000 1000000 50000 200000]; % min(X), max(X), min(Y), max(Y)
par.x_star3 = 800000;
par.y_star3 = 3000;
par.state3_pref = [310, 20,10, 20]; % mean migration date, std deviation, days to average foraging eff over, migration hold

start_month_val = '07'; % Starting month/day: for the whales entering the domain 5 = May
final_month_val = '07'; % Ending month/day for whales entering the domain. 


mu = [22200, 6300, 22200, 6300 ];  % Rate 4 mean step lengths: s1, s2, s3, s4
sigma = [12900, 5820, 12900, 5820]; % Rate 4 std dev step lengths: s1, s2, s3, s4


kappa       = [10, 3, 10, 3];          % Std dev turning angle: s1, s2, s3, s4
kappa_mu    = [0, pi, 0, pi];          % Mean turning angle: s1, s2, s3, s4

par.krill_pref  = [5, 0.3, 0.6];    % Krill preferences: [steepness of transtion, krill threshold]
par.envir_pref  = [-1,-0.2, -16];   % Environmental (SST) preferences for the SDM 

% Turn into vectors: helpful for selection of behavioral states 
par.kappa      = repmat(kappa,par.numWhales,1)     ;
par.kappa_mu   = repmat(kappa_mu,par.numWhales,1)  ;
par.mu         = repmat(mu,par.numWhales,1)        ;
par.sigma      = repmat(sigma,par.numWhales,1)     ;






