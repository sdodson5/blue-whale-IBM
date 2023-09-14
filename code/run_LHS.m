% Latin Hypercube sampling on blue whale IBM simulation
% Stephanie Dodson: sdodson@colby.edu
% September 2023
% For parameter testing and sensitivity analysis: runs one year for all
% parameters before loading in the next year
% Saves JUST output statistics
% Code runs in parallel with parfor loop


%% Setup

close all; clear;

% ROMS data
FILE_NAME = 'add file path here'; %% Add path to ROMS file names
nc_file_names = dir(FILE_NAME);  % Directory of ROMS files - will list in order of year
years = 1990:2010;
mig_strat = 'personal_social'; % Options: personal_social, personal, social

% LHS parameters 
Nsamples = 1001; % number of LHC samples

% Intervals for parameters values
b1_bounds = [3,10]; b2_bounds = [0.3,0.4]; % Calling behaviors
c1_bounds = [10,30]; c2_bounds = [0,0.3];  % personal southward transition
c3_bounds = [0.1,0.3]; c4_bounds = [-20,0];% social southward transition

addpath Util/  % Add path to the functions this code calls

define_system_parameters;  % Script to define default parameters

%% LHS set-up
% Intervals
c1_interval = linspace(c1_bounds(1),c1_bounds(2),Nsamples);
c2_interval = linspace(c2_bounds(1),c2_bounds(2),Nsamples);
c3_interval = linspace(c3_bounds(1),c3_bounds(2),Nsamples);
c4_interval = linspace(c4_bounds(1),c4_bounds(2),Nsamples);

b1_interval = linspace(b1_bounds(1),b1_bounds(2),Nsamples);
b2_interval = linspace(b2_bounds(1),b2_bounds(2),Nsamples);

% Order to sample intervals
c1_order = datasample(1:Nsamples-1,Nsamples-1,'Replace',false);
c2_order = datasample(1:Nsamples-1,Nsamples-1,'Replace',false);
c3_order = datasample(1:Nsamples-1,Nsamples-1,'Replace',false);
c4_order = datasample(1:Nsamples-1,Nsamples-1,'Replace',false);

b1_order = datasample(1:Nsamples-1,Nsamples-1,'Replace',false);
b2_order = datasample(1:Nsamples-1,Nsamples-1,'Replace',false);


%% Run Simulations for each year and each parameter value

for k = 1:length(nc_file_names) % Loop over years
    
    disp(['Year: ' num2str(years(k))])
    
    myYear = years(k);
    
    par.start_date  = datetime(years(k),str2num(start_month_val),1);
   
    ref_date        = datetime(years(k),01,1);
    par.doy_start   = datenum(par.start_date) - datenum(ref_date);  % Day of the year the simulation starts on
   
    nc_file = [nc_file_names(k).folder '/' nc_file_names(k).name];
    
    % Read data from ncdf files
    sst_data   = ncread(nc_file,'temp'); sst_data = permute(sst_data,[2,1,4,3]);
    krill_data = ncread(nc_file,'Pzooplankton'); krill_data = permute(krill_data,[2,1,4,3]);

    % Setting up resolution information - corresponds the the input data
    % from the ncdf files
    grid_pars.resolution = 3000; 
    grid_pars.xrange  = grid_pars.resolution.*[0;size(sst_data,1)];
    grid_pars.yrange  = grid_pars.resolution.*[0;size(sst_data,2)];
    grid_pars.numX    = size(sst_data,1);
    grid_pars.numY    = size(sst_data,2);
    
    parfor (m = 1:Nsamples-1,4) % Loop over the parameter values
        
        disp(m);

        % Select values from within pre-defined intervals
        state3_pref = [0, 0, 10, 20, 0, 0, 10];
        state3_pref(1) = c1_interval(c1_order(m)) + (c1_interval(c1_order(m)+1) - c1_interval(c1_order(m))).*rand(1);
        state3_pref(2) = c2_interval(c2_order(m)) + (c2_interval(c2_order(m)+1) - c2_interval(c2_order(m))).*rand(1);
        state3_pref(5) = c3_interval(c3_order(m)) + (c3_interval(c3_order(m)+1) - c3_interval(c3_order(m))).*rand(1);
        state3_pref(6) = c4_interval(c4_order(m)) + (c4_interval(c4_order(m)+1) - c4_interval(c4_order(m))).*rand(1);
        
        call_signal_sent = [1;0;0];
        call_signal_sent(2) = b1_interval(b1_order(m)) + (b1_interval(b1_order(m)+1) - b1_interval(b1_order(m))).*rand(1);
        call_signal_sent(3) = b2_interval(b2_order(m)) + (b2_interval(b2_order(m)+1)  -b2_interval(b2_order(m))).*rand(1);
        
        % Handle parameter structure and parallel loop
        myPar = par;
        myPar.state3_pref = state3_pref;
        myPar.call_signal_sent = call_signal_sent;

       [~, Y, sst, krill,s , meanSumCalls, ~] = blueWhale_IBM(myPar,sst_data,krill_data,grid_pars,mig_strat); % Run IBM with given parameter values
    
    
        % Quick post-processing: remove the last location of the whales that went out of bounds because this location is out of bounds 
        for j = 1:par.numWhales 
         if nnz(Y(j,:)) < size(Y,2)  % If there is at least one zero
             Y(j,nnz(Y(j,:))) = 0;
         end 
        end
   
       % Compute statistics to save
       [full_pop_stats,mb_stats,above36_stats,mean_call_sum] = compute_migration_stats(s,Y,myPar);
      
      
       % Save data
       outpath_tmp = [mig_strat '_year_' num2str(myYear) '_LHS_' num2str(m) '.mat'];
       parsave_stats(outpath_tmp,myPar,full_pop_stats,mb_stats,above36_stats)
    
       
    end
    
end










