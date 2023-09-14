% Run blue whale IBM simulation for specified year(s)
% Null model: Gaussian day of year
% Stephanie Dodson: sdodson@colby.edu
% September 2023

%% Setup

close all; clear;


% ROMS data
FILE_NAME = 'add file path here';
nc_file_names = dir(FILE_NAME);  % Directory of ROMS files - will list in order of year
years = 1990:2010;

% Run Parameters
mig_strat = 'doy'; % Options: 'no_mig' (no migration) & 'doy' (normally distributed migration dates based on yearday)
numRuns = 1; % Number of runs/year
parallel_flag = 0;  % 0 - serial, 1 - parallel

addpath Util/  % Add path to the functions this code calls

% Define default system parameteres
system_parameters_doy; 
%% Run Simulations for each year

for k = 1:length(nc_file_names) 
    
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

    if parallel_flag % Run code in parallel
        parfor (m = 1:numRuns,8)
            disp(m)
        
        [~, Y, sst, krill, s, forage_eff] = blueWhale_IBM_nullModels(par,sst_data,krill_data,grid_pars,mig_strat);
        
        % Quick post-processing: remove the last location of the whales that went out of bounds because this location is out of bounds 
        for j = 1:par.numWhales
           
            if nnz(Y(j,:)) < size(Y,2)  % If there is at least one zero
                Y(j,nnz(Y(j,:))) = 0;
            end
            
        end
        
        outfile = [mig_strat '_simulation_' num2str(myYear) '_' num2str(m) '.mat'];
        parsave_slim(outfile,Y,krill,s,par,grid_pars,forage_eff); % Save only most crucial information (for space saving)

        
        end

    else % Serial code
        for m = 1:numRuns

            [~, Y, sst, krill, s, forage_eff] = blueWhale_IBM_nullModels(par,sst_data,krill_data,grid_pars,mig_strat);

            % Quick post-processing: remove the last location of the whales that went out of bounds because this location is out of bounds 
            for j = 1:par.numWhales
               
                if nnz(Y(j,:)) < size(Y,2)  % If there is at least one zero
                    Y(j,nnz(Y(j,:))) = 0;
                end
                
            end
            
            outfile = [mig_strat '_simulation_' num2str(myYear) '_' num2str(m) '.mat'];
            save(outfile,'Y','s','krill','par','grid_pars','forage_eff','-v7.3'); % Save a Matlab file with only most crucial information (for space saving)


        end
    end % End if parallel flag
    
end















