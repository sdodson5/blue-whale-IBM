function [full_pop_stats,mb_stats,above36_stats,mean_call_sum] = compute_migration_stats(s,Y,par,meanSumCalls)

% Output gives: [median migration date, number whales migrating, IQR(1),
% IQR(2)]
% for full population, whales that migrate from MB, and whales that
% migrate from lats 36+

full_pop_stats = zeros(4,1);
mb_stats = zeros(4,1);
above36_stats = zeros(4,1);

%% Find full migration distribution
    s3 = zeros(size(s)); s3(s==3) = 1;
    
    numWhales = sum(s>=1); % Number of whales in the entire domain (full population)

    % Find first date transition to s3
    s3_transDate = zeros(size(s3,1),1);
    for j = 1:size(s3,1)

        tmp = find(s3(j,:)==1,1);

        if ~isempty(tmp)
            s3_transDate(j) = tmp;     
        end

    end

   %% Group the migration times based on the approximate migration latitudes

    rad_earth = 6378.137;  % Radius of the earth (KM)
    MM = 360/(2*pi*rad_earth*1000); % Scaling factor
    lat_ref = 32;

    % Map coordinates to lat-long
    Y_lat  = lat_ref + Y.*MM;

    s3_idx = find(s3_transDate > 0);  % Find indices of whales that did transition to southward migration.
    transLat = zeros(size(Y,1),1);
    
    % Find the latitudes at time of migration
    for jj = 1:length(s3_idx)

        idx = s3_idx(jj);
        transLat(idx) = Y_lat(idx,s3_transDate(idx));

    end

    % Remove the whales that did not migrate
    transLat = transLat(s3_idx); % Latitudes the whales migrated from

    tmp = ceil(s3_transDate/4) + par.doy_start;
    doyMig = tmp(s3_idx);  % Day of year migrated: full population

    [~,~,latBins] = histcounts(transLat,32:42);

    % Find the distributions for 36+ lat
    monterey_mig = doyMig(latBins==5);          % Lats 36-67
    mig_above_36lat = doyMig(latBins>=5); % Lats 36+
  
 %% Stats
    % Full population stats
    full_pop_stats(1) = median(doyMig);
    full_pop_stats(2) = length(doyMig);
    full_pop_stats(3:4) = quantile(doyMig,[0.25,0.75]);
    
    
    % MB stats
    mb_stats(1) = median(monterey_mig);
    mb_stats(2) = length(monterey_mig);
    mb_stats(3:4) = quantile(monterey_mig,[0.25,0.75]);
    
    % Lats 36+ stats
    above36_stats(1) = median(mig_above_36lat);
    above36_stats(2) = length(mig_above_36lat);
    above36_stats(3:4) = quantile(mig_above_36lat,[0.25,0.75]);    
    
    
    %% Mean call sum (optional)
    
if nargin == 4
% Mean call sum
    meanSumCalls = meanSumCalls(s3_idx,:);
    
    mean_call_sum = zeros(3,size(s,2));
    mean_call_sum(1,:) = sum(meanSumCalls)./numWhales; % Full population
    
    tmp = s(latBins==5,:);
    tmp(tmp>=1) = 1;
    numMB = sum(tmp);
    mean_call_sum(2,:) = sum(meanSumCalls(latBins==5,:))./numMB;
    
    tmp = s(latBins>=5,:);
    tmp(tmp>=1) = 1;
    num36p = sum(tmp);
    mean_call_sum(3,:) = sum(meanSumCalls(latBins>=5,:))./num36p;
    
else
    mean_call_sum = 0;
    
end

