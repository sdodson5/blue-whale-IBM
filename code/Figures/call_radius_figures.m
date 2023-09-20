%% Analysis of the impact of call radius: 
%
% (1) Plot migration dates as a function of the call radius: aggregate data  over all years
% (2) Plot krill intake as a function of call radius: aggregate data over all years


close all; clear;

% File names: call radii
data = cell(10,1);
data{1} = 'summary_data/FE_only_summary/migDates_krillIntake_summary_FEOnly_*.mat'; % No calls (call radius = 0)
data{2} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_05km_*.mat'; % 5km call distance
data{3} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_10km_*.mat'; % 10km call distance
data{4} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_25km_*.mat'; % 25km call distance
data{5} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_50km_*.mat'; % 50km call distance
data{6} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_100km_*.mat'; % 100km call distance
data{7} = 'summary_data/socialFE_summary/migDates_krillIntake_summary_socialFE_*.mat'; % 125km call distance (default) 
data{8} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_200km_*.mat'; % 200km call distance
data{9} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_300km_*.mat'; % 300km call distance
data{10} = 'summary_data/socialFE_callDist_summary/migDates_krillIntake_summary_socialFE_500km_*.mat'; % 500 km call distance

call_dist = [0,5,10,25,50,100,125,200,300,500];
myXTicks = {'0','5','10','25','50','100','125','200','300','500'};
totYears = 1990:2010;

% Non-migratory population krill intake
nonMigKrill = readmatrix('paper_figures/nonMigratory_totalKrillIntake_year_avg_median.csv');
nonMigKrill = nonMigKrill(:,[1,3]);

IQR = zeros(length(data),length(totYears));
migStats = zeros(length(data),5);


myColors = ["#33BBEE","#90d3cb", "#7bcac2","#67c2b8","#4db8ac","#1aa394", "#009988","#008a7a", "#007a6d","#006b5f"];

personal_krill = 0;
my_medians = zeros(length(data),2); % Store median migration dates and krill intake (in order to put a thin line through the data points). 

for m = 1:length(data) % Loop over call distances

    migDates = [];
    krillIntake = [];
    krillIntake_raw= [];


    files = dir(data{m});

    for k = 1:length(files) % Loop over years
        
        load([files(k).folder '/' files(k).name]);

        migDates = [migDates;mig_hist.FP]; % Aggregate all migration information
        krillIntake = [krillIntake; (totalKrill-nonMigKrill(k,2))./nonMigKrill(k,2)];
        krillIntake_raw = [krillIntake_raw; totalKrill];


       IQR(m,k) = median_migDates(1,3) - median_migDates(1,1);

    end


    % migration statistics: all , low, average, high krill years
    q  = quantile(migDates,[0.25,0.5,0.75]);
   
    migStats(m,:) = [call_dist(m), q, length(migDates)];

    % krill intake statistics: all , low, average, high krill years
    KI  = quantile(krillIntake,[0.25,0.5,0.75]);
    KI_raw = quantile(krillIntake_raw,[0.25,0.5,0.75]);

    if m == 1
        personal_krill = median(krillIntake);
        krill_diff = (KI - personal_krill)./abs(personal_krill);
    else

        krill_diff = (KI - personal_krill)./abs(personal_krill);

    end

    % Plot migration data - all years
    figure(1);
    subplot(2,1,1); hold on;
    plot([call_dist(m),call_dist(m)],[q(1),q(3)],'-','LineWidth',3,'Color',myColors(m));
    plot(call_dist(m), q(2),'s','MarkerSize',10,'MarkerFaceColor',myColors(m),'MarkerEdgeColor',myColors(m))

    % Social - personal krill
    %figure(2); hold on;
    subplot(2,1,2); hold on;
    plot([call_dist(m),call_dist(m)],[krill_diff(1),krill_diff(3)],'-','LineWidth',3,'Color',myColors(m));
    plot(call_dist(m), krill_diff(2),'s','MarkerSize',10,'MarkerFaceColor',myColors(m),'MarkerEdgeColor',myColors(m))

    my_medians(m,:) = [q(2), krill_diff(2)];


end


% Figure labels
subplot(2,1,1);
plot(call_dist, my_medians(:,1),'-','linewidth',1.5,'Color',[0.4,0.4,0.4])
for m = 1:7 %length(data)
     plot(call_dist(m), my_medians(m,1),'s','MarkerSize',10,'MarkerFaceColor',myColors(m),'MarkerEdgeColor',myColors(m))
  
end
ylabel('Migration Onset (yearday)');
set(gca,'fontsize',16,'linewidth',2);
%xlim([-20,520])
ylim([200,320]);
%xticks(1:length(data));
%xticklabels(myXTicks)
xlim([-2,call_dist(end)+2]);
xticks(call_dist)
xticklabels(myXTicks)


subplot(2,1,2);
plot(call_dist, my_medians(:,2),'-','linewidth',1.5,'Color',[0.4,0.4,0.4])
for m = 1:7%length(data)
     plot(call_dist(m), my_medians(m,2),'s','MarkerSize',10,'MarkerFaceColor',myColors(m),'MarkerEdgeColor',myColors(m))
  
end
xlabel('Max. Call Radius (km)');
ylabel('Krill Intake Ratio');
set(gca,'fontsize',16,'linewidth',2);
%xlim([-20,520])
ylim([-1,1.5]);
xlim([-2,call_dist(end)+2]);
xticks(call_dist)
xticklabels(myXTicks)












