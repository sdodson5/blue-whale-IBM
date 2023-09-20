% Plot boxplot of migration dates for low, medium, and high krill years 
% Plot boxplot of krill intake - compared to non-migratory null model
% intake
%% Set up
close all; clear;

extra_models = 0; % 0 produces Figure 3 in the main text, 1 produces Figure 7 in supplementary text

% Plotting colors
FE_color = "#33BBEE";
socialFE_color = "#009988"; 
FEminKrill_color = "#0077BB";
social_color = "#EE7733"; 

doy_color = "#CC3311";
noMig_color = [0.64,0.08,0.18];

%% Load Data

% Artifact of how I initially saved data and strategies together

% Strategy order: doy, FE (personal), FE+social (personal & social), FE+min
% krill (personal & min Krill)
low_mig = readmatrix('migration_dates_lowKrill.csv'); 
med_mig = readmatrix('migration_dates_medKrill.csv'); 
high_mig= readmatrix('migration_dates_highKrill.csv'); 

% Social Only strategy - summary data is in its own file
low_mig_social = readmatrix('migration_dates_socialOnly_low.csv');
med_mig_social = readmatrix('migration_dates_socialOnly_avg.csv');
high_mig_social = readmatrix('migration_dates_socialOnly_high.csv');


% Krill
low_krill = readmatrix('krillIntake_low.csv'); 
med_krill = readmatrix('krillIntake_avg.csv'); 
high_krill = readmatrix('krillIntake_high.csv'); 

% Social-only krill
low_krill_social = readmatrix('krillIntake_socialOnly_low.csv');
med_krill_social = readmatrix('krillIntake_socialOnly_avg.csv');
high_krill_social = readmatrix('krillIntake_socialOnly_high.csv');


if extra_models 
    low_mig = [low_mig(1:2,:);low_mig(4,:); low_mig(3,:)];
    med_mig = [med_mig(1:2,:);med_mig(4,:); med_mig(3,:)];
    high_mig = [high_mig(1:2,:);high_mig(4,:); high_mig(3,:)];

    myColors = [doy_color, FE_color, FEminKrill_color, socialFE_color, social_color];

    % Original order = doy, FE, FE+social, FE+min Krill, noMig. Add social only
    low_krill = [low_krill(1:2,:); low_krill(4,:);low_krill(3,:); low_krill_social; low_krill(end,:)];
    med_krill = [med_krill(1:2,:); med_krill(4,:);med_krill(3,:); med_krill_social; med_krill(end,:)];
    high_krill = [high_krill(1:2,:); high_krill(4,:); high_krill(3,:); high_krill_social; high_krill(end,:)];

    % Append social- only strategy to matrices
    % New order: doy, FE, FE+social, Social Only (better for plotting)
    low_mig = [low_mig; low_mig_social];
    med_mig = [med_mig; med_mig_social];
    high_mig = [high_mig; high_mig_social];

else % need to remove the doy &  FE + min krill & social 
    low_mig = low_mig(2:end-1,:);
    med_mig = med_mig(2:end-1,:);
    high_mig = high_mig(2:end-1,:);

    % Removing the doy & FE + min Krill & social only: order = doy, FE, FE+social, social only,
% nonMig
    low_krill = [low_krill(2:end-2,:); low_krill(end,:)];
    med_krill = [med_krill(2:end-2,:);  med_krill(end,:)];
    high_krill = [high_krill(2:end-2,:); high_krill(end,:)];

    myColors = [FE_color, socialFE_color];
end



str_cell= {'Low'; 'Average'; 'High'};

%% Map krill intake to proportion of the non-migratory population

low_krill_dev = 100 + 100*(low_krill - low_krill(end,2))./low_krill(end,2); 
med_krill_dev = 100 + 100*(med_krill - med_krill(end,2))./med_krill(end,2); 
high_krill_dev = 100 + 100*(high_krill - high_krill(end,2))./high_krill(end,2); 

figure; subplot(2,1,1); hold on; 

ff = 0;
subplot(2,1,2);
xlim([-4,30])
plot(xlim,[100,100],'k-','linewidth',1)

 % Low years
 for k = 1:size(low_mig,1)

    subplot(2,1,1);
    plot(ff,low_mig(k,2),'s','MarkerSize',12,'MarkerEdgeColor', myColors(k),'MarkerFaceColor',myColors(k));
    plot([ff,ff],[low_mig(k,1),low_mig(k,3)],'linewidth',3,'Color',myColors(k) );

    subplot(2,1,2); hold on;
    plot(ff,low_krill_dev(k,2),'s','MarkerSize',12,'MarkerEdgeColor', myColors(k),'MarkerFaceColor',myColors(k));
    plot([ff,ff],[low_krill_dev(k,1),low_krill_dev(k,3)],'linewidth',3,'Color',myColors(k) );
    
    ff = ff+0.25;
 end


 ff = ff+1;

 % Average years
 for k = 1:size(med_mig,1)

    subplot(2,1,1);
    plot(ff,med_mig(k,2),'s','MarkerSize',12,'MarkerEdgeColor', myColors(k),'MarkerFaceColor',myColors(k));
    plot([ff,ff],[med_mig(k,1),med_mig(k,3)],'linewidth',3,'Color',myColors(k) );

    subplot(2,1,2);
    plot(ff,med_krill_dev(k,2),'s','MarkerSize',12,'MarkerEdgeColor', myColors(k),'MarkerFaceColor',myColors(k));
    plot([ff,ff],[med_krill_dev(k,1),med_krill_dev(k,3)],'linewidth',3,'Color',myColors(k) );
    ff = ff+0.25;
 end


  ff = ff+1;

 % High years
 for k = 1:size(high_mig,1)

    subplot(2,1,1);
    plot(ff,high_mig(k,2),'s','MarkerSize',12,'MarkerEdgeColor', myColors(k),'MarkerFaceColor',myColors(k));
    plot([ff,ff],[high_mig(k,1),high_mig(k,3)],'linewidth',3,'Color',myColors(k) );

    subplot(2,1,2);
    plot(ff,high_krill_dev(k,2),'s','MarkerSize',12,'MarkerEdgeColor', myColors(k),'MarkerFaceColor',myColors(k));
    plot([ff,ff],[high_krill_dev(k,1),high_krill_dev(k,3)],'linewidth',3,'Color',myColors(k) );
    ff = ff+0.25;
 end

 % Need to adjust x-ticks and x-limits if plotting all strategies
 % (supplementary figure 7).
subplot(2,1,1);
set(gca,'fontsize',16,'XTick',[]);
ylabel('Yearday'); xlabel('');
xlim([-0.5,3.75])
ylim([200,330])
set(gca,'linewidth',2);
xticks([0.125 1.625 3.125])
%xticks([2, 12, 22]);
xticklabels('')

subplot(2,1,2);
set(gca,'fontsize',16,'XTick',[]);
ylabel('Relative Krill Intake (%)'); xlabel('');
xticks([0.125 1.625 3.125])
%xticks([2, 12, 22]);
xlim([-0.5,3.75])
ylim([20,120])
set(gca,'linewidth',2);
xticklabels(str_cell);
box off;


