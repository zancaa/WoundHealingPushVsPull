%% Script to reproduce figures for Zanca et al. 2022, 'Push or pull? Cell proliferation and migration during wound healing'
% Note: You must be in the WoundHealing folder created to store outputs
% from the simulations

%% Model 1: Continual growth of proliferative hub
% Initialise variables
seeds = 1:10;
forces = 0:10;
qvf = 0.6:0.05:1.2;
dt = 0.0005;
t = 0:0.5:48;
tissue_length_continual = zeros(length(t),length(forces),length(qvf),length(seeds));
tissue_length_free = zeros(length(t),length(forces));
end_times_continual = zeros(length(forces),length(qvf),length(seeds));
tissue_length_ave = zeros(length(t),length(forces),length(qvf));
grads_continual = zeros(length(forces),length(qvf));
grads_free = zeros(length(forces));
tissue_vel_continual = zeros(length(t),length(forces),length(qvf),length(seeds));
tissue_vel_ave = zeros(length(t),length(forces),length(qvf));
tissue_vel_free = zeros(length(forces),1);

for i4 = 1:length(forces)
    for i3 = 1:length(qvf)       
        for i2 = 1:length(seeds)
            % Begin by loading data
%             dir_path = strcat('LeadingEdgeForce_',num2str(forces(i4)),'/',...
%                 'MeanQuiescentVolumeFraction_',num2str(qvf(i3)),'/',...
%                 'Run_',num2str(i2),'/results_from_time_0/');
            dir_path = strcat('2DSweep10by10Dt0.0005/AddMaintainFreeBoundaryForce_',num2str(forces(i4)),'/',...
            'MeanQuiescentVolumeFraction_',num2str(qvf(i3)),'/',...
            'Run_',num2str(i2),'/results_from_time_0/');
            % Note: boundarynodevelocities give displacement NOT velocity
            node_vel = importdata(strcat(dir_path,'boundarynodevelocities.dat'));
    %             divisions = importdata('divisions.dat');
    %             divisions = cell2mat(divisions);
    %             k = find(divisions(1,:) > 10,1);
    %             division_ages{i3,i1,i2} = divisions(4,k:end);
    %             division_heights{i3,i1,i2} = divisions(3,k:end);
            max_y = zeros(1,length(t));
            end_times_continual(i4,i3,i2) = t(length(node_vel));
            for i = 1:length(node_vel)
                ns = node_vel{i};
                ns = str2num(ns);
                nv = ns(6:5:end);
                ns = ns(4:5:end);
                [max_y(i),max_ind] = max(ns);
                tissue_length_continual(i,i4,i3,i2) = max_y(i);
                tissue_vel_continual(i,i4,i3,i2) = nv(1)/dt;
            end
            if length(node_vel) < length(t)
                tissue_length_continual(i+1:end,i4,i3,i2) = NaN(length(i+1:length(t)),1);
                tissue_vel_continual(i+1:end,i4,i3,i2) = NaN(length(i+1:length(t)),1);
            end
        end
    end
    
    tissue_vel_mean_last_12_continual = squeeze(mean(tissue_vel_continual(t>=36,:,:,:),1,'omitnan'));
    tissue_vel_overall_mean_continual = squeeze(mean(tissue_vel_mean_last_12_continual,3,'omitnan'));
    tissue_vel_mean_10_quantile_continual = squeeze(quantile(tissue_vel_mean_last_12_continual,0.1,3));
    tissue_vel_mean_90_quantile_continual = squeeze(quantile(tissue_vel_mean_last_12_continual,0.9,3));
    
    % Load free tissue data
%     dir_path = strcat(strcat('FreeTissue/',...
%     'LeadingEdgeForce_',num2str(forces(i4)),'/',...
%     'results_from_time_0/'));
    dir_path = strcat(strcat('NoAttach2/NoAttach/',...
    'AddMaintainFreeBoundaryForce_',num2str(forces(i4)),'/',...
    'MeanQuiescentVolumeFraction_0/IC10/dt_0.0005/Run_1/results_from_time_0/'));
    node_vel = importdata(strcat(dir_path,'boundarynodevelocities.dat'));
    max_y = zeros(1,length(t));
    for i = 1:length(node_vel)
        ns = node_vel{i};
        ns = str2num(ns);
        if i == length(node_vel)
            nv = ns(6:5:end);
        end
        ns = ns(4:5:end);
        [max_y(i),max_ind] = max(ns);
        max_y(i) = max_y(i);
        tissue_length_free(i,i4) = max_y(i);  
    end
    tissue_vel_free(i4) = nv(1)/dt;
end

%% Produce plots of tissue length and velocity and region heatmap
% Note that the simulation snapshots in Figures 2(a), 2(c), 4(b)-(d) were 
% produced in Paraview.

% Figure 2(b)
figure
hold on
i4 = 4;
i3 = 7;
for i2 = 1:length(seeds)
    p1 = plot(t,tissue_length_continual(:,i4,i3,i2),'Color',[0 0 1 0.4]);
end
p2 = plot(t,mean(tissue_length_continual(:,i4,i3,:),4,'omitnan'),'Color',[0 0 1],'LineWidth',1);
p3 = plot(t,tissue_length_free(:,i4),'k--','LineWidth',1);
set(gca,'FontName','Times New Roman')
xlabel('Time (hrs)','Interpreter','Latex')
ylabel('Tissue length (CDs)','Interpreter','Latex')
axis([0 48 0 70])
xticks([0 12 24 36 48])
ax = gca;
ax.XGrid = 'on';
legend([p1 p2 p3],{'Individual simulations','Average of all simulations','Free tissue'},...
    'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')

% Figure 2(d)
figure
hold on
i4 = 4;
i3 = 9;
for i2 = 1:length(seeds)
    p1 = plot(t,tissue_length_continual(:,i4,i3,i2),'Color',[1 0 0 0.4]);
end
p2 = plot(t,mean(tissue_length_continual(:,i4,i3,:),4,'omitnan'),'Color',[1 0 0],'LineWidth',1);
p3 = plot(t,tissue_length_free(:,i4),'k--','LineWidth',1);
set(gca,'FontName','Times New Roman')
xlabel('Time (hrs)','Interpreter','Latex')
ylabel('Tissue length (CDs)','Interpreter','Latex')
axis([0 48 0 70])
xticks([0 12 24 36 48])
ax = gca;
ax.XGrid = 'on';
legend([p1 p2 p3],{'Individual simulations','Average of all simulations','Free tissue'},...
    'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')

% Figure 2(e)
figure
hold on
for i = 1:length(0.6:0.05:1.2)
    p1 = plot(t,mean(tissue_length_continual(:,i4,i,:),4,'omitnan'),'Color',[0 0 0 0.4]);
    if i == 7
        p2 = plot(t,mean(tissue_length_continual(:,i4,i,:),4,'omitnan'),'Color',[0 0 1],'LineWidth',1.5);
    end
    if i == 9
        p3 = plot(t,mean(tissue_length_continual(:,i4,i,:),4,'omitnan'),'Color',[1 0 0],'LineWidth',1.5);
    end
    if i > 11
        p4 = plot(t,mean(tissue_length_continual(:,i4,i,:),4,'omitnan'),'Color',[0.15 0.4 0 0.6]);
    end
end
p5 = plot(t,tissue_length_free(:,i4),'k--','LineWidth',1.5);
set(gca,'FontName','Times New Roman')
xlabel('Time (hrs)','Interpreter','Latex')
ylabel('Tissue length (CDs)','Interpreter','Latex')
legend([p1 p2 p3 p4 p5],{'Average tissue length','$\mu_\mathrm{QVF} = 0.9$','$\mu_\mathrm{QVF} = 1.0$','No growth','Free tissue'},'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')
axis([0 48 0 70])
xticks([0 12 24 36 48])
ax = gca;
ax.XGrid = 'on';

% Figure 3(a)
figure
hold on
i3 = 7;
for i2 = 1:10
    p1 = plot(t,tissue_vel_continual(:,i4,i3,i2),'Color',[0 0 1 0.4]);
end
axis([0 48 0 2])
p2 = plot(t,tissue_vel_free(i4)*ones(length(t),1),'k--');
fill([36 48 48 36 36],[0 0 2 2 0],[1 1 1],'FaceAlpha',0)
for i2 = 1:10
    p3 = plot(t(t>=36),mean(tissue_vel_continual(t>=36,i4,i3,i2),1)*ones(length(t(t>=36)),1),'k-');
end
set(gca,'FontName','Times New Roman')
legend([p1 p2 p3],{'Instantaneous','Free tissue','Average over last 12 hours'},'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')
xlabel('Time (hrs)','Interpreter','Latex')
ylabel('Wound edge velocity (CDs/hr)','Interpreter','Latex')
xticks([0 12 24 36 48])

% Figure 3(b)
figure
hold on
i3 = 9;
for i2 = 1:10
    p1 = plot(t,tissue_vel_continual(:,i4,i3,i2),'Color',[1 0 0 0.4]);
end
axis([0 48 0 2])
p2 = plot(t,tissue_vel_free(i4)*ones(length(t),1),'k--');
fill([36 48 48 36 36],[0 0 2 2 0],[1 1 1],'FaceAlpha',0)
for i2 = 1:10
    p3 = plot(t(t>=36),mean(tissue_vel_continual(t>=36,i4,i3,i2),1)*ones(length(t(t>=36)),1),'k-');
end
set(gca,'FontName','Times New Roman')
legend([p1 p2 p3],{'Instantaneous','Free tissue','Average over last 12 hours'},'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')
xlabel('Time (hrs)','Interpreter','Latex')
ylabel('Wound edge velocity (CDs/hr)','Interpreter','Latex')
xticks([0 12 24 36 48])

% Figure 3(c)
figure 
hold on
wound_edge_vel_ub = max(tissue_vel_overall_mean_continual(i4,:)) + 0.1*max(tissue_vel_overall_mean_continual(i4,:));
prolif_enhanced_ub = qvf(find(tissue_vel_overall_mean_continual(i4,:)<=tissue_vel_free(i4),1));
fill([min(qvf),prolif_enhanced_ub,prolif_enhanced_ub,min(qvf),min(qvf)],...
    [0,0,wound_edge_vel_ub,wound_edge_vel_ub,0],[0 0 1],'FaceAlpha',0.3,'EdgeColor','None')
prolif_inhib_ub = qvf(find(tissue_vel_overall_mean_continual(i4,:)<=0.0025,1));
fill([prolif_enhanced_ub,prolif_inhib_ub,prolif_inhib_ub,prolif_enhanced_ub,prolif_enhanced_ub],...
    [0,0,wound_edge_vel_ub,wound_edge_vel_ub,0],[1 0 0],'FaceAlpha',0.3,'EdgeColor','None')
fill([prolif_inhib_ub,max(qvf),max(qvf),prolif_inhib_ub,prolif_inhib_ub],...
    [0,0,wound_edge_vel_ub,wound_edge_vel_ub,0],[0 0 0],'FaceAlpha',0.3,'EdgeColor','None')
p1 = errorbar(qvf,tissue_vel_overall_mean_continual(i4,:),...
    tissue_vel_overall_mean_continual(i4,:)-tissue_vel_mean_10_quantile_continual(i4,:),...
    tissue_vel_mean_90_quantile_continual(i4,:)-tissue_vel_overall_mean_continual(i4,:),...
    'Color',[0 0 0],'CapSize',3);
p2 = plot(qvf,tissue_vel_free(i4)*ones(length(qvf),1),'k--');
plot(0.9,tissue_vel_overall_mean_continual(i4,7),'bo','MarkerFaceColor','b')
plot(1.0,tissue_vel_overall_mean_continual(i4,9),'ro','MarkerFaceColor','r')
text(0.746,0.8726,{'Proliferation','enhanced'},'Interpreter','Latex','HorizontalAlignment','center','Rotation',90)
text(1.05,0.8726,{'Proliferation','inhibited'},'Interpreter','Latex','HorizontalAlignment','center','Rotation',90)
text(1.13,0.8726,{'No growth'},'Interpreter','Latex','HorizontalAlignment','center','Rotation',90)
set(gca,'FontName','Times New Roman')
axis([min(qvf) max(qvf) 0 wound_edge_vel_ub])
legend([p1 p2],{'Average','Free tissue'},'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')
xlabel('$\mu_\mathrm{QVF}$','Interpreter','Latex')
ylabel('Wound edge velocity (CDs/hr)','Interpreter','Latex')

% Figure 4(a)
grads_diff_continual = tissue_vel_overall_mean_continual - tissue_vel_free;
figure
hold on
surf((qvf),(forces),grads_diff_continual,'EdgeColor','interp','FaceColor','interp')
set(gca,'FontName','Times New Roman')
xlabel('$\mu_\mathrm{QVF}$','Interpreter','Latex')
ylabel('$F_\mathrm{active}$','Interpreter','Latex')
c = colorbar;
c.Label.String = 'Average wound edge velocity - free tissue velocity';
[gradblu,~]=colorGradient([1 1 1],[0 0 1],72);
[gradred,im]=colorGradient([1 0 0],[1 1 1],72);
cmap = [gradred; gradblu(1:48,:)];
colormap(cmap)
axis([min(qvf) max(qvf) min(forces) max(forces)])
% Add proliferation enhanced/inhibited boundary (needs to be plotted in 3D
% so the contour appears above the surface)
M2 = contour(qvf,forces,grads_diff_continual,[0 0],'EdgeColor','None','LineWidth',2);
plot3(M2(1,2:20),M2(2,2:20),(max(max(grads_diff_continual))+0.1)*ones(length(M2(1,2:20)),1),'k-','LineWidth',2)
% Add no growth patch
M = contour(qvf,forces,tissue_vel_overall_mean_continual,[0.0025 0.0025],'Color','None');
tmp = smooth(M(1,2:end));
tmp2 = smooth(M(2,2:end));
patch([tmp; 1.2; 1],[tmp2; 0; 0],(max(max(grads_diff_continual))+0.1)*ones(length(tmp)+2,1),[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],'LineWidth',2)
% Add markers for simulation snapshots
plot3(0.9,3,(max(max(grads_diff_continual))+0.12),'ko')
plot3(1,3,(max(max(grads_diff_continual))+0.12),'ks')
plot3(1.15,3,(max(max(grads_diff_continual))+0.12),'kx')
% Add text to the different regions
text(0.65,1,(max(max(grads_diff_continual))+0.12),{'Proliferation','enhanced'},...
    'HorizontalAlignment','center','Color','w','Interpreter','Latex')
text(0.85,9,(max(max(grads_diff_continual))+0.12),{'Proliferation','inhibited'},...
    'HorizontalAlignment','center','Interpreter','Latex')
text(1.1,1,(max(max(grads_diff_continual))+0.12),{'No','growth'},...
    'HorizontalAlignment','center','Interpreter','Latex')

%% Model 2: Bounded growth of proliferative hub

% Initialise variables
forces = 0:10;
qvf = 0.6:0.05:1;
max_hub_size = 1:21;
tissue_length_bounded = zeros(length(t),length(max_hub_size),length(forces),length(qvf),length(seeds));
tissue_length_ave_bounded = zeros(length(t),length(max_hub_size),length(forces),length(qvf));
end_times_bounded = zeros(length(max_hub_size),length(forces),length(qvf),length(seeds));
grads_bounded = zeros(length(max_hub_size),length(forces),length(qvf));
grads_diff_bounded = zeros(length(max_hub_size),length(forces),length(qvf));
tissue_vel_bounded = zeros(length(t),length(max_hub_size),length(forces),length(qvf),length(seeds));

% Load data

for i5 = 1:length(max_hub_size)
    for i4 = 1:length(forces)
        for i3 = 1:length(qvf)   
            for i2 = 1:length(seeds)
                % Begin by loading data
%                 dir_path = strcat('dmax_',num2str(max_hub_size(i5)),'/',...
%                     'LeadingEdgeForce_',num2str(forces(i4)),'/',...
%                     'MeanQuiescentVolumeFraction_',num2str(qvf(i3)),'/',...                    
%                     'Run_',num2str(i2),'/results_from_time_0/');
                dir_path = strcat('Window/',...
                    'AddMaintainFreeBoundaryForce_',num2str(forces(i4)),'/',...
                    'MeanQuiescentVolumeFraction_',num2str(qvf(i3)),'/',...
                    'BottomOfProlifHub_',num2str(max_hub_size(i5)),'/',...
                    'Run_',num2str(i2),'/results_from_time_0/');
                % Note: the velocities are actually a displacement, NOT
                % velocity
                node_vel = importdata(strcat(dir_path,'boundarynodevelocities.dat'));
                max_y = zeros(1,length(t));
                % Calculate and save the cell and tissue areas and number of
                % cells at each time t.
                end_times_bounded(i5,i4,i3,i2) = t(length(node_vel));
                for i = 1:length(node_vel)
                    ns = node_vel{i};
                    ns = str2num(ns);
                    nv = ns(6:5:end);
                    ns = ns(4:5:end);
                    [max_y(i),max_ind] = max(ns);
                    tissue_length_bounded(i,i5,i4,i3,i2) = max_y(i);
                    tissue_vel_bounded(i,i5,i4,i3,i2) = nv(1)/dt;
                end
                if length(node_vel) < length(t)
                    tissue_length_bounded(i+1:end,i5,i4,i3,i2) = NaN(length(i+1:length(t)),1);
                    tissue_vel_bounded(i+1:end,i5,i4,i3,i2) = NaN(length(i+1:length(t)),1);
                end
            end
        end
    end
end

tissue_vel_mean_last_12_bounded = squeeze(mean(tissue_vel_bounded(t>=36,:,:,:,:),1,'omitnan'));
tissue_vel_overall_mean_bounded = squeeze(mean(tissue_vel_mean_last_12_bounded,4,'omitnan'));
tissue_vel_mean_10_quantile_bounded = squeeze(quantile(tissue_vel_mean_last_12_bounded,0.1,4));
tissue_vel_mean_90_quantile_bounded = squeeze(quantile(tissue_vel_mean_last_12_bounded,0.9,4));

for i5 = 1:length(max_hub_size)
    grads_diff_bounded(i5,:,:) = squeeze(tissue_vel_overall_mean_bounded(i5,:,:)) - tissue_vel_free;
end

%% Produce plots of tissue length and velocity
% Note: Simulation plots in Figure 5(a) were produced in Paraview
% Figure 5(b)
figure
hold on
i5 = 5;
i4 = 4;
i3 = 7;
for i2 = 1:length(seeds)
    p3 = plot(t,tissue_length_bounded(:,i5,i4,i3,i2),'Color',[1 0 0 0.4]);
end
p4 = plot(t,mean(tissue_length_bounded(:,i5,i4,i3,:),5,'omitnan'),'Color',[1 0 0],'LineWidth',1);
p1 = plot(t,mean(tissue_length_continual(:,i4,i3,:),4,'omitnan'),'Color',[0 0 1],'LineWidth',1);
p2 = plot(t,tissue_length_free(:,i4),'k--','LineWidth',1);
set(gca,'FontName','Times New Roman')
xlabel('Time (hrs)','Interpreter','Latex')
ylabel('Tissue length (CDs)','Interpreter','Latex')
axis([0 48 0 70])
xticks([0 12 24 36 48])
ax = gca;
ax.XGrid = 'on';
legend([p1 p2 p3 p4],{'Continual growth model','Free tissue',...
    'Individual runs bounded growth model','Average of all runs bounded growth model'},...
    'Location','NorthWest','EdgeColor','None','Color','None','Interpreter','Latex')

% Figure 5(c)
figure
hold on
p1 = fill([min(max_hub_size) max(max_hub_size) max(max_hub_size) min(max_hub_size) min(max_hub_size)],...
    [tissue_vel_mean_10_quantile_continual(i4,i3) tissue_vel_mean_10_quantile_continual(i4,i3) ...
    tissue_vel_mean_90_quantile_continual(i4,i3) tissue_vel_mean_90_quantile_continual(i4,i3) ...
    tissue_vel_mean_10_quantile_continual(i4,i3)],'b','FaceAlpha',0.4,'EdgeColor','None');
p2 = plot(max_hub_size,tissue_vel_overall_mean_continual(i4,i3)*ones(length(max_hub_size),1),'b-');
p3 = plot(max_hub_size,tissue_vel_free(i4)*ones(length(max_hub_size),1),'k--');
p4 = errorbar(max_hub_size,squeeze(tissue_vel_overall_mean_bounded(:,i4,i3)),...
    squeeze(tissue_vel_overall_mean_bounded(:,i4,i3))-squeeze(tissue_vel_mean_10_quantile_bounded(:,i4,i3)),...
    squeeze(tissue_vel_mean_90_quantile_bounded(:,i4,i3))-squeeze(tissue_vel_overall_mean_bounded(:,i4,i3)),...
    'Color',[0 0 0],'CapSize',3);
plot(max_hub_size(i5),tissue_vel_overall_mean_bounded(i5,i4,i3),'ro','MarkerFaceColor','r')
set(gca,'FontName','Times New Roman')
xlabel('Proliferative hub size','Interpreter','Latex')
ylabel('Wound edge velocity (CDs/hr)','Interpreter','Latex')
axis([1 21 0 0.85])
legend([p4 p2 p3],{'Bounded growth model','Continual growth model','Free tissue'},'Location','SouthEast','EdgeColor','None','Interpreter','Latex')

% Figure 5(d)
figure
hold on
for i5 = 1:length(max_hub_size)
    if i5 == 5
        [M,p3] = contour(qvf,forces,squeeze(grads_diff_bounded(i5,:,:)),[0 0],'Color',[1 0 0]);
    else
        [M,p2] = contour(qvf,forces,squeeze(grads_diff_bounded(i5,:,:)),[0 0],'Color',[0.8 0.8 0.8]);
    end
end
[M2,p1] = contour(qvf,forces,grads_diff_continual(:,1:length(qvf)),[0 0],'Color','k','LineWidth',2);
plot(0.9,2.2,'ro','MarkerFaceColor','r')
set(gca,'FontName','Times New Roman')
xlabel('$\mu_\mathrm{QVF}$','Interpreter','Latex')
ylabel('$F_\mathrm{active}$','Interpreter','Latex')
axis([min(qvf) max(qvf) min(forces) max(forces)])
legend([p1 p2 p3],{'Continual growth model','Bounded growth model','$d_\mathrm{max} = 5$'},'Color','None','EdgeColor','None','Interpreter','Latex')
xticks([0.6 0.7 0.8 0.9 1])
yticks([0 2 4 6 8 10])