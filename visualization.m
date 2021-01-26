% Plot results for OCM
% yzy, Glasgow College, UESTC

%% Initialization
%%% Read Results
close all; clear;
load output/target50/results.mat
state = num_user*num_situation;
user = num_user;
situation = num_situation;

%% A. Flexible Load (i, omega)
x = zeros(24, user, situation);
y = var.pf;
scenario = 1;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,i,s) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,scenario), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Flexible Load Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'user 1', 'user 2', 'user 3', 'user 4', 'user 5', 'user 6', 'user 7', 'user 8', 'user 9', 'user 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pf.pdf')

%% B. Curtailed Load (i, omega)
x = zeros(24, user, situation);
y = var.pc;
scenario = 1;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,i,s) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,scenario), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Curtailed Load Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'user 1', 'user 2', 'user 3', 'user 4', 'user 5', 'user 6', 'user 7', 'user 8', 'user 9', 'user 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pc.pdf')

%% C. HVAC (i, omega)
%%% HVAC Power
x = zeros(24, user, situation);
y = var.pac;
scenario = 1;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,i,s) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,scenario), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('HVAC Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'user 1', 'user 2', 'user 3', 'user 4', 'user 5', 'user 6', 'user 7', 'user 8', 'user 9', 'user 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pac.pdf')
%%% Indoor Temperature
x = zeros(24, user, situation);
y = var.tin;
scenario = 1;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,i,s) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,scenario), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Indoor Temperature (\circC)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'user 1', 'user 2', 'user 3', 'user 4', 'user 5', 'user 6', 'user 7', 'user 8', 'user 9', 'user 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/tin.pdf')

%% D. Inflexible Load (input only)

%% E. Renewable Energy (omega)
x = zeros(24, situation);
y = var.pre;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Renewable Energy (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pre.pdf')

%% F. Power Grid (omega)
x = zeros(24, situation);
y = var.pg;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Power Grid (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pg.pdf')

%% G. Battery (omega)
%%% Battery
x = zeros(24, situation);
y = var.eb;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Battery Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/eb.pdf')
%%% Charge
x = zeros(24, situation);
y = var.pch;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Battery Charge (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pch.pdf')
%%% Discharge
x = zeros(24, situation);
y = var.pdis;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Battery Discharge (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pdis.pdf')

%% H. Aggregate Supply (omega)
x = zeros(24, situation);
y = var.pa;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Aggregate Supply (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/pa.pdf')

%% I. Operator Cost (omega)
x = zeros(24, situation);
y = var.po;
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
x(:,s) = y(start_idx:end_idx);
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x, 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Costly Power (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'scenario 1', 'scenario 2', 'scenario 3', 'scenario 4', 'scenario 5', 'scenario 6', 'scenario 7', 'scenario 8', 'scenario 9', 'scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, './figs/target50/po.pdf')

%% J. Investment (operator)
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
x = categorical({'Solar Investment','Wind Investment','Battery Investment'});
x = reordercats(x,{'Solar Investment','Wind Investment','Battery Investment'});
y = invest;
bar(x,y)
saveas(fig, './figs/target50/invest.pdf')