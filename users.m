function [] = users(var ,num_user, num_situation, chosen_user, location)

%% Initialization
user = num_user;
situation = num_situation;

%% A. Flexible Load (i, omega)
x = zeros(24, situation, user);
y = var.pf;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,s,i) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,chosen_user), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Flexible Load Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', 'Scenario 6', 'Scenario 7', 'Scenario 8', 'Scenario 9', 'Scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, ['./figs/', location, '/pf.pdf'])

%% B. Curtailed Load (i, omega)
x = zeros(24, situation, user);
y = var.pc;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,s,i) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,chosen_user), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Curtailed Load Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', 'Scenario 6', 'Scenario 7', 'Scenario 8', 'Scenario 9', 'Scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, ['./figs/', location, '/pc.pdf'])

%% C. HVAC (i, omega)
%%% HVAC Power
x = zeros(24, situation, user);
y = var.pac;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,s,i) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,chosen_user), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('HVAC Consumption (kWh)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', 'Scenario 6', 'Scenario 7', 'Scenario 8', 'Scenario 9', 'Scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5); 
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, ['./figs/', location, '/pac.pdf'])
%%% Indoor Temperature
x = zeros(24, situation, user);
y = var.tin;
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
x(:,s,i) = y(start_idx:end_idx);
end
end
lw =2; ms = 1;
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
plot(x(:,:,chosen_user), 'LineWidth', lw)
xlim([1 24]);
xticks([1:12:24,24]);
xticklabels({'6 SEP','13:00','24:00'});
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
xlabel('Date', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Indoor Temperature (\circC)','FontName', 'Helvetica', 'FontSize', 14)
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
h = legend({'Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', 'Scenario 6', 'Scenario 7', 'Scenario 8', 'Scenario 9', 'Scenario 10'} ...
    , 'Location', 'northoutside', 'Orientation','horizontal',"NumColumns",5);  
set(h,'FontName', 'Helvetica', 'FontSize', 12);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
saveas(fig, ['./figs/', location, '/tin.pdf'])
end

