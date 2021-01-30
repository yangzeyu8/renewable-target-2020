% Plot results for renewable target
% yzy, Glasgow College, UESTC

%% Initialization
%%% 0 percent
invest_compare = zeros(8,3);
load output/target0/results.mat
invest_compare(1,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target0/', 'user', num2str(2)])
close all; 
%%% 10 percent
load output/target10/results.mat
invest_compare(2,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target10/', 'user', num2str(2)])
close all;
%%% 20 percent
load output/target20/results.mat
invest_compare(3,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target20/', 'user', num2str(2)])
close all; 
%%% 30 percent
load output/target30/results.mat
invest_compare(4,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target30/', 'user', num2str(2)])
close all; 
%%% 40 percent
load output/target40/results.mat
invest_compare(5,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target40/', 'user', num2str(2)])
close all; 
%%% 50 percent
load output/target50/results.mat
invest_compare(6,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target50/', 'user', num2str(2)])
close all; 
%%% 60 percent
load output/target60/results.mat
invest_compare(7,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target60/', 'user', num2str(2)])
close all; 
%%% 70 percent
load output/target70/results.mat
invest_compare(8,:) = invest;
variable = var;
users(variable, num_user, num_situation, 2, ['target70/', 'user', num2str(2)])
close all; 

%% Plot
fig = figure();
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
set(gca,'Box','on','TickDir','in','TickLength',[.01 .01],'XMinorTick','on', ...
    'YMinorTick','off','YGrid','on','XGrid','off','LineWidth',1);
set(gca ,'FontSize',12,  'FontName', 'Helvetica' );
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 2:167;
x = categorical({'No Target','Target 10%','Target 20%', 'Target 30%', 'Target 40%', 'Target 50%', 'Target 60%', 'Target 70%'});
x = reordercats(x,{'No Target','Target 10%','Target 20%', 'Target 30%', 'Target 40%', 'Target 50%', 'Target 60%', 'Target 70%'});
y = invest_compare';
xlabel('Renewable Energy Target', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Investment (USD)','FontName', 'Helvetica', 'FontSize', 14)
b = bar(x,y,'stacked');
set(b, {'DisplayName'}, {'Solar Generator','Wind Turbine','Battery'}')
legend('Location', 'north', 'Orientation','horizontal')
saveas(fig, './figs/investment.pdf')


