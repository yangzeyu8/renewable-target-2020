% Inflexible load data: Austin data 2016/9/1~2016/9/7
% Temperature data: 9/6 2016 at Austin
% Renewable energy generationg data: 10 points from Hong Kong
% yzy, Glasgow College, UESTC

close all; clear; format compact; rng(2342);

%% Initialization
%%% Read Hourly Usage Data of 10 Users
hourly_usage_file = './data/hourly_usage_201609.csv';  % hourly usage
hu_table = readtable(hourly_usage_file);
hu = table2array(hu_table(:, 3));
%%% 10 Users & 10 Situations
num_user = 3;
num_situation = 10;

%% Input Parameters
%%% Initialization
Pf_ref = zeros(24,num_user,7);
Pc_ref = zeros(24,num_user,7);
Pil = zeros(24,num_user,7);
for n = 1:num_user
    for d = 1:7
        for t = 1:24
            Pil(t, n, d) = hu((n-1)*7*24 + (d-1)*24 + t);
        end
        h = Pil(:,n,d);
        avg = mean(h);
        for t = 1:24
            if h(t) < avg
                Pf_ref(t, n, d) = 0;
                Pc_ref(t, n, d) = 0;
            else
                Pf_ref(t,n,d) = (0.2+0.25*rand())*h(t);
                Pc_ref(t,n,d) = (0.2+0.25*rand())*h(t);
            end
        end
    end
end
Pil = Pil - Pf_ref - Pc_ref;
param.prob = [0.2000 0.0140 0.0440 0.12 0.025 0.0470 0.2490 0.2080 0.0680 0.0250];
%%% A. Flexible Load (i, omega): param, Pf_ref
param.gamma_pf = 0.25; 
Pf_ref = mean(Pf_ref,3);
%%% B. Curtailed Load (i, omega): param, Pc_ref
param.gamma_pc = 0.25; 
param.LBpc = 0; 
param.UBpc = 7; 
Pc_ref = mean(Pc_ref,3);
%%% C. HVAC (i, omega): param, Tin_ref, Tout, Tin0
param.gamma_pac = 0.25; 
param.LBtin = 15; 
param.UBtin = 32;
param.C = 3.3;
param.R = 1.35;
param.eta_pac = 7;
Tin_ref = [22; 23; 24];
Tout = temperature(); % Read temperature data: 9/6~9/12 2016 at Austin
Tin0 = 23;
%%% D. Inflexible Load (Input Only): Pil
Pil = mean(Pil,3);
%%% E. Renewable Energy (omega): param
param.r_S = solar(); % solar power per capacity
param.r_S = param.r_S(:,1:num_situation);
param.r_W = wind(); % wind power per capacity
param.r_W = param.r_W(:,1:num_situation);
%%% F. Power Grid (omega): param
param.UBpg = 5*num_user; 
param.LBpg = 0; 
%%% G. Battery (omega): param, Eb0, D
param.eta_ch = 0.9;
param.eta_dis = 0.9; 
param.eb_rate_min = 0.1;
param.eb_rate_max = 0.9;
param.pch_rate_max = 0.5;
param.pdis_rate_max = 0.5;
param.eb_rate_life = 100000;
Eb0 = 0;
D = 365*10;
%%% H. Aggregate Supply (omega): Pil
%%% I. Operator Cost (omega): param
param.pi1 = 0.3; 
param.pi2 = 1.2;
%%% J. Investment (operator): param, theta, Theta, M
param.c_S = 42000/10;
param.c_W = 7000/5;
param.c_B = 5000/14;
theta = 0;
Theta = 0.025;
M = num_user*10000000;
    
%% Total Simulated 10 Users
[var, cost, invest] = solver(num_user, num_situation, param, Pf_ref, Pc_ref, Tin_ref, Tout, Tin0, Pil, Eb0, D, theta, Theta, M);
save 'output/target0/results.mat'
visualization
compare
