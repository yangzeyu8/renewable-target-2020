% Inflexible load data: Austin data 2016/9/1~2016/9/7
% Temperature data: 9/6 2016 at Austin
% Renewable energy generationg data: 10 points from Hong Kong
% yzy, Glasgow College, UESTC

close all; clear; format compact; rng(2342);

%% Initialization

% Read hourly usage data of 10 users
hourly_usage_file = './data/hourly_usage_201609.csv';  % hourly usage
hu_table = readtable(hourly_usage_file);
hu = table2array(hu_table(:, 3));

num_user = 1;

%% Input Parameters

% 1. param % parameters for optimization
param.gamma_pac = 0.5; % discomfort coefficient for tin
param.UBpac = 10; % upper boundary of pac
param.LBpac = 0; % lower boundary of pac

param.UBtin = 32; % upper boundary of tin
param.LBtin = 15; % lower boundary of tin 

param.gamma_pf = 0.3; % discomfort coefficient for pf
param.UBpf = 7; % upper boundary flexible load per hour pf
param.LBpf = 0; % lower boundary flexible load per hour pf

param.gamma_pc = 0.4; % discomfort coefficient for pc
param.UBpc = 7; % upper boundary curtailed load per hour pc
param.LBpc = 0; % lower boundary curtailed load per hour pc

param.UBpg = 20; % upper boundary grid load per hour pg
param.LBpg = 0; % loweer boundary grid load per hour pg

param.UBpch = 10; % upper boundary battery charge per capacity pch
param.LBpch = 0; % lower boundary battery charge per capacity pch

param.UBpdis = 10; % upper boundary battery discharge per capacity pdis
param.LBpdis = 0; % lower boundary battery discharge per capacity pdis

param.UBeb = 10; % upper boundary battery storage per capacity eb
param.LBeb = 0; % lower boundary battery storage per capacity eb

param.r_S = solar(); % solar power per capacity r_S, r_S is 24*1
param.c_S = 5000; % solar cost per capacity c_S

param.r_W = wind(); % wind power per capacity r_W, r_W is 24*1
param.c_W = 2000; % wind cost per capacity c_W

param.L_B = 100000; % battery life span L_B
param.c_B = 3000;  % battery cost per capacity c_B

param.gamma_po = 1; % discomfort coefficient for money

param.prob = [0.2 0.01369863 0.043835616 0.120547945 0.024657534 0.046575342 0.249315068 0.208219178 0.068493151 0.024657534];

% 2. Tin_ref % preferred temperature
Tin_ref = 24;

% 3. Tout % outdoor temperature of one day
Tout = temperature(); % Read temperature data: 9/6~9/12 2016 at Austin

% 4. Tin0 % indoor initial temperature
Tin0 = 23*ones(num_user, 1);

% 5. Pf_ref % preferred flexible load
% 6. Pc_ref % preferred curtailed load
% 7. Pil % inflexible load
Pf_ref = zeros(24,10,7);
Pc_ref = zeros(24,10,7);
Pil = zeros(24,10,7);
for n = 1:10
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
Pil = mean(mean(Pil,3),2);
Pf_ref = mean(mean(Pf_ref,3),2);
Pc_ref = mean(mean(Pc_ref,3),2);

% 8. theta % renewable target percentage
theta = 0.5;

% 9. Eb0 % initial battery energy
Eb0 = 0.1;

% 10. D % how many days do we observe
D = 365*10;

% 11. M % the investment budgets
M = 5*10000;

%% total simulated 1 user

[var, cost] = OCM(param, Tin_ref, Tout, Tin0, Pf_ref, Pc_ref, Pil, theta, Eb0, D, M);
[varE, costE] = OCME(param, Tin_ref, Tout, Tin0, Pf_ref, Pc_ref, Pil, theta, Eb0, D, M);

save 'output/results.mat'


