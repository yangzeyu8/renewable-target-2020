% solver-based quadratic programming for Operating Cost Minimization (OCM)
% This function returns struct var and struct cost
% var has 11 attributes
% cost has 4 attributes
% lb_tin and ub_tin are 24*1, gamma_pac is 24*24
% yzy 2020/12/21

function [var, cost] = OCM(param, Tin_ref, Tout, Tin0, Pf_ref, Pc_ref, UBpre, Pil, Eb0)

%% User Model Variable Index
idx = [];

% A. Air Conditioner 
idx_tin = [1:24]; idx = [idx, idx_tin];
idx_tin0 = 1 + length(idx); idx = [idx, idx_tin0];
idx_pac = [1:24] + length(idx); idx = [idx, idx_pac];

% B. Flexible Load
idx_pf = [1:24] + length(idx); idx = [idx, idx_pf];

% C. Curtailed Load
idx_pc = [1:24] + length(idx); idx = [idx, idx_pc];

%% Operator Model Variable Index

% D. Renewable Power Supply
idx_pre = [1:24] + length(idx); idx = [idx, idx_pre];

% E. Grid Power Supply
idx_pg = [1:24] + length(idx); idx = [idx, idx_pg];

% F. Aggregate Power Supply
idx_pa = [1:24] + length(idx); idx = [idx, idx_pa];
idx_pamax = 1 + length(idx); idx = [idx, idx_pamax];

% G. Battery
idx_eb = [1:24] + length(idx); idx = [idx, idx_eb];
idx_eb0 = 1 + length(idx); idx = [idx, idx_eb0];
idx_pch = [1:24] + length(idx); idx = [idx, idx_pch];
idx_pdis = [1:24] + length(idx); idx = [idx, idx_pdis];

% H. Solar Energy Investment
idx_beta_S = 1 + length(idx); idx = [idx, idx_beta_S];

% I. Wind Energy Investment
idx_beta_W = 1 + length(idx); idx = [idx, idx_beta_W];

% J. Battery Energy Investment
idx_beta_B = 1 + length(idx); idx = [idx, idx_beta_B];

% Final Variable Length
idx_len = length(idx);

%% Input Parameters

% A. Air Conditioner Input Parameters

% 1. param % parameters for optimization
% 2. Tin_ref % preferred temperature
% 3. Tout % outdoor temperature of one day
% 4. Tin0 % indoor initial temperature

gamma_pac = param.gamma_pac;
UBpac = param.UBpac;
LBpac = param.LBpac;

LBtin = param.LBtin;
UBtin = param.UBtin;

% B. Flexible Load Input Parameters

% 5. Pf_ref % preferred flexible load

gamma_pf = param.gamma_pf;
UBpf = param.UBpf;
LBpf = param.LBpf;

% C. Curtailed Load Input Parameters

% 6. Pc_ref % preferred curtailed load

gamma_pc = param.gamma_pc;
UBpc = param.UBpc;
LBpc = param.LBpc;

% D. Renewable Power Supply Input Parameters

% 7. UBpre % upper boundary of renewable energy
% 8. Pil % inflexible load

% E. Grid Power Supply Input Parameters

gamma_pg = param.gamma_pg;
UBpg = param.UBpg;
LBpg = param.LBpg;

% F. Aggregate Power Supply

% G. Battery Input Parameters

% 9. Eb0 % initial battery energy

UBpch = param.UBpch;
LBpch = param.LBpch;

UBpdis = param.UBpdis;
LBpdis = param.LBpdis;

UBeb = param.UBeb;
LBeb = param.LBeb;

% H. Solar Energy Investment Input Parameters

r_S = param.r_S; % power per unit of invested soalr capacity
c_S = param.c_S; % investment cost of solar power per kW

% I. Wind Energy Investment Input Parameters

r_W = param.r_W; % power per unit of invested wind capacity
c_W = param.c_W; % investment cost of wind power per kW

% J. Battery Energy Investment Input Parameters

L_B = param.L_B; % life span energy per unit of invested battery capacity
c_B = param.c_B; % investment cost of battery power per kW

%% Air Conditioner Calculation

% Cost Function
C = 3.3;  R = 1.35;  eta_pac = 7;
a = (1 - 1/(C*R));  b = 1/(C*R);  c = eta_pac/C;

H_pac = zeros(idx_len, idx_len);

% Here gamma_pac is 24*24
H_pac(idx_tin,idx_tin) = 2*gamma_pac; % param.gamma_pac

f_pac = zeros(idx_len, 1);
f_pac(idx_tin,1) = -2*gamma_pac*Tin_ref*ones(24,1); % Tin_ref

% Constraint: relation among pac, tin and tin0
beq_pac = b*Tout; % Tout

Aeq_pac = zeros(24, idx_len);
Aeq_pac(:,idx_tin) = eye(24) + diag(-a*ones(1,23),-1);
Aeq_pac(:,idx_pac) = eye(24)*c;
Aeq_pac(1,idx_tin0) = -a;

beq_tin0 = Tin0; % Tin0

Aeq_tin0 = zeros(1,idx_len);
Aeq_tin0(1,idx_tin0) = 1;

% Constraint: boundaries of pac
lb_pac = ones(24,1)*LBpac; % param.LBpac
ub_pac = ones(24,1)*UBpac; % param.UBpac

% Constraint: boundaries of tin
% Here lb_tin and ub_tin are 24*1
lb_tin = LBtin; % param.LBtin
ub_tin = UBtin; % param.UBtin

%% Flexible Load Calculation

% Cost Function
H_pf = zeros(idx_len, idx_len);
H_pf(idx_pf, idx_pf) = 2*gamma_pf*eye(24); % param.gamma_pf

f_pf = zeros(idx_len, 1);
f_pf(idx_pf, 1) = -2*gamma_pf*Pf_ref; % Pf_ref

% Constraint: sum of flexible load
Df = sum(Pf_ref);
beq_pf = Df;

Aeq_pf = zeros(1,idx_len);
Aeq_pf(1,idx_pf) = ones(1,24);

% Constraint: boundaries of flexible load
lb_pf = ones(24,1)*LBpf; % param.LBpf
ub_pf = ones(24,1)*UBpf; % param.UBpf

%% Curtailed Load Calculation

% Cost Function
H_pc = zeros(idx_len, idx_len);
H_pc(idx_pc, idx_pc) = 2*gamma_pc*eye(24); % param.gamma_pc

f_pc = zeros(idx_len, 1);
f_pc(idx_pc, 1) = -2*gamma_pc*Pc_ref; % Pc_ref

% Constraint: boundaries of curtailed load
lb_pc = ones(24,1)*LBpc; % param.LBpc
ub_pc = ones(24,1)*UBpc; % param.UBpc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make changes
%% Aggregate Power Supply Calculation

% Cost Function
pi1 = 0.2; pi2 = 0.8; % grid power price
f_pg = zeros(idx_len,1);
f_pg(idx_pg,1) = pi1*ones(24,1);
f_pg(idx_pgmax,1) = pi2;

% Constraint (7)
lb_pg = zeros(24,1);
ub_pg = ones(24,1)*UBpg; % param.UBpg

% Constraint for idx_pgmax
b_pg = zeros(24,1);

A_pg = zeros(24, idx_len);
A_pg(:,idx_pg) = eye(24);
A_pg(:,idx_pgmax) = -ones(24,1);

%% Grid Cost

% Cost Function
pi1 = 0.2; pi2 = 0.8; % grid power price
f_pg = zeros(idx_len,1);
f_pg(idx_pg,1) = pi1*ones(24,1);
f_pg(idx_pgmax,1) = pi2;

% Constraint (7)
lb_pg = zeros(24,1);
ub_pg = ones(24,1)*UBpg; % param.UBpg

% Constraint for idx_pgmax
b_pg = zeros(24,1);

A_pg = zeros(24, idx_len);
A_pg(:,idx_pg) = eye(24);
A_pg(:,idx_pgmax) = -ones(24,1);

%% Battery Cost Calculation

% Cost Function
eta_ch = 0.9;
eta_dis = 0.9;

H_pch = zeros(idx_len, idx_len);
H_pch(idx_pch, idx_pch) = 2*gamma_eb*eye(24); % param.gamma_eb

H_pdis = zeros(idx_len, idx_len);
H_pdis(idx_pdis, idx_pdis) = 2*gamma_eb*eye(24);

% f_eb = zeros(idx_len,1);
% f_eb(idx_pch,1) = gamma_eb*ones(24,1);
% f_eb(idx_pdis,1) = gamma_eb*ones(24,1);

% Constraint (10)
beq_eb = zeros(24,1);

Aeq_eb = zeros(24, idx_len);
Aeq_eb(:,idx_eb) = eye(24) + diag(-ones(1,23), -1);
Aeq_eb(:,idx_pch) = -eta_ch*eye(24);
Aeq_eb(:,idx_pdis) = 1/eta_dis*eye(24);
Aeq_eb(1,idx_eb0) = -1;

% initial battery storage
beq_eb0 = Eb0; % Eb0

Aeq_eb0 = zeros(1,idx_len);
Aeq_eb0(1,idx_eb0) = 1;

% Constraint (11)
lb_eb = zeros(24,1);
ub_eb = ones(24,1)*UBeb; % UBeb

% Constraint (12)
lb_pch = zeros(24,1);
ub_pch = ones(24,1)*UBpch; % param.UBpch

% Constraint (13)
lb_pdis = zeros(24,1);
ub_pdis = ones(24,1)*UBpdis; % param.UBpdis

%% Renewable Energy

% Constraint (8)
lb_pre = zeros(24,1);
ub_pre = UBre; % UBre

% Constraint (17)
beq_pre = Pil; % Pil, the inflexible load

Aeq_pre = zeros(24, idx_len);
Aeq_pre(:,idx_pre) = eye(24);
Aeq_pre(:,idx_pg) = eye(24);
Aeq_pre(:,idx_pdis) = eye(24);
Aeq_pre(:,idx_pac) = -eye(24);
% Aeq_pre(:,idx_pf) = -eye(24);
Aeq_pre(:,idx_pch) = -eye(24);

%% Optimization Setup: Eq(18)
H = H_pch + H_pdis;
f = f_z + f_pg + f_DR;

Aeq = [Aeq_pac; Aeq_tin0; Aeq_DR; Aeq_eb; Aeq_eb0; Aeq_pre];
beq = [beq_pac; beq_tin0; beq_DR; beq_eb; beq_eb0; beq_pre];

A = [A_z1; A_z2; A_pg];
b = [b_z1; b_z2; b_pg];

ub = ones(idx_len,1) * 500;
ub(idx_pac, 1) = ub_pac;
ub(idx_tin, 1) = ub_tin;
ub(idx_tref, 1) = ub_tref;
% ub(idx_pf, 1) = ub_pf;
ub(idx_pg, 1) = ub_pg;
ub(idx_eb, 1) = ub_eb;
ub(idx_pch, 1) = ub_pch;
ub(idx_pdis, 1) = ub_pdis;
ub(idx_pre, 1) = ub_pre;

lb = zeros(idx_len, 1);
lb(idx_pac, 1) = lb_pac;
lb(idx_tin, 1) = lb_tin;
lb(idx_tref, 1) = lb_tref;
lb(idx_z, 1) = lb_z;
% lb(idx_pf, 1) = lb_pf;
lb(idx_pg, 1) = lb_pg;
lb(idx_eb, 1) = lb_eb;
lb(idx_pch, 1) = lb_pch;
lb(idx_pdis, 1) = lb_pdis;
lb(idx_pre, 1) = lb_pre;

options = optimoptions('quadprog','Display','off');
[p_user, fval, exitflag, output] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
if exitflag ~= 1
    fprintf('[Warning] emp does not converge\n');
end

%% Calculate the Var and Cost

% Air Conditioner Var
tin = p_user(idx_tin); var.tin = tin;
tref = p_user(idx_tref); var.tref = tref;
pac = p_user(idx_pac); var.pac = pac;
tin0 = p_user(idx_tin0); var.tin0 = tin0;
z = p_user(idx_z); var.z = z;
pdr = p_user(idx_pdr); var.pdr = pdr;

% Flexible Load Var
% pf = p_user(idx_pf); var.pf = pf;

% Grid Cost Var
pg = p_user(idx_pg); var.pg = pg;
pgmax = p_user(idx_pgmax); var.pgmax = pgmax;

% Battery Var
eb = p_user(idx_eb); var.eb = eb;
eb0 = p_user(idx_eb0); var.eb0 = eb0;
pch = p_user(idx_pch); var.pch = pch;
pdis = p_user(idx_pdis); var.pdis = pdis;

% Renewable Energy Var
pre = p_user(idx_pre); var.pre = pre;

% Four Kinds of Cost
cost_pac = sum(gamma_pac*z);
% cost_pf = gamma_pf*(pf - Pref)'*(pf - Pref);
cost_pg = gamma_pg*(pi1*sum(pg) + pi2*max(pg));
cost_eb = gamma_eb*(sum(pch)+sum(pdis));
reduction = -pi_DR*sum(pdr);

% cost.pac = cost_pac;
% cost.pf  = cost_pf;
% cost.pg = cost_pg;
% cost.eb = cost_eb;

cost = [cost_pac; cost_pg; cost_eb; reduction];



