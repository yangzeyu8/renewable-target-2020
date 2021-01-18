% Solver-based quadratic programming for Operating Cost Minimization (OCM)
% Shape of param.r_S and param.r_W is 24*10
% yzy, Glasgow College, UESTC

function [var, cost] = OCME(param, Tin_ref, Tout, Tin0, Pf_ref, Pc_ref, Pil, theta, Eb0, D, M)

%% Initialization

prob = param.prob;

AeqE = [];
beqE = [];

AE = [];
bE = [];

ubE = [];
lbE = [];

HE = [];
fE = [];

for s = 1:10
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

% K. Operator Cost
idx_po = [1:24] + length(idx); idx = [idx, idx_po];
idx_pomax = 1 + length(idx); idx = [idx, idx_pomax];
idx_z = [1:24] + length(idx); idx = [idx, idx_z];
idx_zmax = 1 + length(idx); idx = [idx, idx_zmax];

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

% 7. Pil % inflexible load
% 8. theta % renewable target percentage

% E. Grid Power Supply Input Parameters

UBpg = param.UBpg;
LBpg = param.LBpg;

% F. Aggregate Power Supply

% G. Battery Input Parameters

% 9. Eb0 % initial battery energy
% 10. D % how many days do we observe

UBpch = param.UBpch;
LBpch = param.LBpch;

UBpdis = param.UBpdis;
LBpdis = param.LBpdis;

UBeb = param.UBeb;
LBeb = param.LBeb;

% H. Solar Energy Investment Input Parameters

r_S = param.r_S(:,s); % power per unit of invested soalr capacity
c_S = param.c_S; % investment cost of solar power per kW

% I. Wind Energy Investment Input Parameters

r_W = param.r_W(:,s); % power per unit of invested wind capacity
c_W = param.c_W; % investment cost of wind power per kW

% J. Battery Energy Investment Input Parameters

L_B = param.L_B; % life span energy per unit of invested battery capacity
c_B = param.c_B; % investment cost of battery power per kW

% K. Operator Cost

% 11. M % the investment budget

gamma_po = param.gamma_po;

%% Air Conditioner Calculation

% Cost Function
C = 3.3;  R = 1.35;  eta_pac = 7;
a = (1 - 1/(C*R));  b = 1/(C*R);  c = eta_pac/C;

H_pac = zeros(idx_len, idx_len);
H_pac(idx_tin,idx_tin) = 2*gamma_pac*eye(24); % param.gamma_pac

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

%% Aggregate Power Supply Calculation

% Cost Function
pi1 = 0.2; pi2 = 0.8; % grid power price
f_po = zeros(idx_len,1);
f_po(idx_z,1) = pi1*ones(24,1);
f_po(idx_zmax,1) = pi2;

% Constraint: calculation of pa
beq_pa = zeros(24,1);

Aeq_pa = zeros(24, idx_len);
Aeq_pa(:,idx_pa) = eye(24);
Aeq_pa(:,idx_pre) = -eye(24);
Aeq_pa(:,idx_pg) = -eye(24);


% Constraint: calculation of po
beq_po = zeros(24,1);

Aeq_po = zeros(24, idx_len);
Aeq_po(:,idx_po) = eye(24);
Aeq_po(:,idx_pa) = -eye(24);
% Here r_S and r_W are 24*1
Aeq_po(:,idx_beta_S) = r_S;
Aeq_po(:,idx_beta_W) = r_W;

% Constraint: limit of pomax and po
b_po = zeros(24,1);

A_po = zeros(24, idx_len);
A_po(:,idx_po) = eye(24);
A_po(:,idx_pomax) = -ones(24,1);

% Constraint: equivalent form z versus po, zmax versus pomax
% po-z <= 0
b_z = zeros(24,1);

A_z = zeros(24, idx_len);
A_z(:,idx_po) = eye(24);
A_z(:,idx_z) = -eye(24);

% z >= 0
lb_z = zeros(24,1); 

% pomax-zmax <= 0
b_zmax = zeros(1,1);

A_zmax = zeros(1, idx_len);
A_zmax(:,idx_pomax) = eye(1);
A_zmax(:,idx_zmax) = -eye(1);

% zmax >= 0
lb_zmax = zeros(1,1); 

% Constraint: boundaries of pg
lb_pg = ones(24,1)*LBpg; % param.LBpg
ub_pg = ones(24,1)*UBpg; % param.UBpg

% Constraint: calculation of battery storage eb
eta_ch = 0.9;
eta_dis = 0.9;

beq_eb = zeros(24,1);

Aeq_eb = zeros(24, idx_len);
Aeq_eb(:,idx_eb) = eye(24) + diag(-ones(1,23), -1);
Aeq_eb(:,idx_pch) = -eta_ch*eye(24);
Aeq_eb(:,idx_pdis) = 1/eta_dis*eye(24);
Aeq_eb(1,idx_eb0) = -1;

% Constraint: initial battery storage eb0
beq_eb0 = Eb0; % Eb0

Aeq_eb0 = zeros(1,idx_len);
Aeq_eb0(1,idx_eb0) = 1;

% Constraint: boundaries of battery storage eb
b_eb = zeros(24,1);

A_eb = zeros(24, idx_len);
A_eb(:, idx_eb) = eye(24);
A_eb(:, idx_beta_B) = -UBeb*ones(24,1); % param.UBeb

lb_eb = ones(24,1)*LBeb; % LBeb

% Constraint: boundaries of pch
b_pch = zeros(24,1);

A_pch = zeros(24, idx_len);
A_pch(:, idx_pch) = eye(24);
A_pch(:, idx_beta_B) = -UBpch*ones(24,1); % param.UBpch

lb_pch = ones(24,1)*LBpch; % param.LBpch

% Constraint: boundaries of pdis
b_pdis_beta = zeros(24,1);

A_pdis_beta = zeros(24, idx_len);
A_pdis_beta(:, idx_pdis) = eye(24);
A_pdis_beta(:, idx_beta_B) = -UBpdis*ones(24,1); % param.UBpis

lb_pdis = ones(24,1)*LBpdis; % param.LBpdis

% Constraint: battery life span pdis
b_pdis = zeros(1,1);

A_pdis = zeros(1, idx_len);
A_pdis(:,idx_pdis) = D*ones(24,1); % D, the number of days we observe
A_pdis(:,idx_beta_B) = -L_B;

% Constraint: boundaries of pre
lb_pre = zeros(24,1);
% ub_pre = UBpre; % UBpre

b_pre_beta = zeros(24,1);

A_pre_beta = zeros(24, idx_len);
A_pre_beta(:, idx_pre) = eye(24);
A_pre_beta(:, idx_beta_S) = -r_S;
A_pre_beta(:, idx_beta_W) = -r_W;

% Constraint: calculation of pre
beq_pre = Pil; % Pil, the inflexible load

Aeq_pre = zeros(24, idx_len);
Aeq_pre(:,idx_pre) = eye(24);
Aeq_pre(:,idx_pg) = eye(24);
Aeq_pre(:,idx_pdis) = eye(24);
Aeq_pre(:,idx_pac) = -eye(24);
Aeq_pre(:,idx_pf) = -eye(24);
Aeq_pre(:,idx_pc) = -eye(24);
Aeq_pre(:,idx_pch) = -eye(24);

% Constraint: renewable target pre
b_pre = zeros(1,1);

A_pre = zeros(1, idx_len);
A_pre(:,idx_pre) = -ones(1,24);
A_pre(:,idx_pa) = theta*ones(1,24);

% Constraint: bugdet limit
b_beta = M; % M, the total investment budget

A_beta = zeros(1, idx_len);
A_beta(:,idx_beta_S) = c_S*ones(1,1);
A_beta(:,idx_beta_W) = c_W*ones(1,1);
A_beta(:,idx_beta_B) = c_B*ones(1,1);

%% Calculation of Investment of Renewable Energy 

% Cost Function: Solar Energy
f_beta_S = zeros(idx_len, 1);
f_beta_S(idx_beta_S, 1) = c_S;

% Cost Function: Wind Energy
f_beta_W = zeros(idx_len, 1);
f_beta_W(idx_beta_W, 1) = c_W;

% Cost Function: Battery Energy
f_beta_B = zeros(idx_len, 1);
f_beta_B(idx_beta_B, 1) = c_B;

%% Optimization Setup

Aeq = [Aeq_tin0; Aeq_pac; Aeq_pf; Aeq_pre; Aeq_eb; Aeq_eb0; Aeq_pa; Aeq_po];
beq = [beq_tin0; beq_pac; beq_pf; beq_pre; beq_eb; beq_eb0; beq_pa; beq_po]; 

A = [A_beta; A_eb; A_pch; A_pdis_beta; A_pdis; A_po; A_z; A_zmax; A_pre_beta; A_pre];
b = [b_beta; b_eb; b_pch; b_pdis_beta; b_pdis; b_po; b_z; b_zmax; b_pre_beta; b_pre];

ub = ones(idx_len,1) * 500000;
ub(idx_pac, 1) = ub_pac;
ub(idx_tin, 1) = ub_tin;
ub(idx_pf, 1) = ub_pf;
ub(idx_pc, 1) = ub_pc;
ub(idx_pg, 1) = ub_pg;

lb = zeros(idx_len, 1);
lb(idx_pac, 1) = lb_pac;
lb(idx_tin, 1) = lb_tin;
lb(idx_pf, 1) = lb_pf;
lb(idx_pc, 1) = lb_pc;
lb(idx_pg, 1) = lb_pg;
lb(idx_eb, 1) = lb_eb;
lb(idx_pch, 1) = lb_pch;
lb(idx_pdis, 1) = lb_pdis;
lb(idx_pre, 1) = lb_pre;
lb(idx_z, 1) = lb_z;
lb(idx_zmax, 1) = lb_zmax;

H = D*(H_pac + H_pc + H_pf);

f = D*(f_pac + f_pc + f_pf + f_po) + f_beta_S + f_beta_W + f_beta_B;

AeqE = blkdiag(AeqE, Aeq);
beqE = [beqE; beq];

AE = blkdiag(AE, A);
bE = [bE; b];

ubE = [ubE; ub];
lbE = [lbE; lb];

HE = blkdiag(HE, prob(s)*H);

fE = [fE; prob(s)*f];
end

options = optimoptions('quadprog','Display','off');
[p_userE, ~, exitflag, ~] = quadprog(HE,fE,AE,bE,AeqE,beqE,lbE,ubE,[],options);
if exitflag ~= 1
    fprintf('[Warning] OCME does not converge\n');
end

%% Calculate the Variable

for i=1:10
    
p_user = p_userE(1+295*(i-1):295*i);

% A. Air Conditioner Var
pac = p_user(idx_pac); var.pac(:, i) = pac;
tin = p_user(idx_tin); var.tin(:, i) = tin;
tin0 = p_user(idx_tin0); var.tin0(:, i) = tin0;

% B. Flexible Load Var
pf = p_user(idx_pf); var.pf(:, i) = pf;

% C. Curtailed Load Var
pc = p_user(idx_pc); var.pc(:, i) = pc;

% D. Renewable Power Supply Var
pre = p_user(idx_pre); var.pre(:, i) = pre;

% E. Grid Power Supply Var
pg = p_user(idx_pg); var.pg(:, i) = pg;

% F. Aggregate Power Supply
pa = p_user(idx_pa); var.pa(:, i) = pa;

% G. Battery Var
eb = p_user(idx_eb); var.eb(:, i) = eb;
eb0 = p_user(idx_eb0); var.eb0(:, i) = eb0;
pch = p_user(idx_pch); var.pch(:, i) = pch;
pdis = p_user(idx_pdis); var.pdis(:, i) = pdis;

% H. Solar Energy Investment Var
beta_S = p_user(idx_beta_S); var.beta_S(:, i) = beta_S;

% I. Wind Energy Investment Var
beta_W = p_user(idx_beta_W); var.beta_W(:, i) = beta_W;

% J. Battery Energy Investment Var
beta_B = p_user(idx_beta_B); var.beta_B(:, i) = beta_B;

% K. Operator Cost
po = p_user(idx_po); var.po(:, i) = po;
pomax = p_user(idx_pomax); var.pomax(:, i) = pomax;
z = p_user(idx_z); var.z(:, i) = z;
zmax = p_user(idx_zmax); var.zmax(:, i) = zmax;

%% Calculate the Four Kinds of Cost

cost_pac = D*gamma_pac*(tin - ones(24,1)*Tin_ref)'*(tin - ones(24,1)*Tin_ref);
cost_pf = D*gamma_pf*(pf - Pf_ref)'*(pf - Pf_ref);
cost_pc = D*gamma_pc*(pc - Pc_ref)'*(pc - Pc_ref);
cost_po = D*gamma_po*(pi1*sum(z) + pi2*max(zmax));
cost_beta_S = gamma_po*c_S*beta_S;
cost_beta_W = gamma_po*c_W*beta_W;
cost_beta_B = gamma_po*c_B*beta_B;

cost(:,i) = [cost_pac; cost_pf; cost_pc; cost_po; cost_beta_S; cost_beta_W; cost_beta_B];
end
