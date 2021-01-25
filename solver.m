% Solver-based quadratic programming for Operating Cost Minimization (OCM)
% We have 10 situations in HK
% yzy, Glasgow College, UESTC

function [var, cost, invest] = solver(num_user, num_situation, param, Pf_ref, Pc_ref, Tin_ref, Tout, Tin0, Pil, Eb0, D, theta, Theta, M)

%% User Model Variable Index
%%% Initialization
state = num_user*num_situation;
user = num_user;
situation = num_situation;
idx = [];
%%% A. Flexible Load (i, omega)
idx_pf = [1:24*state] + length(idx); idx = [idx, idx_pf];
%%% B. Curtailed Load (i, omega)
idx_pc = [1:24*state] + length(idx); idx = [idx, idx_pc];
%%% C. HVAC (i, omega)
idx_tin = [1:24*state] + length(idx); idx = [idx, idx_tin];
idx_tin0 = [1:state] + length(idx); idx = [idx, idx_tin0];
idx_pac = [1:24*state] + length(idx); idx = [idx, idx_pac];
%%% D. Inflexible Load (input only)
%%% E. Renewable Energy (omega)
idx_pre = [1:24*situation] + length(idx); idx = [idx, idx_pre];
%%% F. Power Grid (omega)
idx_pg = [1:24*situation] + length(idx); idx = [idx, idx_pg];
%%% G. Battery (omega)
idx_eb = [1:24*situation] + length(idx); idx = [idx, idx_eb];
idx_eb0 = [1:situation] + length(idx); idx = [idx, idx_eb0];
idx_pch = [1:24*situation] + length(idx); idx = [idx, idx_pch];
idx_pdis = [1:24*situation] + length(idx); idx = [idx, idx_pdis];
%%% H. Aggregate Supply (omega)
idx_pa = [1:24*situation] + length(idx); idx = [idx, idx_pa];
%%% I. Operator Cost (omega)
idx_po = [1:24*situation] + length(idx); idx = [idx, idx_po];
idx_pomax = [1:situation] + length(idx); idx = [idx, idx_pomax];
idx_z = [1:24*situation] + length(idx); idx = [idx, idx_z];
idx_zmax = [1:situation] + length(idx); idx = [idx, idx_zmax];
%%% J. Investment (operator)
idx_beta_S = 1 + length(idx); idx = [idx, idx_beta_S];
idx_beta_W = 1 + length(idx); idx = [idx, idx_beta_W];
idx_beta_B = 1 + length(idx); idx = [idx, idx_beta_B];
%%% Finalization
idx_len = length(idx);

%% Calculations
%%% Read Probability
prob = param.prob;

%% A. Flexible Load (i, omega)
%%% Function Input: param, Pf_ref
gamma_pf = param.gamma_pf; % coefficient of cost function
%%% Constraint Equation (1): Sum of Flexible Load
Aeq_pf = zeros(state,idx_len);
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
Aeq_pf(c,idx_pf(start_idx:end_idx)) = ones(1,24);
end
end
Df = sum(Pf_ref)'; % all users
beq_pf = repmat(Df,[situation,1]); % all users in each situation
%%% Cost Function Equation (2): Discomfort Cost of Flexible Load
H_pf = zeros(idx_len, idx_len);
f_pf = zeros(idx_len, 1);
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
H_pf(idx_pf(start_idx:end_idx), idx_pf(start_idx:end_idx)) = 2*gamma_pf*eye(24)*prob(s);
f_pf(idx_pf(start_idx:end_idx), 1) = -2*gamma_pf*Pf_ref(:,i)*prob(s);
end
end

%% B. Curtailed Load (i, omega)
%%% Function Input: param, Pc_ref
gamma_pc = param.gamma_pc;
LBpc = param.LBpc;
UBpc = param.UBpc;
%%% Constraint Equation (3): Boundaries of Curtailed Load
lb_pc = LBpc;
ub_pc = UBpc;
%%% Cost Function Equation (4): Discomfort Cost of Curtailed Load
H_pc = zeros(idx_len, idx_len);
f_pc = zeros(idx_len, 1);
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
H_pc(idx_pc(start_idx:end_idx), idx_pc(start_idx:end_idx)) = 2*gamma_pc*eye(24)*prob(s);
f_pc(idx_pc(start_idx:end_idx), 1) = -2*gamma_pc*Pc_ref(:,i)*prob(s);
end
end

%% C. HVAC (i, omega)
%%% Function Input: param, Tin_ref, Tout, Tin0
gamma_pac = param.gamma_pac;
LBtin = param.LBtin; % 15
UBtin = param.UBtin; % 32
C = param.C; % 3.3
R = param.R; % 1.35
eta_pac = param.eta_pac; % 7
%%% Constraint Euqation (5): Indoor, Outdoor Temperature and HVAC
a = (1 - 1/(C*R));  b = 1/(C*R);  d = eta_pac/C;
Aeq_pac = zeros(24*state,idx_len);
for s = 1:situation
for i = 1:user
c = user*(s-1)+i; 
start_idx = 1+24*(c-1);
end_idx = 24*c;
Aeq_pac(start_idx:end_idx,idx_tin(start_idx:end_idx)) = eye(24) + diag(-a*ones(1,23),-1);
Aeq_pac(start_idx:end_idx,idx_pac(start_idx:end_idx)) = d*eye(24);
Aeq_pac(start_idx,idx_tin0(c)) = -a;
end
end
beq_pac = b*repmat(Tout,[state,1]);
Aeq_tin0 = zeros(1*state,idx_len);
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
Aeq_tin0(c,idx_tin0(c)) = 1;
end
end
beq_tin0 = repmat(Tin0,[state,1]); % Tin0
%%% Constraint Equation (6): Boundaries of Indoor Tempature
lb_tin = LBtin; 
ub_tin = UBtin; 
%%% Cost Function Euqation (7): Discomfort Cost of HVAC
H_pac = zeros(idx_len, idx_len);
f_pac = zeros(idx_len, 1);
for s = 1:situation
for i = 1:user
c = user*(s-1)+i;
start_idx = 1+24*(c-1);
end_idx = 24*c;
H_pac(idx_tin(start_idx:end_idx), idx_tin(start_idx:end_idx)) = 2*gamma_pf*eye(24)*prob(s);
f_pac(idx_tin(start_idx:end_idx), 1) = -2*gamma_pac*Tin_ref*ones(24,1)*prob(s);
end
end

%% D. Inflexible Load (Input Only)
%%% Function Input: Pil

%% E. Renewable Energy (omega)
%%% Function Input: param
r_S = param.r_S;
r_W = param.r_W;
r_S_max = max(r_S);
r_W_max = max(r_W);
%%% Constraint Equation (8~9): Boundaries of Renewable Energy
A_pre_beta = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_pre_beta(start_idx:end_idx, idx_pre(start_idx:end_idx)) = eye(24);
A_pre_beta(start_idx:end_idx, idx_beta_S) = -r_S(:,s);
A_pre_beta(start_idx:end_idx, idx_beta_W) = -r_W(:,s);
end
b_pre_beta = zeros(24*situation,1);

%% F. Power Grid (omega)
%%% Function Input: param
LBpg = param.LBpg;
UBpg = param.UBpg;
%%% Constraint Equation (10): Boundaries of Power Grid
lb_pg = LBpg;
ub_pg = UBpg;

%% G. Battery (omega)
%%% Function Input: param, Eb0, D
eta_ch = param.eta_ch; % 0.9
eta_dis = param.eta_dis; % 0.9
eb_rate_max = param.eb_rate_max;
pch_rate_max = param.pch_rate_max;
pdis_rate_max = param.pdis_rate_max;
eb_rate_life = param.eb_rate_life;
%%% Constraint Equation (11): Energy Storage versus Charge/Discharge
Aeq_eb = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
Aeq_eb(start_idx:end_idx, idx_eb(start_idx:end_idx)) = eye(24) + diag(-ones(1,23), -1);
Aeq_eb(start_idx:end_idx, idx_pch(start_idx:end_idx)) = -eta_ch*eye(24);
Aeq_eb(start_idx:end_idx, idx_pdis(start_idx:end_idx)) = 1/eta_dis*eye(24);
Aeq_eb(start_idx, idx_eb0(s)) = -1;
end
beq_eb = zeros(24*situation,1);
Aeq_eb0 = zeros(situation,idx_len);
for s = 1:situation
Aeq_eb0(s,idx_eb0(s)) = 1;
end
beq_eb0 = repmat(Eb0,[situation,1]); 
%%% Constraint Equation (12): Boundaries of Battery Storage
A_eb = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_eb(start_idx:end_idx, idx_eb(start_idx:end_idx)) = eye(24);
A_eb(start_idx:end_idx, idx_beta_B) = -eb_rate_max;
end
b_eb = zeros(24*situation,1);
%%% Constraint Equation (13): Boundaries of Charging Power
A_pch = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_pch(start_idx:end_idx, idx_pch(start_idx:end_idx)) = eye(24);
A_pch(start_idx:end_idx, idx_beta_B) = -pch_rate_max;
end
b_pch = zeros(24*situation,1);
%%% Constraint Equation (14): Boundaries of Discharging Power
A_pdis = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_pdis(start_idx:end_idx, idx_pdis(start_idx:end_idx)) = eye(24);
A_pdis(start_idx:end_idx, idx_beta_B) = -pdis_rate_max;
end
b_pdis = zeros(24*situation,1);
%%% Constraint Equation (15): Battery Lifespan Limit
A_eb_L  = zeros(1, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_eb_L(1, idx_pdis(start_idx:end_idx)) = ones(1,24)*prob(s)*D;
end
A_eb_L(1, idx_beta_B) = -eb_rate_life;
b_eb_L = zeros(1,1);

%% H. Aggregate Supply (omega)
%%% Function Input: Pil
%%% Constraint Equation (16): Power Balance
Aeq_balance = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
start_user_idx = 1+24*(s-1)*user;
end_user_idx = 24*s*user;
Aeq_balance(start_idx:end_idx, idx_pre(start_idx:end_idx)) = eye(24);
Aeq_balance(start_idx:end_idx, idx_pg(start_idx:end_idx)) = eye(24);
Aeq_balance(start_idx:end_idx, idx_pdis(start_idx:end_idx)) = eye(24);
Aeq_balance(start_idx:end_idx, idx_pch(start_idx:end_idx)) = -eye(24);
Aeq_balance(start_idx:end_idx, idx_pf(start_user_idx:end_user_idx)) = -repmat(eye(24),[1,user]);
Aeq_balance(start_idx:end_idx, idx_pc(start_user_idx:end_user_idx)) = -repmat(eye(24),[1,user]);
Aeq_balance(start_idx:end_idx, idx_pac(start_user_idx:end_user_idx)) = -repmat(eye(24),[1,user]);
end
beq_balance = repmat(sum(Pil, 2),[situation,1]);
%%% Constraint Equation (17~18): Aggregate Supply Calculation
Aeq_pa = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
Aeq_pa(start_idx:end_idx,idx_pa(start_idx:end_idx)) = eye(24);
Aeq_pa(start_idx:end_idx,idx_pre(start_idx:end_idx)) = -eye(24);
Aeq_pa(start_idx:end_idx,idx_pg(start_idx:end_idx)) = -eye(24);    
end
beq_pa = zeros(24*situation,1);

%% I. Operator Cost (omega)
%%% Function Input: param
pi1 = param.pi1; % 0.2
pi2 = param.pi2; % 0.8
%%% Invisible Constraint: Calculation of Operation Power
Aeq_po = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
Aeq_po(start_idx:end_idx, idx_po(start_idx:end_idx)) = eye(24);
Aeq_po(start_idx:end_idx, idx_pa(start_idx:end_idx)) = -eye(24);
Aeq_po(start_idx:end_idx, idx_beta_S) = r_S(:,s);
Aeq_po(start_idx:end_idx, idx_beta_W) = r_W(:,s);
end
beq_po = zeros(24*situation,1);
%%% Invisible Constraint: Maximum of Operation Power
A_po = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_po(start_idx:end_idx,idx_po(start_idx:end_idx)) = eye(24);
A_po(start_idx:end_idx,idx_pomax(s)) = -ones(24,1);
end
b_po = zeros(24*situation,1);
%%% Invisible Constraint: Equivilent Form of Operation Power
%%%% po-z <= 0 and z>=0
A_z = zeros(24*situation, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_z(start_idx:end_idx,idx_po(start_idx:end_idx)) = eye(24);
A_z(start_idx:end_idx,idx_z(start_idx:end_idx)) = -eye(24);
end
b_z = zeros(24*situation,1); 
%%%% pomax-zmax <= 0 and zmax>=0
A_zmax = zeros(1*situation, idx_len);
for s = 1:situation
A_zmax(s,idx_pomax(s)) = eye(1);
A_zmax(s,idx_zmax(s)) = -eye(1);
end
b_zmax = zeros(1*situation,1); 
%%% Cost Function Equation (19): Cost of Operator
f_po = zeros(idx_len,1);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
f_po(idx_z(start_idx:end_idx), 1) = pi1*ones(24,1)*prob(s);
f_po(idx_zmax(s), 1) = pi2*prob(s);
end

%% J. Investment (operator)
%%% Function Input: param, theta, Theta, M
c_S = param.c_S;
c_W = param.c_W;
c_B = param.c_B;
%%% Constraint Equation (20): Renewable Energy Target
A_pre = zeros(1, idx_len);
for s = 1:situation
start_idx = 1+24*(s-1);
end_idx = 24*s;
A_pre(1,idx_pa(start_idx:end_idx)) = theta*ones(1,24)*prob(s);
A_pre(1,idx_pre(start_idx:end_idx)) = -ones(1,24)*prob(s);
A_pre(1,idx_pg(start_idx:end_idx)) = -Theta*ones(1,24)*prob(s);
end
b_pre = zeros(1,1);
%%% Cost Function Equation (21): Investment of Renewable Energy
f_beta = zeros(idx_len,1);
f_beta(idx_beta_S, 1) = c_S;
f_beta(idx_beta_W, 1) = c_W;
f_beta(idx_beta_B, 1) = c_B;
%%% Constraint Equation (24~25): Budget Limit
A_beta = zeros(1, idx_len);
A_beta(1, idx_beta_S) = c_S*ones(1,1);
A_beta(1, idx_beta_W) = c_W*ones(1,1);
A_beta(1, idx_beta_B) = c_B*ones(1,1);
b_beta = M; 

%% Optimization Setup
%%% Relations
Aeq = [Aeq_pf; Aeq_pac; Aeq_tin0; Aeq_eb; Aeq_eb0; Aeq_balance; Aeq_pa; Aeq_po];
beq = [beq_pf; beq_pac; beq_tin0; beq_eb; beq_eb0; beq_balance; beq_pa; beq_po]; 
A = [A_pre_beta; A_eb; A_pch; A_pdis; A_eb_L; A_po; A_z; A_zmax; A_pre; A_beta];
b = [b_pre_beta; b_eb; b_pch; b_pdis; b_eb_L; b_po; b_z; b_zmax; b_pre; b_beta];
%%% Upper Boundaries
ub = ones(idx_len,1) * 500;
ub(idx_pc, 1) = ub_pc;
ub(idx_tin, 1) = ub_tin;
ub(idx_pg, 1) = ub_pg;
%%% Lower Boundaries
lb = zeros(idx_len, 1);
lb(idx_pc, 1) = lb_pc;
lb(idx_tin, 1) = lb_tin;
lb(idx_pg, 1) = lb_pg;
%%% Objective Function
H = D*(H_pf + H_pc + H_pac);
f = D*(f_pf + f_pc + f_pac + f_po) + f_beta;
%%% Start Optimization
options = optimoptions('quadprog','Display','off');
[p_user, ~, exitflag, ~] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
if exitflag ~= 1
    fprintf('[Warning] Solver does not converge\n');
end

%% Calculate the Variable
%%% Initialization
%%% A. Flexible Load (i, omega)
pf = p_user(idx_pf); var.pf = pf;
%%% B. Curtailed Load (i, omega)
pc = p_user(idx_pc); var.pc = pc;
%%% C. HVAC (i, omega)
pac = p_user(idx_pac); var.pac = pac;
tin = p_user(idx_tin); var.tin = tin;         
tin0 = p_user(idx_tin0); var.tin0 = tin0;
%%% D. Inflexible Load (input only)
%%% E. Renewable Energy (omega)
pre = p_user(idx_pre); var.pre = pre;
%%% F. Power Grid (omega)
pg = p_user(idx_pg); var.pg = pg;
%%% G. Battery (omega)
eb = p_user(idx_eb); var.eb = eb;
eb0 = p_user(idx_eb0); var.eb0 = eb0;
pch = p_user(idx_pch); var.pch = pch;
pdis = p_user(idx_pdis); var.pdis = pdis;
%%% H. Aggregate Supply (omega)
pa = p_user(idx_pa); var.pa = pa;
%%% I. Operator Cost (omega)
po = p_user(idx_po); var.po = po;
pomax = p_user(idx_pomax); var.pomax = pomax;
z = p_user(idx_z); var.z = z;
zmax = p_user(idx_zmax); var.zmax = zmax;
%%% J. Investment (operator)
beta_S = p_user(idx_beta_S); var.beta_S = beta_S;
beta_W = p_user(idx_beta_W); var.beta_W = beta_W;
beta_B = p_user(idx_beta_B); var.beta_B = beta_B;
%%% Finalization

%% Calculate Cost
cost_pf = zeros(1, situation);
cost_pc = zeros(1, situation);
cost_pac = zeros(1, situation);
cost_po = zeros(1, situation);
for s = 1:situation
cost_pf(s) = D*gamma_pf*prob(s)*(pf((1+24*user*(s-1)):24*user*s) - Pf_ref(:))'*(pf((1+24*user*(s-1)):24*user*s) - Pf_ref(:));
cost_pc(s) = D*gamma_pc*prob(s)*(pc((1+24*user*(s-1)):24*user*s) - Pc_ref(:))'*(pc((1+24*user*(s-1)):24*user*s) - Pc_ref(:));
cost_pac(s) = D*gamma_pac*prob(s)*(tin((1+24*user*(s-1)):24*user*s) - ones(24*user,1)*Tin_ref)'*(tin((1+24*user*(s-1)):24*user*s) - ones(24*user,1)*Tin_ref);
cost_po(s) = D*gamma_pac*prob(s)*(pi1*sum(z((1+24*(s-1)):24*s)) + pi2*sum(zmax(s)));
end    
invest = c_S*beta_S + c_W*beta_W + c_B*beta_B;
cost = [cost_pf; cost_pc; cost_pac; cost_po];

