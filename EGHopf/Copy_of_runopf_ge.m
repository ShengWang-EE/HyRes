function [solution, information] = runopf_ge(mpc)
mpc0 = mpc;
%% initialization
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);

nGb = size(mpc.Gbus,1);
nGs = size(mpc.Gsou,1);
nGl = size(mpc.Gline,1);
iGd = find(mpc.Gbus(:,3)~=0);
nGd = size(iGd,1);
nb = size(mpc.bus,1);
ng = size(mpc.gen,1);
id = find(mpc.bus(:,3)~=0);
nd = size(id,1);
i_gpp = find(mpc.gen_extra.fuel=="Natural Gas");
n_gpp = size(i_gpp,1);
i_pipe = find(string(mpc.Gline_extra.Topology) == "Pipeline");
i_comp_reg = find(string(mpc.Gline_extra.Topology) ~= "Pipeline");
i_comp = find(contains(string(mpc.Gline_extra.Topology), "compressor", 'IgnoreCase', true) );
i_reg = find(contains(string(mpc.Gline_extra.Topology), "Pressure regulator",'IgnoreCase', true));
n_pipe = size(i_pipe,1);
n_comp = size(i_comp,1);

[GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF] = initializeParameters();
k_comp = 0.03;
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc);
upf = mpc.branch(il, RATE_A) - Pfinj(il);
mpc.Gbus(:,3) = mpc.Gbus(:,3) *3.6*1e6/GCV.ng;
mpc.Gsou(:,3) = mpc.Gsou(:,3) *3.6*1e6/GCV.ng;
mpc.Gsou(:,4) = mpc.Gsou(:,4) *3.6*1e6/GCV.ng;
%% ------------------------pre solve---------------------------------------
%% state variable
prs_sq = sdpvar(nGb,1); % nodal gas pressure
Qgs = sdpvar(nGs,1); % gas production of gas source
gf = sdpvar(nGl,1); % gas flow in the pipeline
Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)

% alpha = binvar(nGl,1); % direction of gas flow, 0-1        
% alpha(i_direction) = mpc.Gline_extra.Predirection(i_direction);
% gamma = (alpha-0.5)*2;          % 1,-1
delta = sdpvar(n_pipe,1);
LCe = 0;        
LCg = sdpvar(nGd,1);
% compressor
comp_ratio_sq = sdpvar(n_comp,1); % square of compression ratio

%% upper and lower bounds
prs_sq_min = mpc.Gbus(:,5).^2;
prs_sq_max = mpc.Gbus(:,6).^2;
gs_min = mpc.Gsou(:,3);              
gs_max = mpc.Gsou(:,4); % convert GWh/day to Mm3/day
gf_min = -mpc.Gline(:,5)/10;             gf_max = mpc.Gline(:,5)/10; 
if ~isempty(mpc.ptg)
    Qptg_min = mpc.ptg.Qptg_min;        Qptg_max = mpc.ptg.Qptg_max;
end
Va_min = -Inf(nb,1);                    Va_max = Inf(nb,1); 
Pg_min = mpc.gen(:, PMIN);                  Pg_max = mpc.gen(:, PMAX);
line_max = mpc.branch(il,RATE_A);
% compressor
comp_ratio_max = 2; % [1.5,2]都合理
comp_ratio_min = 1;
%% cons
% ptg cons
if ~isempty(mpc.ptg)
    n_ptg = size(mpc.ptg,1);
    Qptg = sdpvar(n_ptg,1); % gas (hydrogen) production of PTG
    Pptg = Qptg/24/3600 * GCV.ng * eta.electrolysis;
    ptg_cons = [    Qptg_min <= Qptg <= Qptg_max;];
else
    Pptg = 0;
    Qptg = 0;
end
% box cons
boxCons = [
    prs_sq_min <= prs_sq <= prs_sq_max;
    gs_min <= Qgs <= gs_max;
    Pg_min <= Pg <= Pg_max;
    LCg >= 0;
    ];
electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,LCe,mpc,id) == 0;];
electricityBranchFlowConsDC = [- line_max <= Bf(il,:)*Va <= line_max;];

% gas balance
Pgpp = Pg(i_gpp);
Qgpp = Pgpp / GCV.ng * 3600 * 24 * eta.GFU;
q_comp = k_comp .* gf(i_comp) .* (comp_ratio_sq-1); % gas consumption of comp, assume all gas driven
gasNodalBalanceCons = [consfcn_gasNodalBalance(gf,Qgs,Qgpp,LCg,Qptg,q_comp,mpc,iGd,i_gpp,i_comp,nGb,nGl,nGs,nGd,n_gpp) == 0;];

% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
C = mpc.Gline(:,3);
gasFlowCapacityCons = [
    - gf_max <= gf <= gf_max;
    ];
%
% compressor cons
fb_comp = FB(i_comp); tb_comp = TB(i_comp);
compCons = [
    comp_ratio_min.^2 <= comp_ratio_sq <= comp_ratio_max.^2;
    gf(i_comp) >= 0;
    prs_sq(fb_comp) .* comp_ratio_sq == prs_sq(tb_comp);
    ];
% pipe cons
fb_pipe = FB(i_pipe); tb_pipe = TB(i_pipe);
gasFlowCons = [
    gf(i_pipe).^2 <= C(i_pipe).^2 .* delta;
    delta >= prs_sq(fb_pipe) - prs_sq(tb_pipe);
    delta >= prs_sq(tb_pipe) - prs_sq(fb_pipe);
    0 <= delta <= (prs_sq_max(1) - prs_sq_min(1));
    ];
% regulator cons
fb_reg = FB(i_reg); tb_reg = TB(i_reg);
regCons = [
    gf(i_reg) >= 0;
    prs_sq(tb_reg) <= prs_sq(fb_reg);
    ];
%% solve
cons = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasFlowCapacityCons;
    gasFlowCons;
    compCons;
    regCons;
    ];
penalty1 = 1e6; penalty2 = 1e1;
objfcn = objfcn_IEGSoperatingCost(Pg,Qgs,mpc) + 1e-2*sum(gf.^2) + penalty1*sum(LCg) + penalty2*sum(delta) + 1e-4*sum((prs_sq-65^2).^2); % closer to mid level

opts1 = sdpsettings('solver','gurobi');
opts1.gurobi.NonConvex = 2;

information = optimize(cons,objfcn,opts1);

%% first results
prs_sq_val   = value(prs_sq);
Qgs_val   = value(Qgs);
gf_val    = value(gf);
Va_val    = value(Va);
Pg_val    = value(Pg);
delta_val = value(delta);
LCg_val   = value(LCg);
%% decide direction
% gap = C(i_pipe).^2 .* delta_val - gf_val(i_pipe).^2;
% 
% eps_f = 1e-6;
% r_gap = gap ./ max(gf_val(i_pipe).^2, eps_f);
% 
% tau_f = 0.03;
% tau_r = 0.10;
% 
% % 你现在用的是全网统一压力界 -> 每条管道同一个 sqrt(dpi)
% dpi_global = max(prs_sq_max(1) - prs_sq_min(1), 0);
% gf_max_numerical = abs(C) .* sqrt(dpi_global);
% 
% lockable = (abs(gf_val) >= tau_f .* gf_max_numerical) & (r_gap <= tau_r);
dir = sign(gf_val(i_pipe));                 % +1: FB->TB, -1: TB->FB, 0: near 0
% % 加上pressure regulator和compressor的方向
% lockable(i_comp_reg) = 1;
% dir(i_comp_reg) = 1;

% 只对 lockable 且 dir!=0 的边锁方向
% idx_pos = find(lockable & (dir > 0));
% idx_neg = find(lockable & (dir < 0));
idx_pos = find((dir >= 0));
idx_neg = find((dir < 0));

gasDirCons = [
    gf(i_pipe(idx_pos)) >= 0;
    gf(i_pipe(idx_neg)) <= 0;
];
%% second stage

% For dir = +1:  pi_FB - pi_TB = gf^2 / C^2
% For dir = -1:  pi_TB - pi_FB = gf^2 / C^2
% Also enforce consistent pressure ordering (helps solver).
gasWeymouthEqCons = [
    prs_sq(fb_pipe(idx_pos)) - prs_sq(tb_pipe(idx_pos)) == gf(i_pipe(idx_pos)).^2 ./ (C(i_pipe(idx_pos)).^2);
    prs_sq(fb_pipe(idx_pos)) >= prs_sq(tb_pipe(idx_pos));

    prs_sq(tb_pipe(idx_neg)) - prs_sq(fb_pipe(idx_neg)) == gf(i_pipe(idx_neg)).^2 ./ (C(i_pipe(idx_neg)).^2);
    prs_sq(tb_pipe(idx_neg)) >= prs_sq(fb_pipe(idx_neg));
];

% ----- 3) Build full constraint set for stage2 -----
cons2 = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasFlowCapacityCons;     % -gf_max <= gf <= gf_max (keep!)
    gasDirCons;
    gasWeymouthEqCons;
    compCons;
    regCons;
];

% Optional: you can drop penalty2*sum(delta) since delta is no longer used.
objfcn2 = objfcn_IEGSoperatingCost(Pg,Qgs,mpc) + 1e-2*sum(gf.^2) + penalty1*sum(LCg) + 1e-4*sum((prs_sq-65^2).^2);

% ----- 4) Warm start from stage1 solution (highly recommended) -----
assign(prs_sq, prs_sq_val);  % stage1 numeric pi = p^2
assign(Qgs,    Qgs_val);
assign(gf,     gf_val);
assign(Va,     Va_val);
assign(Pg,     Pg_val);
assign(LCg,    LCg_val);
if ~isempty(mpc.ptg)
    assign(Qptg, Qptg_val);
end

% ----- 5) Solve with nonconvex enabled -----
opts2 = sdpsettings('solver','gurobi');
opts2.usex0 = 1;
opts2.gurobi.NonConvex = 2;
% (optional) sometimes helps:
% opts3.gurobi.MIPFocus = 1;  % not MIP but sometimes ok to omit
% opts3.gurobi.NumericFocus = 3;

information2 = optimize(cons2, objfcn2, opts2);

% ----- 6) Collect results -----
prs2_val   = sqrt(value(prs_sq));
Qgs2_val   = value(Qgs);
gf2_val    = value(gf);
Va2_val    = value(Va);
Pg2_val    = value(Pg);
LCg2_val   = value(LCg);

% Pack
solution.stage2.prs = prs2_val;
solution.stage2.Qgs = Qgs2_val;
solution.stage2.gf  = gf2_val;
solution.stage2.Va  = Va2_val;
solution.stage2.Pg  = Pg2_val;
solution.stage2.LCg = LCg2_val;
dir_all = zeros(nGl,1);
dir_all(i_pipe) = sign(gf_val(i_pipe));
dir_all(i_comp) = 1;   % compressor fixed
dir_all(i_reg)  = 1;   % regulator fixed
solution.stage2.dir = dir_all;

information.stage2 = information2;
%% calculate linepack
% ---------- System linepack calculation (pipelines only) ----------
% Inputs you already have:
% prs_sq_val : nGb x 1  (pressure^2, in bar^2 if mpc.Gbus pressure is in bar)
% FB, TB     : nGl x 1  (from/to gas bus indices)
% i_pipe     : indices of pipelines within Gline
% Need: length (m) and diameter (m) for each pipeline in i_pipe

% ---- 1) get p (Pa) at nodes ----
prs_Pa  = prs2_val * 1e5;                    % Pa  (1 bar = 1e5 Pa)

% ---- 2) choose how to approximate pressure along a pipe ----
% simplest: average of end pressures (Pa)
p_avg = 0.5 * (prs_Pa(FB(i_pipe)) + prs_Pa(TB(i_pipe)));  % Pa

% ---- 3) pipe geometry ----
% Replace these two lines with your actual fields:
L = mpc.Gline_extra.Length(i_pipe) * 1000;       % m
D = mpc.Gline_extra.Diameter(i_pipe);     % m

A = pi/4 * D.^2;                            % m^2
Vpipe = A .* L;                             % m^3 (geometric volume)

% ---- 4) thermodynamics to convert to standard cubic meters (Sm^3) ----
R = 8.314462618;                            % J/(mol*K)
T  = 288.15;                                 % K (assumed gas temp)
Z  = 0.95;                                   % compressibility factor (assumed)

p0 = 1.01325e5;                              % Pa (standard)
T0 = 288.15;                                 % K (standard)

% moles in pipe at operating conditions: n = pV/(ZRT)
n_mol = (p_avg .* Vpipe) ./ (Z * R * T);     % mol

% convert moles to standard volume: Vstd = n*R*T0/p0
Vstd_Sm3 = n_mol .* R .* T0 ./ p0;           % Sm^3
Vstd_mmsm3 = Vstd_Sm3 / 1e6;                     % million standard m^3

% ---- 5) sum system linepack ----
LP_mmsm3 = sum(Vstd_mmsm3);                      % Sm^3


%% plot to check
% plot_UK_gas_network(mpc);
% % 标注压强
% GBlat = mpc.Gbus_extra.Lat; GBlon = mpc.Gbus_extra.Lon;
% for i = 1:size(mpc.Gbus,1)
%     text(GBlat(i), GBlon(i), sprintf('%.2f', prs2_val(i)), 'FontSize', 10, 'Color', 'k', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% end
end

