function [solution, info] = runopf_ge(mpc)
%RUNOPF_GE  Integrated electricity-gas OPF (DC power flow + gas network)
% Two-stage solve:
%   Stage 1: convex relaxation (delta-based)
%   Stage 2: fix pipe directions + nonconvex Weymouth equalities (warm start)

%% --------------------------- indices & sets -----------------------------
% MATPOWER indices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
 VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; %#ok<ASGLU>
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
 MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
 QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen; %#ok<ASGLU>
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
 TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
 ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch; %#ok<ASGLU>

% active electric branches for DC thermal limits
i_line_lim = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);

% dimensions
n_gb = size(mpc.Gbus, 1);
n_gs = size(mpc.Gsou, 1);
n_gl = size(mpc.Gline, 1);
n_b  = size(mpc.bus,  1);
n_g  = size(mpc.gen,  1);

% gas demand buses (your convention: Gbus(:,3) is demand)
i_gas_dem_bus = find(mpc.Gbus(:, 3) ~= 0);
n_gas_dem     = numel(i_gas_dem_bus);

% electric load buses (your convention: bus(:,3)=PD)
i_elec_load_bus = find(mpc.bus(:, 3) ~= 0);

% gas-fired generators
i_gas_fired_gen = find(mpc.gen_extra.fuel == "Natural Gas");
n_gas_fired_gen = numel(i_gas_fired_gen);

% gas components by topology
topo = string(mpc.Gline_extra.Topology);
i_pipe = find(topo == "Pipeline");
i_comp = find(contains(topo, "compressor", "IgnoreCase", true));
i_reg  = find(contains(topo, "Pressure regulator", "IgnoreCase", true));

n_pipe = numel(i_pipe);
n_comp = numel(i_comp);

% gas link endpoints & coefficients
fb = mpc.Gline(:, 1);
tb = mpc.Gline(:, 2);
C  = mpc.Gline(:, 3);

%% --------------------------- parameters ---------------------------------
[GCV, ~, ~, ~, ~, ~, ~, ~, ~, eta, ~] = initializeParameters();

k_comp = 0.03;          % compressor gas consumption coefficient
comp_ratio_min = 1.0;
comp_ratio_max = 2.0;

% DC power flow matrices
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc); %#ok<ASGLU>
line_max = mpc.branch(i_line_lim, RATE_A);

% unit conversion (keep your original scaling)
mpc.Gbus(:, 3) = mpc.Gbus(:, 3) * 3.6e6 / GCV.ng;
mpc.Gsou(:, 3) = mpc.Gsou(:, 3) * 3.6e6 / GCV.ng;
mpc.Gsou(:, 4) = mpc.Gsou(:, 4) * 3.6e6 / GCV.ng;

%% --------------------------- decision vars ------------------------------
% gas
pi_sq   = sdpvar(n_gb, 1);            % squared pressure (p^2)
q_src   = sdpvar(n_gs, 1);            % source injection
q_flow  = sdpvar(n_gl, 1);            % link flow
lc_gas  = sdpvar(n_gas_dem, 1);       % gas load curtailment at demand nodes

% electricity (DC)
va      = sdpvar(n_b, 1);             % voltage angles
p_gen   = sdpvar(n_g, 1);             % generator outputs
lc_elec = 0;                          % keep your original (scalar 0)

% pipe relaxation variable
delta   = sdpvar(n_pipe, 1);

% compressor
comp_ratio_sq = sdpvar(n_comp, 1);    % squared compression ratio

% PTG (optional)
has_ptg = isfield(mpc, "ptg") && ~isempty(mpc.ptg);
if has_ptg
    n_ptg = size(mpc.ptg, 1);
    q_ptg = sdpvar(n_ptg, 1);
    p_ptg = q_ptg/24/3600 * GCV.ng * eta.electrolysis;
else
    q_ptg = 0;
    p_ptg = 0;
end

%% --------------------------- bounds -------------------------------------
pi_sq_min = mpc.Gbus(:, 5).^2;
pi_sq_max = mpc.Gbus(:, 6).^2;

q_src_min = mpc.Gsou(:, 3);
q_src_max = mpc.Gsou(:, 4);

q_flow_max =  mpc.Gline(:, 5)/10;
q_flow_min = -mpc.Gline(:, 5)/10;

p_gen_min = mpc.gen(:, PMIN);
p_gen_max = mpc.gen(:, PMAX);

if has_ptg
    q_ptg_min = mpc.ptg.Qptg_min;
    q_ptg_max = mpc.ptg.Qptg_max;
end

%% --------------------------- common constraints -------------------------
cons_box = [
    pi_sq_min <= pi_sq <= pi_sq_max
    q_src_min <= q_src <= q_src_max
    p_gen_min <= p_gen <= p_gen_max
    lc_gas >= 0
];

if has_ptg
    cons_ptg = [q_ptg_min <= q_ptg <= q_ptg_max];
else
    cons_ptg = [];
end

% electricity constraints (DC)
cons_elec_balance = [consfcn_electricityNodalBalance(va, p_gen, p_ptg, lc_elec, mpc, i_elec_load_bus) == 0];
cons_elec_branch  = [-line_max <= Bf(i_line_lim, :) * va <= line_max];

% gas balance
p_gas_fired = p_gen(i_gas_fired_gen);
q_gas_fired = p_gas_fired * 3600 * 24 / eta.GFU / GCV.ng; % MW to Mm3/day

q_comp = k_comp .* q_flow(i_comp) .* (comp_ratio_sq - 1);

cons_gas_balance = [
    consfcn_gasNodalBalance( ...
        q_flow, q_src, q_gas_fired, lc_gas, q_ptg, q_comp, ...
        mpc, i_gas_dem_bus, i_gas_fired_gen, i_comp, ...
        n_gb, n_gl, n_gs, n_gas_dem, n_gas_fired_gen) == 0
];

% gas flow bounds
cons_gas_cap = [q_flow_min <= q_flow <= q_flow_max];

%% --------------------------- compressor & regulator ----------------------
fb_comp = fb(i_comp);
tb_comp = tb(i_comp);
cons_comp = [
    comp_ratio_min^2 <= comp_ratio_sq <= comp_ratio_max^2
    q_flow(i_comp) >= 0
    pi_sq(fb_comp) .* comp_ratio_sq == pi_sq(tb_comp)
];

fb_reg = fb(i_reg);
tb_reg = tb(i_reg);
cons_reg = [
    q_flow(i_reg) >= 0
    pi_sq(tb_reg) <= pi_sq(fb_reg)
];

%% ============================= stage 1 ==================================
fb_pipe = fb(i_pipe);
tb_pipe = tb(i_pipe);

delta_max = max(pi_sq_max(1) - pi_sq_min(1), 0); % keep your global bound style

cons_pipe_relax = [
    q_flow(i_pipe).^2 <= (C(i_pipe).^2) .* delta
    delta >= pi_sq(fb_pipe) - pi_sq(tb_pipe)
    delta >= pi_sq(tb_pipe) - pi_sq(fb_pipe)
    0 <= delta <= delta_max
];

cons_stage1 = [
    cons_box
    cons_ptg
    cons_elec_balance
    cons_elec_branch
    cons_gas_balance
    cons_gas_cap
    cons_pipe_relax
    cons_comp
    cons_reg
];

pen_lc_gas = 1e6;
pen_delta  = 1e1;

obj_stage1 = objfcn_IEGSoperatingCost(p_gen, q_src, mpc) ...
           + 1e-2 * sum(q_flow.^2) ...
           + pen_lc_gas * sum(lc_gas) ...
           + pen_delta  * sum(delta) ...
           + 1e-4 * sum((pi_sq - 65^2).^2);

opts_stage1 = sdpsettings('solver', 'gurobi');
opts_stage1.gurobi.NonConvex = 2;

info.stage1 = optimize(cons_stage1, obj_stage1, opts_stage1);

% stage1 values (warm start)
pi_sq_1  = value(pi_sq);
q_src_1  = value(q_src);
q_flow_1 = value(q_flow);
va_1     = value(va);
p_gen_1  = value(p_gen);
lc_gas_1 = value(lc_gas);

if has_ptg
    q_ptg_1 = value(q_ptg); % <-- fix: needed for warm-start
end

%% ============================= stage 2 ==================================
% direction from stage1 (pipes only)
dir_pipe = sign(q_flow_1(i_pipe)); % +: fb->tb, -: tb->fb, 0: ~0

idx_pos = find(dir_pipe >= 0);
idx_neg = find(dir_pipe < 0);

cons_gas_dir = [
    q_flow(i_pipe(idx_pos)) >= 0
    q_flow(i_pipe(idx_neg)) <= 0
];

% Weymouth equalities (direction-consistent)
cons_weymouth_eq = [
    pi_sq(fb_pipe(idx_pos)) - pi_sq(tb_pipe(idx_pos)) == (q_flow(i_pipe(idx_pos)).^2) ./ (C(i_pipe(idx_pos)).^2)
    pi_sq(fb_pipe(idx_pos)) >= pi_sq(tb_pipe(idx_pos))

    pi_sq(tb_pipe(idx_neg)) - pi_sq(fb_pipe(idx_neg)) == (q_flow(i_pipe(idx_neg)).^2) ./ (C(i_pipe(idx_neg)).^2)
    pi_sq(tb_pipe(idx_neg)) >= pi_sq(fb_pipe(idx_neg))
];

cons_stage2 = [
    cons_box
    cons_ptg
    cons_elec_balance
    cons_elec_branch
    cons_gas_balance
    cons_gas_cap
    cons_gas_dir
    cons_weymouth_eq
    cons_comp
    cons_reg
];

obj_stage2 = objfcn_IEGSoperatingCost(p_gen, q_src, mpc) ...
           + 1e-2 * sum(q_flow.^2) ...
           + pen_lc_gas * sum(lc_gas) ...
           + 1e-4 * sum((pi_sq - 65^2).^2);

% warm start
assign(pi_sq,  pi_sq_1);
assign(q_src,  q_src_1);
assign(q_flow, q_flow_1);
assign(va,     va_1);
assign(p_gen,  p_gen_1);
assign(lc_gas, lc_gas_1);
if has_ptg
    assign(q_ptg, q_ptg_1);
end

opts_stage2 = sdpsettings('solver', 'gurobi', 'usex0', 1);
opts_stage2.gurobi.NonConvex = 2;

info.stage2 = optimize(cons_stage2, obj_stage2, opts_stage2);

%% --------------------------- collect results ----------------------------
pi_sq_2  = value(pi_sq);
q_src_2  = value(q_src);
q_flow_2 = value(q_flow);
va_2     = value(va);
p_gen_2  = value(p_gen);
lc_gas_2 = value(lc_gas);

p_bar_2 = sqrt(pi_sq_2);

% direction for all gas links (pipes from stage1 sign; comp/reg fixed +1)
dir_all = zeros(n_gl, 1);
dir_all(i_pipe) = sign(q_flow_1(i_pipe));
dir_all(i_comp) = 1;
dir_all(i_reg)  = 1;

% linepack (pipelines only)
linepack_mmsm3 = compute_linepack_mmsm3(mpc, p_bar_2, i_pipe, fb, tb);

solution.stage2.pressure_bar   = p_bar_2;
solution.stage2.pressure_sq    = pi_sq_2;
solution.stage2.q_src          = q_src_2;
solution.stage2.q_flow         = q_flow_2;
solution.stage2.va             = va_2;
solution.stage2.p_gen          = p_gen_2;
solution.stage2.lc_gas         = lc_gas_2;
solution.stage2.dir            = dir_all;
solution.stage2.linepack_mmsm3 = linepack_mmsm3;

end

%% ========================= local helper function =========================
function linepack_mmsm3 = compute_linepack_mmsm3(mpc, pressure_bar, i_pipe, fb, tb)
%COMPUTE_LINEPACK_MMSM3  System linepack (pipelines only), in million Sm^3

% node pressure (Pa)
p_pa = pressure_bar * 1e5;

% average pressure per pipe (Pa)
p_avg = 0.5 * (p_pa(fb(i_pipe)) + p_pa(tb(i_pipe)));

% geometry
L = mpc.Gline_extra.Length(i_pipe) * 1000; % km -> m (keep your convention)
D = mpc.Gline_extra.Diameter(i_pipe);      % m
A = (pi/4) * D.^2;
V = A .* L;                                 % m^3

% thermo (keep your constants)
R  = 8.314462618;  % J/(mol*K)
T  = 288.15;       % K
Z  = 0.95;         % -
p0 = 1.01325e5;    % Pa
T0 = 288.15;       % K

n_mol = (p_avg .* V) ./ (Z * R * T);
Vstd_sm3 = n_mol .* R .* T0 ./ p0;

linepack_mmsm3 = sum(Vstd_sm3) / 1e6;
end
