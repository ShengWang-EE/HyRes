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
n_miss_direction = size(mpc.Gline_extra.Predirection ==0,1); %认为非管道的元件都事先给定方向了
i_direction = find(mpc.Gline_extra.Predirection ~=0);

[GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF] = initializeParameters();
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc);
upf = mpc.branch(il, RATE_A) - Pfinj(il);
mpc.Gbus(:,3) = mpc.Gbus(:,3) *3.6*1e6/GCV.ng;
mpc.Gsou(:,3) = mpc.Gsou(:,3) *3.6*1e6/GCV.ng;
mpc.Gsou(:,4) = mpc.Gsou(:,4) *3.6*1e6/GCV.ng;
%% ------------------------pre solve---------------------------------------
%% state variable
prs = sdpvar(nGb,1); % nodal gas pressure
Qgs = sdpvar(nGs,1); % gas production of gas source
gf = sdpvar(nGl,1); % gas flow in the pipeline
Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)

LCe = 0;        LCg = 0;
%% upper and lower bounds
prs_min = mpc.Gbus(:,5);
prs_max = mpc.Gbus(:,6);
gs_min = mpc.Gsou(:,3);              
gs_max = mpc.Gsou(:,4); % convert GWh/day to Mm3/day
gf_min = -mpc.Gline(:,5)/10;             gf_max = mpc.Gline(:,5)/10; 
if ~isempty(mpc.ptg)
    Qptg_min = mpc.ptg.Qptg_min;        Qptg_max = mpc.ptg.Qptg_max;
end
Va_min = -Inf(nb,1);                    Va_max = Inf(nb,1); 
Pg_min = mpc.gen(:, PMIN);                  Pg_max = mpc.gen(:, PMAX);
line_max = mpc.branch(il,RATE_A);
%% contraints
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
%
boxCons_pre = [
    prs_min <= prs <= prs_max;
    gs_min <= Qgs <= gs_max;
    Pg_min <= Pg <= Pg_max;
    ];
electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,LCe,mpc,id) == 0;];
electricityBranchFlowConsDC = [- line_max <= Bf(il,:)*Va <= line_max;];
%
Pgpp = Pg(i_gpp);
Qgpp = Pgpp / GCV.ng * 3600 * 24 * eta.GFU;
gasNodalBalanceCons = [consfcn_gasNodalBalance(gf,Qgs,Qgpp,LCg,Qptg,mpc,iGd,i_gpp,nGb,nGl,nGs,nGd,n_gpp) == 0;];

%
gasPipelineFlowCons_pre = [
    -gf_max <= gf <= gf_max;
    % gf(i_direction) >= 0;
    ];

% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
C = mpc.Gline(:,3);
gasFlowCons_pre = [
        prs(FB) - prs(TB) == gf ./ (C*4);
        ];
%% solotion
cons_pre = [
    boxCons_pre;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasPipelineFlowCons_pre;
    gasFlowCons_pre;
    ];
objfcn = objfcn_IEGSoperatingCost(Pg,Qgs,mpc) + 1e-2*sum(gf.^2); % closer to mid level
yalmipOptions = sdpsettings('solver','gurobi');
yalmipOptions.gurobi.MIPgap = 1e-2;
yalmipOptions.gurobi.MIPFocus = 1;
yalmipOptions.gurobi.Heuristics = 0.2;
information_pre = optimize(cons_pre,objfcn,yalmipOptions);

gamma = sign(value(gf));
yalmip('clear');
%% -----------------solve----------------------------
%% state variable again
prs = sdpvar(nGb,1); % nodal gas pressure
Qgs = sdpvar(nGs,1); % gas production of gas source
gf = sdpvar(nGl,1); % gas flow in the pipeline
Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)
%% contraints again
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
%
boxCons = [
    prs_min <= prs <= prs_max;
    gs_min <= Qgs <= gs_max;
    Pg_min <= Pg <= Pg_max;
    ];
electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,LCe,mpc,id) == 0;];
electricityBranchFlowConsDC = [- line_max <= Bf(il,:)*Va <= line_max;];
%
Pgpp = Pg(i_gpp);
Qgpp = Pgpp / GCV.ng * 3600 * 24 * eta.GFU;
gasNodalBalanceCons = [consfcn_gasNodalBalance(gf,Qgs,Qgpp,LCg,Qptg,mpc,iGd,i_gpp,nGb,nGl,nGs,nGd,n_gpp) == 0;];

% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
C = mpc.Gline(:,3);
gasPipelineFlowCons = [
    (gamma-1) .* gf_max / 2 <= gf <= (gamma+1) .* gf_max / 2;
    ];
%
gasFlowCons = [
    prs(FB).^2 - prs(TB).^2 == gamma .* gf.^2 ./ C.^2;
    ];
%
cons = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasPipelineFlowCons;
    gasFlowCons;
    ];
objfcn = objfcn_IEGSoperatingCost(Pg,Qgs,mpc) + 1e-2*sum(gf.^2); % closer to mid level
information = optimize(cons,objfcn,yalmipOptions);

%% results
prs = value(prs); % nodal gas pressure
Qgs = value(Qgs); % gas production of gas source
gf = value(gf); % gas flow in the pipeline
Qptg = value(Qptg); % gas (hydrogen) production of PTG
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
gamma = value(gamma);
Qgpp = value(Qgpp);
[operatingCost,electricityGenerationCost,gasPurchasingCost] = objfcn_IEGSoperatingCost(Pg,Qgs,mpc);

[solution.prs, solution.PGs, solution.Gf, solution.Qptg, solution.Qgpp, ...
    solution.Va, solution.Pg, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost,  ...
    solution.gasPurchasingCost] = deal(...
    prs, Qgs, gf, Qptg, Qgpp, Va, Pg, gamma, ...
    operatingCost, electricityGenerationCost, gasPurchasingCost);
%% plot
plot_UK_gas_network(mpc);
% 标注压强
GBlat = mpc.Gbus_extra.Lat; GBlon = mpc.Gbus_extra.Lon;
for i = 1:size(mpc.Gbus,1)
    text(GBlat(i), GBlon(i), sprintf('%.2f', prs(i)), 'FontSize', 10, 'Color', 'k', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end
end
