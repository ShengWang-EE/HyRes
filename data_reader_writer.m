function mpc1 = data_reader_writer(option)
% 根据pypsa里面读取数据的文件的框架翻译过来
% pypsa里面的太复杂了，先做出个MVP，再添加细节吧
% 还是按照matpower调整吧，这样方便debug
% 设置多种模式，一种是输出纯电力系统，一种是输出电力天然气系统（主要是gencost）
%%
% clc
% clear

%% electricity demand
% preprocess_UK_power_demand_data(); % 将gridwatch的raw data预处理
gridwatch_data = readtable('data\power demand (from gridwatch)\power demand 2011-2025 (halfhour).csv');

% accroding to date range sort out data (valid year from 2012-2024)
start_date = '2024-1-1'; end_date = '2024-12-30';
selected_date = (gridwatch_data.timestamp >= start_date) & (gridwatch_data.timestamp <= end_date);
electricity_system_data = gridwatch_data(selected_date,:);
electricity_system_data.demand = fillmissing( ...
    electricity_system_data.demand, 'linear');
demandDistributionInfo = readtable('PyPSA-GB/data/demand/Demand_Distribution.csv');

electricity_nodal_demand = demandDistributionInfo.Var2(2:end)/sum(demandDistributionInfo.Var2(2:end)) * electricity_system_data.demand';

mean_electric_demand = mean(electricity_nodal_demand,2);

% powerDemand2024mean = 28420; % MW
% powerDemand2024peak = 47542;
% gasDemand2024mean = 196.8139; % mcm/day
% gasDemand2024peak = 387.5100; 
%% bus
% property needed: BusID	Type	Pd	Qd	Gs	Bs	ESystemID	Vm	Va	BaseKV	Zone	Vmax	Vmin
bus_ref = readtable('data/GB energy network v1.xlsx','Sheet','bus','Range','A1:M30');
busInfo = readtable('PyPSA-GB/data/network/BusesBasedGBsystem/buses.csv'); % 包含联络线，等效成发电机
[busInfo.name, busInfo.carrier] = deal(string(busInfo.name), string(busInfo.carrier));

nb = size(find(busInfo.carrier == "AC"),1);
% bus 
bus = array2table(zeros(nb,13), 'VariableNames',{'BusID','Type','Pd','Qd','Gs','Bs','ESystemID','Vm','Va','BaseKV','Zone','Vmax','Vmin'});
bus.BusID = [1:nb]';
bus.Type = bus_ref.Type;
bus.Pd = mean_electric_demand;
bus.Qd = zeros(nb,1);   % DC opf
bus.Gs = zeros(nb,1);
bus.Bs = zeros(nb,1);
bus.ESystemID = ones(nb,1);
bus.Vm = ones(nb,1);
bus.BaseKV = busInfo.v_nom(1:nb);
bus.Zone = ones(nb,1);
bus.Vmax = 1.1*ones(nb,1);
bus.Vmin = 0.9*ones(nb,1);
% bus extra
bus_extra.name = busInfo.name(1:nb); 
bus_extra.lon = busInfo.x(1:nb); bus_extra.lat = busInfo.y(1:nb);

bus_extra = struct2table(bus_extra);

mpc.bus = bus;                            
mpc.bus_extra = bus_extra;
%% branch/line
% property: FromBus	ToBus	R	X	b	RateA	RateB	RateC	K	Angle	Status	ang_min	ang_max
branch_ref = readtable('GB energy network v1.xlsx','Sheet','branch','Range','A1:M100');
branchInfo = readtable('PyPSA-GB/data/network/BusesBasedGBsystem/lines.csv');
n_branch = size(branchInfo,1);

branch = array2table(zeros(n_branch,13),'VariableNames',{'FromBus','ToBus','r','x','b','RateA','RateB','RateC','K','Angle','Status','ang_min','ang_max'});
[~,branch.FromBus] = ismember(branchInfo.bus0,bus_extra.name);
[~,branch.ToBus] = ismember(branchInfo.bus1,bus_extra.name);
[branch.r,branch.x,branch.b,branch.RateA,branch.RateB,branch.RateC] = ...
    deal(branchInfo.r,branchInfo.x,branchInfo.b,branchInfo.s_nom,branchInfo.s_nom,branchInfo.s_nom);
branch.K = branch_ref.K;
branch.Angle = branch_ref.Angle;
branch.Status = branch_ref.Status;
branch.ang_min = -360*ones(n_branch,1);
branch.ang_max = 360*ones(n_branch,1);

branch_extra.fbName = branchInfo.bus0; branch_extra.tbName = branchInfo.bus1;
mpc.branch = branch;
mpc.branch_extra = struct2table(branch_extra);
%% gen
% gen所在的bus按照距离最近判断
% property: BusID	Pg	Qg	Qmax	Qmin	Vg	mBase	Status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	RampRate	ramp_10	ramp_30	ramp_q	apf

% gen = readtable('GB energy network v1.xlsx','Sheet','gen','Range','A1:U67');
% gencost = readtable('GB energy network v1.xlsx','Sheet','gencost','Range','A1:G67');
genInfo = readtable('PyPSA-GB/data/power stations/power_stations_locations_2020.csv');
n_gen = size(genInfo,1);

gen = array2table(zeros(n_gen,21),'VariableNames',{'BusID','Pg','Qg','Qmax','Qmin','Vg','mBase','Status','Pmax', ...
    'Pmin','Pc1','Pc2','Qc1min','Qc1max','Qc2min','Qc2max','RampRate','ramp_10','ramp_30','ramp_q','apf'});

gen_extra.lon = genInfo.x; gen_extra.lat = genInfo.y;
for i = 1:size(genInfo,1)
    [~, gen.BusID(i,1)] = min(distance(gen_extra.lat(i),gen_extra.lon(i),bus_extra.lat,bus_extra.lon));
end

gen.Pmax = genInfo.InstalledCapacity_MW_; gen.Pmin =zeros(size(genInfo,1),1);
% ramp
genFuelInfo = readtable('PyPSA-GB/data/generator_data_by_fuel.csv');
fuelCostInfo = readtable('PyPSA-GB/data/marginal_cost_data.xlsx');
for i = 1:size(genInfo,1)
    iType = find(genInfo.Fuel(i) == string(genFuelInfo.fuel));
    if isempty(iType)
        iType = find(genInfo.Technology(i) == string(genFuelInfo.fuel));
    end
    gen_extra.min_up_time(i,1) = genFuelInfo.min_up_time(iType); gen_extra.min_down_time(i,1) = genFuelInfo.min_down_time(iType);
    gen_extra.ramp_up(i,1) = genFuelInfo.ramp_limit_up(iType); gen_extra.ramp_down(i,1) = genFuelInfo.ramp_limit_down(iType);
    % cost
    gencost_extra.start(i,1) = genFuelInfo.start_up_cost(iType);
    if string(genFuelInfo.fuel(iType)) == 'Coal'
        gencost_extra.marginal(i,1) = fuelCostInfo.Coal_p_kWh_(end) / 100 * 1e3; % pound/MWh
    elseif string(genFuelInfo.fuel(iType)) == 'Oil'
        gencost_extra.marginal(i,1) = fuelCostInfo.Oil_p_kWh_(end) / 100 * 1e3; % pound/MWh
    elseif ismember(string(genFuelInfo.fuel(iType)), ["CCGT", "OCGT"]) % gas fired units
        if option == 1 % power system data only
            gencost_extra.marginal(i,1) = fuelCostInfo.Gas_p_kWh_(end) / 100 * 1e3; % pound/MWh
        elseif option == 2 % power and gas system data
            gencost_extra.marginal(i,1) = 0; % calculated from gas system side
        end
    end
end
gen_extra.fuel = genInfo.Fuel; gen_extra.tech = genInfo.Technology;
gencost_extra = struct2table(gencost_extra);
% 补充其他gen property
gen.Vg = ones(n_gen,1);
gen.mBase = 100*ones(n_gen,1);
gen.Status = ones(n_gen,1);
mpc.gen_extra = gen_extra;
% gencost
gencost = array2table(zeros(n_gen,7),'VariableNames',{'CostType','StartCost','ShutCost','Order','CostA','CostB','CostC'});
gencost.CostType = 2*ones(n_gen,1);
gencost.StartCost = gencost_extra.start;
gencost.ShutCost = zeros(n_gen,1);
gencost.Order = 3*ones(n_gen,1);
gencost.CostA = zeros(n_gen,1);
gencost.CostB = gencost_extra.marginal;
gencost.CostC = zeros(n_gen,1);

%% renewable (tidal之类的考虑在未来场景里，现在没有）
% solar
solarInfo = readtable('\PyPSA-GB\data\renewables\atlite\inputs\Solar_Photovoltaics\Solar_Photovoltaics_2020.csv');
nSolar = size(solarInfo,1);
for i = 1:nSolar
    [~, solar.bus(i,1)] = min(distance(solarInfo.y(i),solarInfo.x(i),bus_extra.lat,bus_extra.lon));
end
solar.capacity = solarInfo.InstalledCapacity_MWelec_;
% onshore wind 
onshoreWindInfo = readtable('\PyPSA-GB\data\renewables\atlite\inputs\Wind_Onshore\Wind_Onshore_2020.csv');
nOnshoreWind = size(onshoreWindInfo,1);
for i = 1:nOnshoreWind
    [~, onshoreWind.bus(i,1)] = min(distance(onshoreWindInfo.y(i),onshoreWindInfo.x(i),bus_extra.lat,bus_extra.lon));
end
onshoreWind.capacity = onshoreWindInfo.TurbineCapacity_MW_ .* onshoreWindInfo.No_OfTurbines;
% offshore wind
offshoreWindInfo = readtable('\PyPSA-GB\data\renewables\atlite\inputs\Wind_Offshore\Wind_Offshore_2020.csv');
nOffshoreWind = size(offshoreWindInfo,1);
for i = 1:nOffshoreWind
    [~, offshoreWind.bus(i,1)] = min(distance(offshoreWindInfo.y(i),offshoreWindInfo.x(i),bus_extra.lat,bus_extra.lon));
end
offshoreWind.capacity = offshoreWindInfo.TurbineCapacity_MW_ .* offshoreWindInfo.No_OfTurbines;

% merge together to gen and gencost
renewable.lon = [solarInfo.x; onshoreWindInfo.x; offshoreWindInfo.x];
renewable.lat = [solarInfo.y; onshoreWindInfo.y; offshoreWindInfo.y];
renewable.bus = [solar.bus;onshoreWind.bus;offshoreWind.bus];
renewable.Pmax = [solar.capacity; onshoreWind.capacity; offshoreWind.capacity];
renewable.Pmin = zeros(nSolar+nOnshoreWind+nOffshoreWind,1);
renewable.min_up_time = zeros(nSolar+nOnshoreWind+nOffshoreWind,1);
renewable.min_down_time = zeros(nSolar+nOnshoreWind+nOffshoreWind,1);
renewable.ramp_up = renewable.Pmax;
renewable.ramp_down = renewable.Pmax;
renewable.fuel = [repmat("Solar",[nSolar,1]); repmat("Onshore wind",[nOnshoreWind,1]); repmat("Offshore wind",[nOffshoreWind,1]);];
renewable.tech = [repmat("Solar",[nSolar,1]); repmat("Onshore wind",[nOnshoreWind,1]); repmat("Offshore wind",[nOffshoreWind,1]);];
renewable = struct2table(renewable);

% update according to renewable profile
year = 2020;
solar_output_data1 = readtable("\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_"+year+"_1.csv");
solar_output_data2 = readtable("\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_"+year+"_2.csv");
solar_output_data3 = readtable("\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_"+year+"_3.csv");
solar_output_data4 = readtable("\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_"+year+"_4.csv");
solar_output_data = [solar_output_data1;solar_output_data2;solar_output_data3;solar_output_data4];
onshorewind_output_data = readtable("\PyPSA-GB\data\renewables\atlite\outputs\Wind_Onshore\Wind_Onshore_"+year+".csv");
offshorewind_output_data = readtable("\PyPSA-GB\data\renewables\atlite\outputs\Wind_Offshore\Wind_Offshore_"+year+".csv");

mean_renewable_capacity = [table2array(mean(solar_output_data(:,2:end),1)),table2array(mean(onshorewind_output_data(:,2:end),1)),table2array(mean(offshorewind_output_data(:,2:end),1))];
renewable.Pmax = mean_renewable_capacity';

% aggregate renewables
renewable_group = groupsummary(renewable, {'bus','fuel', 'tech'}, ...
    {'sum','max','mean'}, ...
    {'Pmax','Pmin','min_up_time','min_down_time','lon','lat','ramp_up','ramp_down'});
renewable_group.fuel = categorical(renewable_group.fuel,["Solar","Onshore wind","Offshore wind"], 'Ordinal',true);
renewable_group = sortrows(renewable_group,{'fuel','bus'});
renewable_agg_extra = table( ...
    renewable_group.mean_lon, renewable_group.mean_lat, renewable_group.bus, ...
    renewable_group.sum_Pmax, renewable_group.sum_Pmin, ...
    renewable_group.max_min_up_time, renewable_group.max_min_down_time, ...
    renewable_group.sum_ramp_up, renewable_group.sum_ramp_down, ...
    string(renewable_group.fuel),  renewable_group.tech, ...
    'VariableNames', {'lon','lat','bus','Pmax','Pmin','min_up_time','min_down_time','ramp_up','ramp_down','fuel','tech'} );

n_renewable_agg = size(renewable_group,1);
renewable_agg = array2table(zeros(n_renewable_agg,21), ...
    'VariableNames',{'BusID','Pg','Qg','Qmax','Qmin','Vg','mBase','Status','Pmax', ...
    'Pmin','Pc1','Pc2','Qc1min','Qc1max','Qc2min','Qc2max','RampRate','ramp_10','ramp_30','ramp_q','apf'});
renewable_agg.BusID = renewable_group.bus;
renewable_agg.Pmax = renewable_group.sum_Pmax;
renewable_agg.Pmin = renewable_group.sum_Pmin;
renewable_agg.Vg = ones(n_renewable_agg,1);
renewable_agg.mBase = 100*ones(n_renewable_agg,1);
renewable_agg.Status = ones(n_renewable_agg,1);

mpc.gen = [gen; renewable_agg];

renewable_agg_cost = array2table(zeros(n_renewable_agg,7),'VariableNames',{'CostType','StartCost','ShutCost','Order','CostA','CostB','CostC'});
renewable_agg_cost.CostType = 2*ones(n_renewable_agg,1);
renewable_agg_cost.StartCost = zeros(n_renewable_agg,1);
renewable_agg_cost.ShutCost = zeros(n_renewable_agg,1);
renewable_agg_cost.Order = 3*ones(n_renewable_agg,1);
renewable_agg_cost.CostA = zeros(n_renewable_agg,1);
renewable_agg_cost.CostB = zeros(n_renewable_agg,1);
renewable_agg_cost.CostC = zeros(n_renewable_agg,1);

mpc.gencost = [gencost; renewable_agg_cost]; % 补足聚合后的成本

%% interconnector/link (2020 scenario)
% interconnectorInfo = readtable('PyPSA-GB/data/interconnectors/links.csv');
% n_interconnector = size(interconnectorInfo,1);
% [~,interconnector.bus] = ismember(string(interconnectorInfo.bus1),string(busInfo.name));
% interconnector.capacity = interconnectorInfo.p_nom;
% mpc.interconnector = struct2table(interconnector);


%% --------------------------------------------gas system--------------------------------------------------
%% gas bus
Gbus = readtable('data/GB energy network v1.xlsx','sheet','Gbus','range','A1:I79');
mpc.Gbus = Gbus(:,1:6);
mpc.Gbus_extra = Gbus(:,7:end);

%% gas demand
filename = {'2022 gas demand.csv','2023 gas demand.csv','2024 gas demand.csv','2025 gas demand.csv'};
gas_demand_table = preprocess_UK_gas_demand_data(filename{1}); % unit, kWh/day
gas_demand_pure = (gas_demand_table.NTSEnergyOfftaken_IndustrialOfftakeTotal + gas_demand_table.NTSEnergyOfftaken_LDZOfftakeTotal ...
    + gas_demand_table.NTSEnergyOfftaken_StorageInjectionTotal)/1e6; % GPP excluded, GWh/day
gas_demand_mean = mean(gas_demand_pure);
gas_demand_distribution_factor = mpc.Gbus.Demand/100;
gas_demand_nodal = gas_demand_mean * gas_demand_distribution_factor; % distributed according to the percentage

mpc.Gbus.Demand = gas_demand_nodal;
%% pipeline
Gline = readtable('data/GB energy network v1.xlsx','sheet','Gline','range','A1:J91');
% 计算C (from chatgpt)
Gline.C = 2.85* Gline.Diameter.^(8/3) ./ sqrt(Gline.Length);
Gline.C(string(Gline.Topology) ~= "Pipeline") = 999;

mpc.Gline = Gline(:,1:5);
mpc.Gline_extra = Gline(:,6:end);
%% gas source
Gsou = readtable('data/GB energy network v1.xlsx','sheet','Gsou','range','A1:E10');
mpc.Gsou = Gsou(:,1:4);
mpc.Gsou_extra = Gsou(:,5);
%% gas cost
GcostInfo = readtable('data/GB energy network v1.xlsx','sheet','Gcost','range','A1:B214');
Gcost.Price = GcostInfo.Price(1);

mpc.Gcost = struct2table(Gcost);
%% gas storage
Gstore = readtable('data/GB energy network v1.xlsx','sheet','Gstore','range','A1:H10');
mpc.Gstore = Gstore;
%% ptg (none in base case)
mpc.ptg = [];
%% gas system location for gpp
ng = size(mpc.gen,1);
mpc.gen_extra.gas_bus = zeros(ng,1);
i_gpp = find(mpc.gen_extra.fuel =="Natural Gas");
n_gpp = size(i_gpp,1);

for i = 1:n_gpp
    i_gen = i_gpp(i);
    [~, mpc.gen_extra.gas_bus(i_gen)] = min(distance(gen_extra.lat(i_gen),gen_extra.lon(i_gen),Gbus.Lat,Gbus.Lon));
end
%% 转化为普通mpc
mpc.baseMVA = 100;
mpc1 = mpc;
mpc1.baseMVA = 100;
mpc1.version = '2';
mpc1.bus = mpc.bus{:,:};
mpc1.gen = mpc.gen{:,:};
mpc1.branch = mpc.branch{:,:};
mpc1.gencost = mpc.gencost{:,:};
mpc1.areas = [1,1];
mpc1.Gbus = mpc.Gbus{:,:};
mpc1.Gline = mpc.Gline{:,:};
mpc1.Gsou = mpc.Gsou{:,:};
mpc1.Gcost = mpc.Gcost{:,:};
mpc1.Gstore = mpc.Gstore{:,:};
% mpc1.ptg = mpc.ptg{:,:};
end

