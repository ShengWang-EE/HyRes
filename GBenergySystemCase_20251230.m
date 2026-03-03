% function mpc = GBenergySystemCase()
% 基础案例的GB电力系统数据来源于.m文件（数据也写入了GB energy network v1.xlsx文件里）；
% 那些PSA的csv，excel文件用于之后添加坐标名字等额外信息，以及构建不同场景的时候使用
% 或者不用.m文件里的，直接用csv里的，这样之后转换成能和原来的python代码更加对应的上一点
% 构建目前的network，再在这基础上构建未来的
% .m文件和csv表格里的数据不一致，比如gen，还是同意采用csv的吧
% 根据pypsa里面读取数据的文件的框架翻译过来
%%
clc
clear
mpc0 = GBreducednetwork;
%% bus
% property needed: 
bus = readtable('GB energy network v1.xlsx','Sheet','bus','Range','A1:M30');
busInfo = readtable('PyPSA-GB/data/network/BusesBasedGBsystem/buses.csv'); % 包含联络线，等效成发电机

busExtra.name = busInfo.name; 
busExtra.lon = busInfo.x; busExtra.lat = busInfo.y;
busExtra = struct2table(busExtra);

nb = 29;                                                                    % 29 bus
mpc.bus = [bus, busExtra(1:29,:)];                            % exclude interconnector here
%% branch/line
branch = readtable('GB energy network v1.xlsx','Sheet','branch','Range','A1:M100');
branchInfo = readtable('PyPSA-GB/data/network/BusesBasedGBsystem/lines.csv');

lineExtra.fbName = branchInfo.bus0; lineExtra.tbName = branchInfo.bus1;
[~,lineExtra.fb] = ismember(lineExtra.fbName,busExtra.name);
[~,lineExtra.tb] = ismember(lineExtra.tbName,busExtra.name);
[lineExtra.r,lineExtra.x,lineExtra.b,lineExtra.capacity] = deal(branchInfo.r,branchInfo.x,branchInfo.b,branchInfo.s_nom);

mpc.branch = branch;
% mpc.branchExtra = struct2table(lineExtra);
%% gen
% gen所在的bus不能按照最近距离算。看看手册
% gen = readtable('GB energy network v1.xlsx','Sheet','gen','Range','A1:U67');
% gencost = readtable('GB energy network v1.xlsx','Sheet','gencost','Range','A1:G67');
genInfo = readtable('PyPSA-GB/data/power stations/power_stations_locations_2020.csv');
n_gen = size(genInfo,1);



genExtra.lon = genInfo.x; genExtra.lat = genInfo.y;
for i = 1:size(genInfo,1)
    [~, genExtra.bus(i,1)] = min(distance(genExtra.lat(i),genExtra.lon(i),busExtra.lat,busExtra.lon));
end
genExtra.Pmax = genInfo.InstalledCapacity_MW_; genExtra.Pmin =zeros(size(genInfo,1),1);
% ramp
genFuelInfo = readtable('PyPSA-GB/data/generator_data_by_fuel.csv');
fuelCostInfo = readtable('PyPSA-GB/data/marginal_cost_data');
for i = 1:size(genInfo,1)
    iType = find(genInfo.Fuel(i) == string(genFuelInfo.fuel));
    if isempty(iType)
        iType = find(genInfo.Technology(i) == string(genFuelInfo.fuel));
    end
    genExtra.min_up_time(i,1) = genFuelInfo.min_up_time(iType); genExtra.min_down_time(i,1) = genFuelInfo.min_down_time(iType);
    genExtra.rampUp(i,1) = genFuelInfo.ramp_limit_up(iType); genExtra.rampDown(i,1) = genFuelInfo.ramp_limit_down(iType);
    % cost
    gencostExtra.start(i,1) = genFuelInfo.start_up_cost(iType);
    if string(genFuelInfo.fuel(iType)) == 'Coal'
        gencostExtra.marginal(i,1) = fuelCostInfo.Coal_p_kWh_(end) / 100 * 1e3; % pound/MWh
    elseif string(genFuelInfo.fuel(iType)) == 'Oil'
        gencostExtra.marginal(i,1) = fuelCostInfo.Oil_p_kWh_(end) / 100 * 1e3; % pound/MWh
    elseif ismember(string(genFuelInfo.fuel(iType)), ["CCGT", "OCGT"]) % gas fired units
        % gencost.marginal(i,1) = fuelCostInfo.Gas_p_kWh_(end) / 100 * 1e3; % pound/MWh
        gencostExtra.marginal(i,1) = 0; % calculated from gas system side
    end
end
genExtra.fuel = genInfo.Fuel; genExtra.tech = genInfo.Technology;
genExtra = struct2table(genExtra);
mpc.gen = gen;
mpc.gencost = gencost;
mpc.genExtra = genExtra;
mpc.gencostExtra = struct2table(gencostExtra);

%% distribution 层面的发电机考虑吗？有数据
%% renewable (tidal之类的考虑在未来场景里，现在没有）
% solar
solarInfo = readtable('\PyPSA-GB\data\renewables\atlite\inputs\Solar_Photovoltaics\Solar_Photovoltaics_2020.csv');
nSolar = size(solarInfo,1);
for i = 1:nSolar
    [~, solar.bus(i,1)] = min(distance(solarInfo.y(i),solarInfo.x(i),busExtra.lat,busExtra.lon));
end
solar.capacity = solarInfo.InstalledCapacity_MWelec_;
% onshore wind 
onshoreWindInfo = readtable('\PyPSA-GB\data\renewables\atlite\inputs\Wind_Onshore\Wind_Onshore_2020.csv');
nOnshoreWind = size(onshoreWindInfo,1);
for i = 1:nOnshoreWind
    [~, onshoreWind.bus(i,1)] = min(distance(onshoreWindInfo.y(i),onshoreWindInfo.x(i),busExtra.lat,busExtra.lon));
end
onshoreWind.capacity = onshoreWindInfo.TurbineCapacity_MW_ .* onshoreWindInfo.No_OfTurbines;
% offshore wind
offshoreWindInfo = readtable('\PyPSA-GB\data\renewables\atlite\inputs\Wind_Offshore\Wind_Offshore_2020.csv');
nOffshoreWind = size(offshoreWindInfo,1);
for i = 1:nOffshoreWind
    [~, offshoreWind.bus(i,1)] = min(distance(offshoreWindInfo.y(i),offshoreWindInfo.x(i),busExtra.lat,busExtra.lon));
end
offshoreWind.capacity = offshoreWindInfo.TurbineCapacity_MW_ .* offshoreWindInfo.No_OfTurbines;

% merge together
renewable.bus = [solar.bus;onshoreWind.bus;offshoreWind.bus];
renewable.capacity = [solar.capacity; onshoreWind.capacity; offshoreWind.capacity];
mpc.renewable = struct2table(renewable);
%% renewable profile
%% interconnector/link (2020 scenario)
interconnectorInfo = readtable('PyPSA-GB/data/interconnectors/links.csv');
n_interconnector = size(interconnectorInfo,1);
[~,interconnector.bus] = ismember(string(interconnectorInfo.bus1),string(busInfo.name));
interconnector.capacity = interconnectorInfo.p_nom;
mpc.interconnector = struct2table(interconnector);

%% electricity demand
systemDemandInfo = readtable('PyPSA-GB/data/demand/DemandData_2017.csv'); % 采用2017年的数据，只有17年是完整的
demandDistributionInfo = readtable('PyPSA-GB/data/demand/Demand_Distribution.csv');

electricityDemand2020 = demandDistributionInfo.Var2(2:end)/sum(demandDistributionInfo.Var2(2:end)) * systemDemandInfo.ND';

mean_electric_demand = mean(electricityDemand2020,2);
mpc.bus.demand = mean_electric_demand;

% powerDemand2024mean = 28420; % MW
% powerDemand2024peak = 47542;
% gasDemand2024mean = 196.8139; % mcm/day
% gasDemand2024peak = 387.5100; 
%% --------------------------------------------gas system--------------------------------------------------
%% gas source
Gbus = readtable('data/GB energy network v1.xlsx','sheet','Gbus','range','A1:I79');
mpc.Gbus = Gbus;
%% pipeline
Gline = readtable('data/GB energy network v1.xlsx','sheet','Gline','range','A1:J91');
mpc.Gline = Gline;
%% gas source
Gsou = readtable('data/GB energy network v1.xlsx','sheet','Gsou','range','A1:E10');
mpc.Gsou = Gsou;
%% gas cost
GcostInfo = readtable('data/GB energy network v1.xlsx','sheet','Gcost','range','A1:B214');
Gcost.Price = GcostInfo.Price(1);

mpc.Gcost = struct2table(Gcost);
%% gas storage
Gstore = readtable('data/GB energy network v1.xlsx','sheet','Gstore','range','A1:H10');
mpc.Gstore = Gstore;
%% ptg (none in base case)
%% gas system location for gpp

ng = size(mpc.genExtra,1);
for i = 1:ng
    [~, mpc.genExtra.gasbus(i,1)] = min(distance(genExtra.lat(i),genExtra.lon(i),Gbus.Lat,Gbus.Lon));
end
