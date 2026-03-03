function [historicalYearData] = readYearlyHistoricalData(iYear,mpc)
%% get large volume historical data and only keep for a specific year （先考虑光伏之类的，容量之后再考虑）

%% electricity demand
if iYear <= 2020 %都采用2017年的数据，只有17年是完整的
    systemDemandInfo = readtable('PyPSA-GB/data/demand/DemandData_2017.csv'); % 采用2017年的数据，只有17年是完整的
    demandDistributionInfo = readtable('PyPSA-GB/data/demand/Demand_Distribution.csv');
    electricityDemand = demandDistributionInfo.Var2(2:end)/sum(demandDistributionInfo.Var2(2:end)) * systemDemandInfo.ND';
    historicalYearData.electricityDemand = electricityDemand(:,1:2:end-1); % reduce the resolution to hourly
else
    error('unavaliable')
end

%% gas demand (英国没有小时分辨率的数据，因此用爱尔兰的参考）
% get uk gas demand (daily)
UKgasDemandInfo = readtable('Gas demand profile_01-01-2021_31-12-2021.csv','Range','A1:I366');

% 出力爱尔兰数据（包括补充一些缺失数据）
gasDemandCurve_NDM = table2array(readtable('Ireland gasconsumption 2023 merged.xlsx', ...
    'sheet','merged','range','D1:D8735'));
gasDemandCurve_LDM = table2array(readtable('Ireland gasconsumption 2023 merged.xlsx', ...
    'sheet','merged','range','J1:J8735'));
gasDemandCurve_power = table2array(readtable('Ireland gasconsumption 2023 merged.xlsx', ...
    'sheet','merged','range','P1:P8735'));
gasDemandTotal_raw = gasDemandCurve_NDM + gasDemandCurve_LDM + gasDemandCurve_power;
IEgasDemandTotal = zeros(8760,1);
IEgasDemandPowerProportion = zeros(8760,1);

dateTime_raw = convertCharsToStrings(table2cell(readtable('Ireland gasconsumption 2023 merged.xlsx', ...
    'sheet','merged','range','C1:C8735')));
for i = 1:size(dateTime_raw,1)
    dateTime_object{i} = datetime(dateTime_raw{i});
    monthValue(i) = month(dateTime_object{i});
    dayValue(i) = day(dateTime_object{i});
    hourValue(i) = hour(dateTime_object{i});
end
gasDemandPowerProportion_raw = gasDemandCurve_power ./ gasDemandTotal_raw;
for i = 1:size(gasDemandTotal_raw,1)
    differenceInHour = hours(dateTime_object{i} - dateTime_object{1}) + 1;
    IEgasDemandTotal(differenceInHour) = gasDemandTotal_raw(i);
    IEgasDemandPowerProportion(differenceInHour) = gasDemandPowerProportion_raw(i);
end
while min(IEgasDemandTotal) == 0
    for i = 1:8760
        if IEgasDemandTotal(i) == 0
            IEgasDemandTotal(i) = IEgasDemandTotal(i-1);
        end
        if IEgasDemandPowerProportion(i) == 0
            IEgasDemandPowerProportion(i) = IEgasDemandPowerProportion(i-1);
        end
    end
end

% 根据爱尔兰的负荷曲线修改UK的
nDay = 365;
gasDemand = zeros(24*nDay,1);
hourlyGasDemandForPowerRatio = zeros(24*nDay,1);
UKdailyGasDemandPowerProportion = UKgasDemandInfo.PowerstationsTotal ./ UKgasDemandInfo.DemandActualNTSD;

for iDay = 1:nDay
    startHour = (iDay-1) * 24 + 1;
    endHour = iDay * 24;
    % demand curve
    demandCurveShape = IEgasDemandTotal(startHour:endHour) / mean(IEgasDemandTotal(startHour:endHour));
    hourlyGasDemand = UKgasDemandInfo.DemandActualNTSD(iDay) * demandCurveShape;
    gasDemand(startHour:endHour) = hourlyGasDemand;
    % gas demand for power generation ratio curve
    ratioCurveShape = IEgasDemandPowerProportion(startHour:endHour) / mean(IEgasDemandPowerProportion(startHour:endHour));
    GasDemandPowerProportion = UKdailyGasDemandPowerProportion(iDay) * ratioCurveShape;
    hourlyGasDemandForPowerRatio(startHour:endHour) = GasDemandPowerProportion;
    
end
gasDemand_power = gasDemand .* hourlyGasDemandForPowerRatio;
gasDemand_nonpower = gasDemand .* (1-hourlyGasDemandForPowerRatio);
% 将系统负荷分解至节点
gasDistributionFactor = mpc.Gbus.demand;
nodalGasDemand_nonpower = gasDistributionFactor * gasDemand_nonpower';

historicalYearData.gasDemand = nodalGasDemand_nonpower;
%% renewables
if iYear <= 2020 %用19年的避免闰年
    % solar
    solarInfo1 = readtable('\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_2019_1.csv');
    solarInfo2 = readtable('\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_2019_2.csv');
    solarInfo3 = readtable('\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_2019_3.csv');
    solarInfo4 = readtable('\PyPSA-GB\data\renewables\atlite\outputs\PV\PV_2019_4.csv');
    % merge together
    solarCapacity = [solarInfo1(:,2:end);solarInfo2(:,2:end);solarInfo3(:,2:end);solarInfo4(:,2:end);];
    % onshore wind
    onshoreWindInfo = readtable('\PyPSA-GB\data\renewables\atlite\outputs\Wind_Onshore\Wind_Onshore_2019.csv');
    onshoreWindCapacity = onshoreWindInfo(:,2:end);
    % offshore wind
    offshoreWindInfo = readtable('\PyPSA-GB\data\renewables\atlite\outputs\Wind_Offshore\Wind_Offshore_2019.csv');
    offshoreWindCapacity = offshoreWindInfo(:,2:end);
end
renewableCapacity = [table2array(solarCapacity),table2array(onshoreWindCapacity),table2array(offshoreWindCapacity)]';
renewableCapacity(isnan(renewableCapacity)) = 0; %有些地方是空白的，就设为0吧
% 把海量renewable按照节点合并
nb = size(mpc.bus,1);
renewableCapacityInBus = zeros(nb,8760);
for i = 1:size(renewableCapacity,1)
    busIndex = mpc.renewable.bus(i);
    renewableCapacityInBus(busIndex,:) = renewableCapacityInBus(busIndex,:) + renewableCapacity(i,:);
end
historicalYearData.renewableCapacity = renewableCapacityInBus;
