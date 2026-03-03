function [newmpc] = getProfile(mpc,iYear,iDay,historicalYearData)
% 这个函数用来确定时间相关的参数，包括长期的（容量扩展，发电机投入退役），短期的
% （负荷曲线，可再生能源）等，并根据这些参数修改mpc
%% ----------------------- long term -----------------------------------
%% power line capacity

%% gas pipeline capacity

%% fossil fuel generators

%% ------------------------short term ----------------------------------
%% electricity demand









baseMVA = 100;
[startHour,endHour,startDay,endDay] = selectScenario(options.season,mpc.electricityDemandCuve,onshoreWindAvaliableCapacity);
nDay = endDay - startDay + 1;
%% add offshore wind
switch options.year
    case 'now'
        OWFcapacity = 0;
    case '2030'
        OWFcapacity = 5; 
    case '2040'
        OWFcapacity = 20;
    case '2050'
        OWFcapacity = 37;
end
unitCapacity = LCOEcurve(2,1);
busCord = table2array(mpc.busName(:,3:4));
GbusCord = table2array(mpc.GbusName(:,2:3));
addGen = []; addPTG = [];
if OWFcapacity ~= 0
    [~,nOWF] = max(LCOEcurve(:,1) > OWFcapacity * 1e3);
    for i = 1:nOWF
        lat = latColumn(i); lon = lonColumn(i);
        distanceToBusE = (lon-busCord(:,1)).^2 + (lat-busCord(:,2)).^2;
        distanceToBusG = (lon-GbusCord(:,1)).^2 + (lat-GbusCord(:,2)).^2;
        [~,OWFbusE(i)] = min(distanceToBusE);
        [~,OWFbusG(i)] = min(distanceToBusG);
        addGen = [addGen;[OWFbusE(i)	unitCapacity	0	unitCapacity	-unitCapacity	1	100	1	unitCapacity	0	0	0	0	0	0	0	unitCapacity	0	0	0	0]];
        ptgUnitCapacity = OWFcapacity / nOWF * 0.5 * 0.2825 * 24;
        addPTG = [addPTG; [OWFbusG(i)	OWFbusE(i)	0.950000000000000	0	ptgUnitCapacity*2]];
    end
end
addGencost = repmat([2	0	0	3	0	0	0],[nOWF,1]);
addGentype = repmat(["Offshore"],[nOWF,1]);
mpc.gen = [mpc.gen; addGen]; mpc.ptg = [mpc.ptg; addPTG]; mpc.gencost = [mpc.gencost; addGencost]; mpc.genType = [mpc.genType; addGentype];
%% modify power/gas demand and transmission capacity (参考J13)
switch options.year
    case '2030'
        % mpc.bus(:,3:4) = mpc.bus(:,3:4) * 1.022^(2030-2023);
        mpc.electricityDemandCuve = mpc.electricityDemandCuve * 1.022^(2030-2023);
        mpc.branch(:,6:8) = mpc.branch(:,6:8) * 1.022^(2030-2023);
        mpc.Gbus(:,3) = mpc.Gbus(:,3) * 1.162805205;
        mpc.Gline(:,3) = mpc.Gline(:,3) * 1.162805205;
    case '2040'
        % mpc.bus(:,3:4) = mpc.bus(:,3:4) * 1.022^(2040-2023);
        mpc.electricityDemandCuve = mpc.electricityDemandCuve * 1.022^(2040-2023);
        mpc.branch(:,6:8) = mpc.branch(:,6:8) * 1.022^(2040-2023);
        mpc.Gbus(:,3) = mpc.Gbus(:,3) * 0.970086967;
        mpc.Gline(:,3) = mpc.Gline(:,3) * 0.970086967;
    case '2050'
        % mpc.bus(:,3:4) = mpc.bus(:,3:4) * 1.022^(2050-2023);
        mpc.electricityDemandCuve = mpc.electricityDemandCuve * 1.022^(2050-2023);
        mpc.branch(:,6:8) = mpc.branch(:,6:8) * 1.022^(2050-2023);
        mpc.Gbus(:,3) = mpc.Gbus(:,3) * 0.822776129;
        mpc.Gline(:,3) = mpc.Gline(:,3) * 0.822776129;
end

%% select curve according to days
NK = 24;
for iDay = 1:nDay
    yalmip('clear');
    periodHour = startHour+(iDay-1)*24:startHour+iDay*24-1;
    onshoreWindCapacity = onshoreWindAvaliableCapacity(periodHour,:);
    solarCapacity = onshoreSolarAvaliableCapacity(periodHour,:);
    hydroCapacity = hydroAvaliableCapacity(periodHour,:);
    interconnectorCapacity = interconnectorAvaliableCapacity(periodHour,:);
    gasDemandCurve = mpc.gasDemandCurve(periodHour);
    electricityDemandCurve = mpc.electricityDemandCuve(periodHour);
    gasDemandPowerProportion = mpc.gasDemandPowerProportion(periodHour);
    offshoreCapacityFactor = offshoreWindAvaliableCapacityCoeff(periodHour,:);
    [solution{iDay}, solution_info{iDay}] = runGEopf_continous(mpc,onshoreWindCapacity,solarCapacity,hydroCapacity, ...
        interconnectorCapacity,offshoreCapacityFactor,electricityDemandCurve,gasDemandCurve,gasDemandPowerProportion,NK,options);
end
nHour = endHour - startHour + 1;
electricityGeneration = zeros(nHour,size(genTypeSet,1));
[onshoreWindCurtailment,offshoreWindCurtailment,gasDemand,GPPgasConsumption] = deal(zeros(nHour,1));
for iDay = 1:nDay
    for iGenType = 1:size(genTypeSet,1)
        typeName = genTypeSet(iGenType);
        genIndex = find(mpc.genType == typeName);
        electricityGeneration((iDay-1)*24+1:iDay*24,iGenType) = sum(solution{iDay}.Pg(:,genIndex'),2) * baseMVA; 
    end
    onshoreWindCurtailment((iDay-1)*24+1:iDay*24) = sum(solution{iDay}.onshoreWindCurtailment,2) * baseMVA; % MW
    offshoreWindCurtailment((iDay-1)*24+1:iDay*24) = sum(solution{iDay}.offshoreWindCurtailment,2) * baseMVA; % MW
    gasDemand((iDay-1)*24+1:iDay*24) = sum(solution{iDay}.energyDemandMultiPeriods,2) * 1e9 / 3600 / 24; % MWh/h
    GPPgasConsumption((iDay-1)*24+1:iDay*24) = sum(solution{iDay}.gasEnergyConsumption,2) * baseMVA; % MW
end
windCurtailment = onshoreWindCurtailment + offshoreWindCurtailment;
interconnectorPower = sum(interconnectorAvaliableCapacity(startHour:endHour,:),2);
electricityDemand = mpc.electricityDemandCuve(startHour:endHour); % MW
end