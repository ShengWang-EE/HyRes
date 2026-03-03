%% 1 load the UK energy system data
clear
clc

mpc_e = data_reader_writer(1); % power system data only
mpc_ge = data_reader_writer(2); % gas system also
save stop1.mat
%% 2 test run models to make sure they work
clear
clc
load stop1.mat
% 2.1 testrun GB power system from mpc provided by other people
GB_power_system = GBreducednetwork();
[test_results] = runopf(GB_power_system);

% 2.2 test run electricity opf with our data
mpopt = mpoption('model','DC');
[test_results] = runopf(mpc_e,mpopt);

% 2.3 test run electricity and gas opf with our data
[test_results_ge] = runopf_ge(mpc_ge);


% [test_results_hge] = runopf_hge(mpc_ge);

%% 3 officially run the models

% 3.1 select a stressful time point in 2025 in winter with high electricity and gas demand



















% % get the load and renewable profiles (so won't have to be infered from
% % weather, and to the node scale
% iYear = 2020;
% historicalYearData = readYearlyHistoricalData(iYear,mpc_e);

% % set hydrogen fraction
% hyFrac = 0:0.1:1;
% nDay = 365;

% for ihy = 1:size(hyFrac,1)
%     for iDay = 1:nDay
%         % get the boundary conditions
%         [dailyProfile] = getDailyProfile(iDay,historicalYearData);
%         % run joint uc
%         [results] = unitCommitmentUK(mpc_e,dailyProfile);
%     end
% end

% % calculate expected system condition

% %