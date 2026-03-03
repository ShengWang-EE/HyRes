function [] = preprocess_UK_power_demand_data()
% this script load UK power demand data from grid watch, and reduce the
% resolution of this data to half hour

data1 = readtable('data\power demand (from gridwatch)\power demand 2011-2015.csv');
data2 = readtable('data\power demand (from gridwatch)\power demand 2016-2020.csv');
data3 = readtable('data\power demand (from gridwatch)\power demand 2021-2025.csv');

data0 = [data1;data2;data3];

% 确保 timestamp 是 datetime 格式
data0.timestamp = datetime(data0.timestamp,'InputFormat','yyyy-MM-dd HH:mm:ss');

% 方法1：使用 retime (需要 timetable)
TT = table2timetable(data0,'RowTimes','timestamp');
TT_30min = retime(TT,'regular','mean','TimeStep',minutes(30));

writetimetable(TT_30min,'data\power demand (from gridwatch)\power demand 2011-2025 (halfhour).csv')
end