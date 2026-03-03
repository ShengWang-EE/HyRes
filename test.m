%% visualise the longitude and latitude of places with name, to check if there is any obvious error
clear
clc
GBcor = readtable('GBenergyNetwork v1','sheet','Gbus','range','H2:I79');
GBname = readtable('GBenergyNetwork v1','sheet','Gbus','range','G1:G79');
EBinfo = readtable('buses.csv');
EBnames = EBinfo.name;
EBcor = [EBinfo.y,EBinfo.x];

GBlat = table2array(GBcor(:,1)); GBlon = table2array(GBcor(:,2));
GBnames = string(GBname.NameOfTheLocation);

EBlat = EBinfo.y; EBlon = EBinfo.x;

figure;
geoplot(GBlat, GBlon, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % 画红点
hold on;
geoplot(EBlat, EBlon, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5); % 画蓝点
geobasemap streets; % 加底图
title('地点分布图');

% 标注地名
% for i = 1:length(GBlat)
%     text(GBlat(i), GBlon(i), GBnames{i}, 'FontSize', 10, 'Color', 'r', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% end
% for i = 1:length(EBlat)
%     text(EBlat(i), EBlon(i), EBnames{i}, 'FontSize', 10, 'Color', 'b', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% end

for i = 1:length(GBlat)
    text(GBlat(i), GBlon(i), num2str(i), 'FontSize', 10, 'Color', 'r', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end
for i = 1:length(EBlat)
    text(EBlat(i), EBlon(i), num2str(i), 'FontSize', 10, 'Color', 'b', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end
%% calculate parameters of pipelines
clear
clc

[GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF,rho_stp] = initializeParameters_J15();

f = 0.0025;
z = 1;
r = R.all(1);

L = readtable('GBenergyNetwork v1','sheet','Gline','range','F2:F91').Var1 * 1000;
D = readtable('GBenergyNetwork v1','sheet','Gline','range','G2:G91').Var1;

Theta = sqrt(64 * T_gas * L * Z * f /pi^2 ./ D.^5) * Prs_stp/T_stp / 1e5;
C = 1./sqrt(Theta.^2 / r) * 0.0864; % 转换单位