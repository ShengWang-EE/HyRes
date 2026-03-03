%% UK power network
colors = orderedcolors("gem12");

EBinfo = readtable('buses.csv');
EBnames = EBinfo.name;
EBcor = [EBinfo.y,EBinfo.x];

EBlat = EBinfo.y; EBlon = EBinfo.x;

fig = figure;
geoplot(EBlat, EBlon, '.', 'MarkerSize', 16, 'MarkerEdgeColor', colors(3,:),'LineWidth', 1.5); % 画蓝点
hold on;
geobasemap streets-light; % 加底图
title('GB power system in 2020');

% add power line
lineInfo = readtable('PyPSA-GB/data/network/BusesBasedGBsystem/lines.csv');
fbName = lineInfo.bus0; tbName = lineInfo.bus1;
[~,fbIndex] = ismember(lineInfo.bus0,EBinfo.name);
[~,tbIndex] = ismember(lineInfo.bus1,EBinfo.name);
for i = 1:size(lineInfo,1)
    geoplot([EBinfo.y(fbIndex(i)) EBinfo.y(tbIndex(i))], [EBinfo.x(fbIndex(i)) EBinfo.x(tbIndex(i))], ...
        '-', 'Color',colors(3,:), 'LineWidth', 1.5);
end
% add interconnectors
interconnectorInfo = readtable('PyPSA-GB/data/interconnectors/links.csv');
fbName = interconnectorInfo.bus0; tbName = interconnectorInfo.bus1;
[~,fbIndex] = ismember(interconnectorInfo.bus0,EBinfo.name);
[~,tbIndex] = ismember(interconnectorInfo.bus1,EBinfo.name);
for i = 1:size(interconnectorInfo,1)
    geoplot([EBinfo.y(fbIndex(i)) EBinfo.y(tbIndex(i))], [EBinfo.x(fbIndex(i)) EBinfo.x(tbIndex(i))], ...
        '-', 'Color',colors(3,:), 'LineWidth', 1.5);
end

% 标注地名
for i = 1:length(EBlat)
    text(EBlat(i), EBlon(i), "EB" + "" + string(i), 'FontSize', 10, 'Color', 'k', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

fig.Position = [100,100,600,800];
exportgraphics(gcf, 'figs/fig UK power system 2020.pdf', 'ContentType', 'vector');
%%