function [] = plot_UK_gas_network(mpc)

colors = orderedcolors("gem12");
iGd = find(mpc.Gbus(:,3)~=0);
iGs = mpc.Gsou(:,1);
GBlat = mpc.Gbus_extra.Lat; GBlon = mpc.Gbus_extra.Lon;

fig = figure;
geoplot(GBlat, GBlon, '.', 'MarkerSize', 16, 'MarkerEdgeColor', colors(1,:),'LineWidth', 1.5); % 画蓝点
hold on;
geoplot(GBlat(iGd), GBlon(iGd), '.', 'MarkerSize', 16, 'MarkerEdgeColor', colors(2,:),'LineWidth', 1.5);
geoplot(GBlat(iGs), GBlon(iGs), '.', 'MarkerSize', 16, 'MarkerEdgeColor', colors(3,:),'LineWidth', 1.5);
geobasemap streets-light; % 加底图
title('GB gas system in 2020');

% add gas pipe
fb_index = mpc.Gline(:,1); tb_index = mpc.Gline(:,2);
for i = 1:size(mpc.Gline,1)
    geoplot([GBlat(fb_index(i)) GBlat(tb_index(i))], [GBlon(fb_index(i)) GBlon(tb_index(i))], ...
        '-', 'Color',colors(3,:), 'LineWidth', 1.5);
end

% % 标注地名
% for i = 1:length(GBlat)
%     text(GBlat(i), GBlon(i), "GB" + "" + string(i), 'FontSize', 10, 'Color', 'k', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% end

fig.Position = [100,100,600,800];
% exportgraphics(gcf, 'figs/fig UK gas system 2020.pdf', 'ContentType', 'vector');

end