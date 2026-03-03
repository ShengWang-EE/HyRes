%% fig 

%%
fig1b = figure;

x = linspace(1,100,100);
y = linspace(prs_i_min,prs_i_max,100);

mesh(x,y,relativePressure');
hold on;
colormap(colors);
xlabel('Hydrogen fraction (%)');
ylabel('Inlet pressure (bar)');
zlabel('Relative outlet pressure');
% view(3);                         % 设置为3D视角
c = colorbar;                        % 显示颜色条
c.Label.String = '';
ax1b = gca;
ax1b.YLim = [prs_i_min,prs_i_max];
x_line1 = linspace(1,100,100);
y_line1 = prs_min * ones(1,100);
z_line1 = relativePressure(:,1);
plot3(x_line1, y_line1, z_line1, 'black-', 'LineWidth', 2);

x_line2 = linspace(1,100,100);
y_line2 = prs_max * ones(1,100);
z_line2 = relativePressure(:,end);
plot3(x_line2, y_line2, z_line2, 'black-', 'LineWidth', 2);
hold off;
print(gcf, 'figs\fig 1b outlet pressure.svg', '-dsvg');

