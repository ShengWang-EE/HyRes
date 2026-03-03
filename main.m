%% 1 single pipeline
clear
clc

[GCV, M, fs, a, R, T_stp, Prs_stp, Z_stp, rho_stp, Z_gas, T_gas, eta, CDF] = initializeParameters();
% 固定出口气压和能量
D = 1; % m
L = 50000; % 50km
f_darcy = 0.012;   % 推荐默认（中性、很常用）
q0 = 20; % mcm/day (i.e., 1e6 Sm3/day)

prs_min = 38; prs_max = 85;% bar
% assuming flow from i to j
prs_j_min = 38; 
prs_j_max = 85; 

n_res = 100;
prs_res = (prs_j_max-prs_j_min) / n_res;

for iHy = 1:100
    hyComp = iHy/100;
    newGCV = [GCV.ng, GCV.hy] * [1-hyComp, hyComp]';
    newR = [R.ng, R.hy] * [1-hyComp, hyComp]';

    for iPrs = 1:n_res
        q(iHy,iPrs) = q0*GCV.ng / newGCV;
        prs_j(iHy,iPrs) = prs_j_min + (iPrs-1) * prs_res; % 38-85 bar    
        [linepack(iHy,iPrs),avaliable_linepack(iHy,iPrs),chargable_linepack(iHy,iPrs),prs_i(iHy,iPrs),K_wey(iHy,iPrs)] = ...
            linepackEnergySingle(D,L,newR,newGCV,Z_gas,Z_stp,Prs_stp,f_darcy,T_stp,T_gas,q(iHy,iPrs),prs_j(iHy,iPrs),prs_min,prs_max);
    end
end
% relative linepack
i_outlimit = (prs_i<prs_min) | (prs_i>prs_max) | (prs_j<prs_min)|(prs_j>prs_max);
% linepack(i_outlimit) = nan;
relativeLinepack = linepack./ repmat(linepack(1,:),[100,1]);
% relative avaliable linepack
% avaliable_linepack(i_outlimit) = nan;
relative_avaliable_linepack = avaliable_linepack./ repmat(avaliable_linepack(1,:),[100,1]);
% relative chargable linepack
% chargable_linepack(i_outlimit) = nan;
relative_chargable_linepack = chargable_linepack./ repmat(chargable_linepack(1,:),[100,1]);

relativePressure = prs_j ./ repmat(prs_j(1,:),[100,1]);

%% 计算不同掺氢比下，若需要恢复到相同的available linepack (气压最低时的水平），气压需要上升到原来的多少倍
%  （在相同energy flow rate的前提下）论证无法通过调节气压使系统恢复
for iHy = 1:100
    for iPrs = 1:n_res
        hyComp = iHy/100;
        min_alp = avaliable_linepack(1,iPrs);
        newGCV = [GCV.ng, GCV.hy] * [1-hyComp, hyComp]';
        newR = [R.ng, R.hy] * [1-hyComp, hyComp]';
        K_SI = (pi/4)*D^2 * (Z_stp*T_stp/Prs_stp) * sqrt( D*newR / (f_darcy*Z_gas*T_gas*L) );  % (Sm3/s)/Pa
        K_wey = K_SI * 1e5 * 0.0864;   % (mcm/day)/bar
        new_q = q0*GCV.ng / newGCV;
    
        [sol_prs_i(iHy,iPrs),sol_prs_j(iHy,iPrs)] = solve_prs_with_constant_linepack(min_alp,K_wey,new_q,newGCV,D,L,Prs_stp,...
        T_stp,T_gas,Z_gas,prs_max,prs_min);
    end
    

end

%% fig 1(b): relative available linepack
figure;
set(gcf, 'Units', 'pixels', 'Position', [120, 120, 620, 620]);
set(gcf, 'Color', 'w');
colors = parula(100);
x = linspace(1,100,100);
hold on;
for i = 1:100
    plot(relative_avaliable_linepack(:,i),'Color',colors(i,:));
end
line([0 20],[relative_avaliable_linepack(20,1),relative_avaliable_linepack(20,1)],'Color', [0.5,0.5,0.5], 'LineStyle', '--');
line([0 20],[relative_avaliable_linepack(20,end),relative_avaliable_linepack(20,end)],'Color', [0.5,0.5,0.5], 'LineStyle', '--');
line([0 100],[relative_avaliable_linepack(100,1),relative_avaliable_linepack(100,1)],'Color', [0.5,0.5,0.5], 'LineStyle', '--');
line([0 100],[relative_avaliable_linepack(100,end),relative_avaliable_linepack(100,end)],'Color', [0.5,0.5,0.5], 'LineStyle', '--');
hold off;
ax1b = gca;
ax1b.YLim = [0.3,1];
ax1b.XLim = [1,100];
main_pos = [0.16, 0.12, 0.62, 0.78];
set(ax1b, 'Position', main_pos);
pbaspect(ax1b, [1 1 1]);
box(ax1b, 'on');
set(ax1b, 'Layer', 'top');
xlabel('Hydrogen fraction (%)');
ylabel('Relative linepack energy');
text(ax1b, -0.12, 1.03, 'b', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Clipping', 'off', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
line(ax1b, [20 20], ax1b.YLim, 'Color', [0.45,0.45,0.45], 'LineStyle', '--', 'LineWidth', 0.9);
line(ax1b, [100 100], ax1b.YLim, 'Color', [0.45,0.45,0.45], 'LineStyle', '--', 'LineWidth', 0.9);

% 放大区域的范围（通过缩小窗口提高放缩比，而不是增大小图尺寸）
xRange_base = [20, 25];
yRange_base = [0.85, 0.9];
zoom_factor = 0.68; % <1 means stronger zoom
x_mid = mean(xRange_base);
y_mid = mean(yRange_base);
x_half = 0.5 * diff(xRange_base) * zoom_factor;
y_half = 0.5 * diff(yRange_base) * zoom_factor;
xRange = [x_mid - x_half, x_mid + x_half];
yRange = [y_mid - y_half, y_mid + y_half];

mainAxes = gca;

hold(mainAxes, 'on');
rectangle(mainAxes, 'Position', [xRange(1), yRange(1), diff(xRange), diff(yRange)], ...
          'EdgeColor', 'black', 'LineWidth', 1, 'LineStyle', '--');
hold(mainAxes, 'off');

inset_pos = [main_pos(1)+0.56*main_pos(3), main_pos(2)+0.58*main_pos(4), ...
             0.34*main_pos(3), 0.272*main_pos(4)]; % keep inset size, only increase zoom ratio
insetAxes = axes('Position', inset_pos); % 小的坐标轴
hold on;
for i = 1:100
    plot(x,relative_avaliable_linepack(:,i),'Color',colors(i,:));
    line([0 20],[0,relative_avaliable_linepack(20,1)],'Color', 'black', 'LineStyle', '--');
end

hold off;
axis(insetAxes, [xRange(1), xRange(2), yRange(1), yRange(2)]);  % 放大到选择的区域
yticks(insetAxes, linspace(yRange(1), yRange(2), 4)); % reduce y-tick density
ytickformat(insetAxes, '%.2f');
box(insetAxes, 'on');
set(insetAxes, 'Layer', 'top');

print(gcf, 'figs\fig 1b relative linepack energy.svg', '-dsvg');
%% fig 1(c)： 画出在一定范围外，不可能通过调节气压使alp恢复
hy_comp = linspace(0,1,100); % x-axis: hydrogen composition
original_outlet_prs = prs_j(1,:)'; % y-axis: original outlet pressure
sol_prs_i(i_outlimit) = nan;

n_points = nnz(~isnan(sol_prs_i));
count = 1;
for i_hy = 1:size(sol_prs_i,1)
    for i_prs = 1:size(sol_prs_i,2)
        if ~isnan(sol_prs_i(i_hy,i_prs))
            plotting_coordinate(count,1) = i_hy/100;
            plotting_coordinate(count,2) = original_outlet_prs(i_prs);
            plotting_coordinate(count,3) = sol_prs_i(i_hy,i_prs);
            count = count + 1;
        end
    end
end

% 画热力图，plotting_coordinate的三个元素分别对应x，y，以及相应位置的颜色
x_plot = plotting_coordinate(:,1); % 氢气掺混比例
y_plot = plotting_coordinate(:,2); % 原始出口压力 (bar)
z_plot = plotting_coordinate(:,3); % 维持相同 ALP 所需入口压力 (bar)

% 使用完整网格：网格内但未被填充的点视为缺失值 (NaN)
x_unique = (1:size(sol_prs_i,1)) / 100;
y_unique = original_outlet_prs(:)';
z_grid = sol_prs_i'; % 行对应y(压力), 列对应x(掺氢比例)

figure;
set(gcf, 'Units', 'pixels', 'Position', [120, 120, 620, 620]); % smaller square figure
imagesc(x_unique * 100, y_unique, z_grid, 'AlphaData', ~isnan(z_grid));
set(gca, 'YDir', 'normal');
set(gca, 'Color', 'w');  % NaN 区域显示为白色
set(gcf, 'Color', 'w');
pbaspect([1 1 1]); % square plotting area
xlim([min(x_unique) * 100, max(x_unique) * 100]);
ylim([min(y_unique), max(y_unique)]);
main_pos_c = [0.16, 0.12, 0.62, 0.78];
set(gca, 'Position', main_pos_c); % keep same main panel size as fig b
% Nature 风格偏好的低饱和顺序色带（浅色，便于后续叠加标注）
c_anchor = [ ...
    1.00 1.00 1.00; ...
    0.94 0.98 0.99; ...
    0.85 0.93 0.95; ...
    0.74 0.87 0.90; ...
    0.63 0.79 0.83];
cmap = interp1(linspace(0,1,size(c_anchor,1)), c_anchor, linspace(0,1,256), 'linear');
colormap(cmap);
% Overlay boundary where required inlet pressure equals prs_max.
z_valid = z_grid(~isnan(z_grid));
x_line = [];
y_line = [];
if ~isempty(z_valid) && min(z_valid) <= prs_max && max(z_valid) >= prs_max
    Cb = contourc(x_unique * 100, y_unique, z_grid, [prs_max prs_max]);
    if ~isempty(Cb)
        idx = 1;
        best_n = 0;
        best_xy = [];
        while idx < size(Cb,2)
            npt = Cb(2,idx);
            xy = Cb(:, idx+1:idx+npt);
            if npt > best_n
                best_n = npt;
                best_xy = xy;
            end
            idx = idx + npt + 1;
        end
        if ~isempty(best_xy)
            xb = best_xy(1,:);
            yb = best_xy(2,:);
            [xb, ord] = sort(xb);
            yb = yb(ord);
            [xb, ia] = unique(xb, 'stable');
            yb = yb(ia);

            hold on;
            if numel(xb) >= 4
                % Smooth boundary for cleaner annotation overlay.
                x_s = linspace(xb(1), xb(end), 1200);
                y_s = interp1(xb, yb, x_s, 'pchip');
                w = max(5, 2 * floor(numel(x_s) / 120) + 1); % odd window
                y_s = smoothdata(y_s, 'movmean', w);
            else
                x_s = xb;
                y_s = yb;
            end

            % Always extend to current axis limits.
            xl = xlim(gca);
            y_lim = ylim(gca);
            if numel(x_s) >= 2
                k = max(2, min(12, numel(x_s)));
                p_l = polyfit(x_s(1:k), y_s(1:k), 1);
                p_r = polyfit(x_s(end-k+1:end), y_s(end-k+1:end), 1);
                y_left = polyval(p_l, xl(1));
                y_right = polyval(p_r, xl(2));
            else
                y_left = y_s(1);
                y_right = y_s(1);
            end
            y_left = min(max(y_left, y_lim(1)), y_lim(2));
            y_right = min(max(y_right, y_lim(1)), y_lim(2));

            x_line = [xl(1), x_s(:).', xl(2)];
            y_line = [y_left, y_s(:).', y_right];
            plot(x_line, y_line, 'Color', [0.15 0.15 0.15], 'LineWidth', 1.6);
            hold off;
        end
    end
end

% Annotate boundary and add hatched infeasible region (upper-right side).
if ~isempty(x_line) && numel(x_line) >= 2
    [x_line, ord_line] = sort(x_line);
    y_line = y_line(ord_line);
    [x_line, ia_line] = unique(x_line, 'stable');
    y_line = y_line(ia_line);

    hold on;
    % Label on boundary: pressure = 85 bar
    i_lab = round(0.72 * numel(x_line));
    i_lab = max(2, min(numel(x_line)-1, i_lab));
    x_lab = x_line(i_lab);
    y_lab = y_line(i_lab);
    y_off = 0.015 * range(ylim(gca));
    text(x_lab, y_lab, sprintf(' Pressure = %.0f bar ', prs_max), ...
        'Color', [0.1 0.1 0.1], 'FontSize', 9, 'Rotation', 0, ...
        'BackgroundColor', 'w', 'Margin', 2, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Position', [x_lab, y_lab + y_off, 0]);

    % Infeasible hatch: whole area above boundary.
    xlim_now = xlim(gca);
    ylim_now = ylim(gca);
    xh = linspace(xlim_now(1), xlim_now(2), 500);
    yb = interp1(x_line, y_line, xh, 'linear', 'extrap');
    yb = min(max(yb, ylim_now(1)), ylim_now(2));

    patch([xh, fliplr(xh)], [yb, ylim_now(2) * ones(size(xh))], ...
        [0.92 0.92 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.12);

    m_h = 0.55;
    b_vals = linspace(ylim_now(1) - m_h * xlim_now(2), ylim_now(2) - m_h * xlim_now(1), 32);
    for ib = 1:numel(b_vals)
        y_h = m_h * xh + b_vals(ib);
        mask = (y_h >= yb) & (y_h <= ylim_now(2));
        if any(mask)
            dmask = diff([false, mask, false]);
            i1 = find(dmask == 1);
            i2 = find(dmask == -1) - 1;
            for iseg = 1:numel(i1)
                plot(xh(i1(iseg):i2(iseg)), y_h(i1(iseg):i2(iseg)), ...
                    'Color', [0.45 0.45 0.45], 'LineWidth', 0.6);
            end
        end
    end

    x_inf = mean(xlim_now);
    yb_mid = interp1(x_line, y_line, x_inf, 'linear', 'extrap');
    y_inf = 0.5 * (yb_mid + ylim_now(2));
    text(x_inf, y_inf, 'Infeasible', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', [0.2 0.2 0.2], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold off;
end
ax = gca;
set(ax, 'Position', main_pos_c); % enforce unchanged main panel size
drawnow;
c = colorbar;
set(ax, 'Position', main_pos_c); % colorbar creation may shrink axes
drawnow;

% Align colorbar with the actual plot-box top/bottom (not outer axes box).
ax.Units = 'normalized';
c.Units = 'normalized';
ax_pos = ax.Position;
pbar = ax.PlotBoxAspectRatio;
pbar_ratio = pbar(1) / pbar(2);
ax_ratio = ax_pos(3) / ax_pos(4);
if ax_ratio > pbar_ratio
    pb_h = ax_pos(4);
    pb_w = pb_h * pbar_ratio;
    pb_x = ax_pos(1) + (ax_pos(3) - pb_w) / 2;
    pb_y = ax_pos(2);
else
    pb_w = ax_pos(3);
    pb_h = pb_w / pbar_ratio;
    pb_x = ax_pos(1);
    pb_y = ax_pos(2) + (ax_pos(4) - pb_h) / 2;
end
cb_gap = 0.035;
cb_w = 0.028;
cbar_pos = [pb_x + pb_w + cb_gap, pb_y, cb_w, pb_h];
set(c, 'Position', cbar_pos);

fs_ax = ax.FontSize;
ylabel(c, 'Required inlet pressure to recover ALP (bar)', 'FontSize', fs_ax);
xlabel('Hydrogen fraction (%)', 'FontSize', fs_ax);
ylabel('Original outlet pressure (bar)', 'FontSize', fs_ax);
text(ax, -0.12, 1.03, 'c', 'Units', 'normalized', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Clipping', 'off', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

print(gcf, 'figs\fig 1c feasible pressure.svg', '-dsvg');
%% fig SI: linepack, ALP, CLP
a=1

