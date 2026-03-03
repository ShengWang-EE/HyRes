% 函数定义：数据空间到归一化图形单位的转换
% 输入参数：
% ax - 坐标轴句柄
% x, y - 数据坐标
% 输出参数：
% norm_x, norm_y - 归一化坐标
function [norm_x, norm_y] = ds2nfu(ax, x, y)
    ax_pos = ax.Position;  % 坐标轴位置
    x_lim = ax.XLim;       % x轴数据范围
    y_lim = ax.YLim;       % y轴数据范围
    
    % 转换
    norm_x = ax_pos(1) + ((x - x_lim(1)) / diff(x_lim)) * ax_pos(3);
    norm_y = ax_pos(2) + ((y - y_lim(1)) / diff(y_lim)) * ax_pos(4);
end
