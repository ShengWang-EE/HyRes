function [GCV, M, fs, a, R, T_stp, Prs_stp, Z_stp, rho_stp, Z_gas, T_gas, eta, CDF] = initializeParameters()
% 1 天然气成分
composition_ng = [90 	4 	1 	0.5 	0 	3 	1.5] / 100;
% 2 热值：CH4, C2H6, C3H8, C4H10, H2, N2, CO2
% 电厂效率，电解槽效率常见用 LHV（净热值）口，天然气计量/结算常见用 HHV（GCV, Gross Calorific
% Value）口径（尤其在英国，"GCV"这个词本身就更偏 HHV）。为了统一，都用HHV，但是在电厂计算效率的时候用基于HHV的效率
% 在天然气系统分析里，"标况"没有全球唯一。英国一般是15 °C + 1.01325 bar（101.325 kPa），real dry gas
% https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir82-2401.pdf
GCV.CH4 = 3.77 * 1e7; 
GCV.C2H6 = 6.60 * 1e7;
GCV.C3H8 = 9.40 * 1e7;
GCV.C4H10 = 12.8 * 1e7;
GCV.hy = 12.10 * 1e6;      % J/m3
GCV.N2 = 0;     % J/m3
GCV.CO2 = 0;
GCV.all = [GCV.CH4, GCV.C2H6, GCV.C3H8, GCV.C4H10, GCV.hy, GCV.N2, GCV.CO2];
% GCV.ng_ref = 41.04 * 1e6;     % J/m3
GCV.ng = GCV.all * composition_ng';
% 
M.CH4 = 16 * 1e-3;
M.C2H6 = 30 * 1e-3;
M.C3H8 = 44 * 1e-3;
M.C4H10 = 58 * 1e-3;
M.hy = 2 * 1e-3;           % kg/mol
M.N2 = 28 * 1e-3;
M.CO2 = 44 * 1e-3;
M.all = [M.CH4, M.C2H6, M.C3H8, M.C4H10, M.hy, M.N2, M.CO2];
M.air = 29 * 1e-3;         % kg/mol
% M.ng_ref = 17.478 * 1e-3;     % kg/mol
M.ng = M.all * composition_ng';
%  ﬂame speed 
fs.CH4 = 148;
fs.C2H6 = 301;
fs.C3H8 = 398;
fs.C4H10 = 513;
fs.hy = 339;           % kg/mol
fs.N2 = 0;
fs.CO2 = 0;
fs.All = [fs.CH4, fs.C2H6, fs.C3H8, fs.C4H10, fs.hy, fs.N2, fs.CO2];
fs.ng = fs.All * composition_ng';
% a
a_CH4 = 0.3;
a_C2H6 = 0.75;
a_C3H8 = 0.9759;
a_C4H10 = 1.0928;
a_hy = 0;         
a_N2 = 0.699;
a_CO2 = 0.9759;
a.All = [a_CH4, a_C2H6, a_C3H8, a_C4H10, a_hy, a_N2, a_CO2];
%
Rgas = 8.31446261815324; % J/(mol*K) 这个数对所与气体是不变的，但是有时候R要化成kg的单位，所以与相对分子质量有关
R.air = Rgas / M.air; % J/(kg*K) 约287
R.all = Rgas ./ M.all;
[R.CH4, R.C2H6, R.C3H8, R.C4H10, R.hy, R.N2, R.CO2] = deal(...
    Rgas./M.CH4, Rgas./M.C2H6, Rgas./M.C3H8, Rgas./M.C4H10, Rgas./M.hy, Rgas./M.N2, Rgas./M.CO2);
R.ng = Rgas / M.ng;
T_stp = 288.15;            % K, (15°C)
Prs_stp = 101325;          % Pa
Z_gas = 0.9;               % dimenssionless
Z_stp = 1;
T_gas = 288.15;            % K
eta.electrolysis = 0.6;    % from the energy perspective, the effciency is about 80%
eta.methanation = 0.8;
eta.GFU = 0.5;             % from the energy perspective, 从1/200换算而来
rho_stp = Prs_stp / (Z_stp * R.ng * T_stp); % kg/m3, not accurate value, Prs_stp / (Z * R_ng * T_stp);
%% CDF
CDF.electricity = 1e4; % MW/hour, 大概数值，从jia文章中拿的
CDF.gas = CDF.electricity * GCV.ng / 3600 / 24;
end