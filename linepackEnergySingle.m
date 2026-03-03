function [linepack,avaliable_linepack,chargable_linepack,Prs_i,K_wey] = linepackEnergySingle...
    (D,L,newR,newGCV,Z_pipe,Z_stp,Prs_stp,f_darcy,T_stp,T_gas,q,Prs_j,prs_min,prs_max)

% --- compute coefficient in SI ---
K_SI = (pi/4)*D^2 * (Z_stp*T_stp/Prs_stp) * sqrt( D*newR / (f_darcy*Z_pipe*T_gas*L) );  % (Sm3/s)/Pa
K_wey = K_SI * 1e5 * 0.0864;   % (mcm/day)/bar

Prs_i = sqrt(Prs_j^2 + (q / K_wey)^2 ); % bar

% 
avrgPrs = 2/3 * (Prs_i + Prs_j - Prs_i * Prs_j/(Prs_i + Prs_j));
V = pi * D^2 / 4 * L;
V_std = V * (avrgPrs / (Prs_stp/1e5)) * (T_stp / T_gas) / Z_pipe;
linepack = V_std * newGCV / 1e12; % TJ

% min linepack energy
prs_j_min = prs_min;
prs_i_min = sqrt( prs_j_min^2 + (q / K_wey)^2 );
avrgPrs = 2/3 * (prs_i_min + prs_j_min - prs_i_min * prs_j_min/(prs_i_min + prs_j_min));
V = pi * D^2 / 4 * L;
V_std = V * (avrgPrs / (Prs_stp/1e5)) * (T_stp / T_gas) / Z_pipe;
linepackEnergy_min = V_std * newGCV / 1e12; % TJ
avaliable_linepack = linepack - linepackEnergy_min;
% max linepack energy
prs_i_max = prs_max;
prs_j_max = sqrt( prs_i_max^2 - (q / K_wey)^2 );
avrgPrs = 2/3 * (prs_i_max + prs_j_max - prs_i_max * prs_j_max/(prs_i_max + prs_j_max));
V = pi * D^2 / 4 * L;
V_std = V * (avrgPrs / (Prs_stp/1e5)) * (T_stp / T_gas) / Z_pipe;
linepackEnergy_max = V_std * newGCV / 1e12; % TJ
chargable_linepack = linepackEnergy_max - linepack;

%
end