function [g] = consfcn_gasNodalBalance(Gf,PGs,Qgpp,LCg,Qptg,q_comp,mpc,iGd,i_gpp,i_comp,nGb,nGl,nGs,nGd,n_gpp)
%%
% 1 supply-demand
PGsbus = mpc.Gsou(:,1); 
Cgs_PGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
Qdbus = iGd; 
Cgs_Qd = sparse(Qdbus, (1:nGd)', 1, nGb, nGd); % connection matrix
g = Cgs_PGs*PGs - Cgs_Qd * mpc.Gbus(iGd,3); % supply-demand, Mm3/day
%2 ptg (因为GB里面有重复的，所以得一个一个累加，不能用向量累加）
if ~isempty(mpc.ptg)
    for i = 1:n_ptg
        GB = mpc.ptg(i,1);
        g(GB) = g(GB) + Qptg(i);
    end
end
% 3 gpp
for i = 1:n_gpp
    GB = mpc.gen_extra.gas_bus(i_gpp(i));
    g(GB) = g(GB)-Qgpp(i);
end
% gas flow
for  m = 1:nGl
    fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
    g(fb) = g(fb) - Gf(m);
    g(tb) = g(tb) + Gf(m);
end
% LCg
g(iGd) = g(iGd) + LCg;
% comp
g(mpc.Gline(i_comp,1)) = g(mpc.Gline(i_comp,1)) - q_comp;

end

