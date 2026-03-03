function [prs_i,prs_j,info] = solve_prs_with_constant_linepack(min_alp,K_wey,new_q,newGCV,D,L,Prs_stp,...
    T_stp,T_gas,Z_pipe,prs_max,prs_min)

% Robust equivalent formulation:
% Eq(1) Weymouth is enforced exactly: prs_i = sqrt(prs_j^2 + (new_q/K_wey)^2)
% Eq(2) Solve one scalar equation in prs_j: available_linepack - min_alp = 0
% If exact root is not bracketed, use bounded least-squares fallback.

prs_i = NaN;
prs_j = NaN;
info = struct( ...
    'method', '', ...
    'exitflag', NaN, ...
    'alp_residual', NaN, ...
    'used_fallback', false, ...
    'bracket', [NaN, NaN], ...
    'g_lb', NaN, ...
    'g_ub', NaN);

if K_wey == 0
    info.method = 'invalid_input';
    return
end

r = abs(new_q / K_wey);

pj_lb = max(prs_min, eps);
pj_ub = max(prs_max, pj_lb + 1);
max_ub = max(500, 6 * max(prs_max, pj_lb));
max_expand = 30;

g_lb = g_safe(pj_lb);
g_ub = g_safe(pj_ub);
is_bracketed = sign(g_lb) * sign(g_ub) <= 0;

for k = 1:max_expand
    if is_bracketed || pj_ub >= max_ub
        break
    end
    pj_ub = min(max_ub, 1.25 * pj_ub + 1.0);
    g_ub = g_safe(pj_ub);
    is_bracketed = sign(g_lb) * sign(g_ub) <= 0;
end

info.bracket = [pj_lb, pj_ub];
info.g_lb = g_lb;
info.g_ub = g_ub;

if is_bracketed
    opts_zero = optimset('Display', 'off', 'TolX', 1e-10);
    [pj_sol, fval, exitflag] = fzero(@g_safe, [pj_lb, pj_ub], opts_zero);
    method = 'fzero_bracketed';
else
    opts_min = optimset('Display', 'off', 'TolX', 1e-10);
    [pj_sol, fval_sq, exitflag] = fminbnd(@(p) g_safe(p).^2, pj_lb, pj_ub, opts_min);
    fval = sign(g_safe(pj_sol)) * sqrt(max(fval_sq, 0));
    method = 'fminbnd_fallback';
    info.used_fallback = true;
end

pi_sol = sqrt(pj_sol^2 + r^2);
prs_i = pi_sol;
prs_j = pj_sol;

info.method = method;
info.exitflag = exitflag;
info.alp_residual = fval;

    function val = g_safe(pj)
        alp = available_linepack_from_pj(pj);
        val = alp - min_alp;
        if ~isfinite(val)
            val = 1e9;
        end
    end

    function alp = available_linepack_from_pj(pj)
        pi_ = sqrt(pj^2 + r^2);
        alp = available_linepack(pi_, pj);
    end

    function alp = available_linepack(pi_, pj)
        lp_now = linepack_TJ(pi_, pj);

        pj_min = pj_lb;
        pi_min = sqrt(pj_min^2 + r^2);
        lp_min = linepack_TJ(pi_min, pj_min);

        alp = lp_now - lp_min;
    end

    function lp = linepack_TJ(pi_, pj)
        % average pressure (bar) — compute safely
        denom = (pi_ + pj);
        if denom <= 0
            lp = NaN;
            return
        end
        p_avg = (2/3) * (pi_ + pj - (pi_*pj)/denom);

        V = (pi * D^2 / 4) * L;   % m^3 (geometric)
        p_std_bar = Prs_stp / 1e5;

        V_std = V * (p_avg / p_std_bar) * (T_stp / T_gas) / Z_pipe; % Sm^3
        lp = V_std * newGCV / 1e12;  % TJ
    end
end
