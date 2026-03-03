function [Bbus, Bf, Pbusinj, Pfinj] = makeBdc_ws(bus, branch)
% makeBdc - Builds the B matrices and phase shift injections for DC power flow.
% ::
%
%   [BBUS, BF, PBUSINJ, PFINJ] = MAKEBDC(MPC)
%   [BBUS, BF, PBUSINJ, PFINJ] = MAKEBDC(BASEMVA, BUS, BRANCH)
%
%   Returns the B matrices and phase shift injection vectors needed for
%   a DC power flow. The bus real power injections are related to bus
%   voltage angles by
%       P = BBUS * Va + PBUSINJ
%   The real power flows at the from end the lines are related to the bus
%   voltage angles by
%       Pf = BF * Va + PFINJ
%   Does appropriate conversions to p.u.
%   Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
%
%   Example:
%       [Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
%
% See also dcpf.

%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus.BusID ~= (1:nb)')
    error('makeBdc: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

%% for each branch, compute the elements of the branch B matrix and the phase
%% shift "quiescent" injections, where
%%
%%      | Pf |   | Bff  Bft |   | Vaf |   | Pfinj |
%%      |    | = |          | * |     | + |       |
%%      | Pt |   | Btf  Btt |   | Vat |   | Ptinj |
%%
stat = branch.Status;                    %% ones at in-service branches
b = stat ./ (branch.b+1e-8);                    %% series susceptance （防止等于0的情况）
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch.K);                       %% indices of non-zero tap ratios
tap(i) = branch.K(i);                        %% assign non-zero tap ratios
b = b ./ tap;

%% build connection matrix Cft = Cf - Ct for line and from - to buses
f = branch.FromBus;                           %% list of "from" buses
t = branch.ToBus;                           %% list of "to" buses
i = [(1:nl)'; (1:nl)'];                         %% double set of row indices
Cft = sparse(i, [f;t], [ones(nl, 1); -ones(nl, 1)], nl, nb);    %% connection matrix

%% build Bf such that Bf * Va is the vector of real branch powers injected
%% at each branch's "from" bus
Bf = sparse(i, [f; t], [b; -b], nl, nb);    % = spdiags(b, 0, nl, nl) * Cft;

%% build Bbus
Bbus = Cft' * Bf;

%% build phase shift injection vectors
Pfinj = b .* (-branch.Angle * pi/180);      %% injected at the from bus ...
    % Ptinj = -Pfinj;                           %% ... and extracted at the to bus
Pbusinj = Cft' * Pfinj;                         %% Pbusinj = Cf * Pfinj + Ct * Ptinj;
