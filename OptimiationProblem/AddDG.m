function Plosscal = AddDG(...
                       x, Pdg, Qdg,...
                       Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql, ...
                       Nl, Nr, Sb,...
                       Ploss0)
%--------------------------------------------------------------------------
% Obtimization function with several inputs.
% Input: 
%     DGType, 1 for P, 2 for P and Q
%     x, node for generation
%     Pdg, Active generation capacity
%     Ybus, Y bus
%     Busn, bus number
%     Btype, bus type
%     V, voltage magnitude in pu
%     del, voltage angle in rad
%     Pg, Active power generated
%     Qg, Reactive power generated
%     Pl, Active load
%     Ql, Reactive Load
%     Nl, from bus array
%     Nr, to bus array 
%     Sb, power base
%     Ploss0, Power loss used as a reference
% Output
%     Plosscal
%--------------------------------------------------------------------------

% Eliminate repeated positions
[xval, xpos] = unique(x);
Pdg = Pdg(xpos);
Qdg = Qdg(xpos);

for k = 1:length(xval)
    xk = round(xval(k));
    Pg(xk) = Pg(xk) + Pdg(k);
    
    Qg(xk) = Qg(xk) + Qdg(k);

end

% Power flow
[V, del] = power_flow(Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql);

% Losses calculation
[Lij] = system_states(V, del, Ybus, Nl, Nr, Sb);

Plossk = real(sum(Lij))*1e6;         % Power loss kW
Ploss0 = Ploss0*1e6;                 % Power loss reference kw   

Plosscal = Plossk/Ploss0*100;        % Taken as reference 100

% Penalization
if not(isempty(V( V > 1.05 | V < 0.95)))
    
    ErrorUp = max(V)-1.04;
    ErrrDown = 0.96-min(V);
    Error = max([ErrorUp, ErrrDown]);
    
    % Penalization calc
    Vpe = (1-1/exp(Error*100))*100;
    
    % Penalization application
    Plosscal = Plosscal + Vpe;
    
end
end