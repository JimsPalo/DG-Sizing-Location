function [xcal, Pgcal, Qgcal] = sizing_opt(...
                       ndg, DGType, MinP, MaxP, MinQ, MaxQ,...
                       Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql, ...
                       Nl, Nr, Sb,...
                       Ploss0)

% bus boundaries
lbBus = ones(1, ndg);
ubBus = length(Busn)*ones(1, ndg);

% active power boundaries
lbP = MinP*ones(1, ndg);
ubP = MaxP*ones(1, ndg);

% active power boundaries
lbQ = MinQ*ones(1, ndg);
ubQ = MaxQ*ones(1, ndg);

if DGType ==1
    
    % DG's type 1
    % Objective function
    ObjFunc = @(x) AddDG(...
                       x(1:ndg), x(ndg+1:end)/Sb/1e3, zeros(1,ndg),...
                       Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql, ...
                       Nl, Nr, Sb,...
                       Ploss0);
    
    % Optimization
    x = particleswarm(ObjFunc, ndg*2, [lbBus, lbP], [ubBus, ubP] );
else
    
    % DG's type 2
    % Objective function
    ObjFunc = @(x) AddDG(...
                       x(1:ndg), x(ndg+1:ndg*2)/Sb/1e3, x(ndg*2+1:end)/Sb/1e3,...
                       Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql, ...
                       Nl, Nr, Sb,...
                       Ploss0);
                   
    % Optimization
    x = particleswarm(ObjFunc, ndg*3, [lbBus, lbP, lbQ], [ubBus, ubP, ubQ] );
end

if DGType ==1
    xcal = x(1:ndg);
    Pgcal = x(ndg+1:end);
    Qgcal = zeros(1,ndg);
else
    xcal = x(1:ndg);
    Pgcal = x(ndg+1:ndg*2);
    Qgcal = x(ndg*2+1:end);
end

end