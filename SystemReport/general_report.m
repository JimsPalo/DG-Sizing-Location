function general_report(...
    DGType, ndg,...
    x, Pdg,Qdg,...
    Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql,...
    Nl, Nr, Sb)
% Generates a graphic and printed report
if x ~= 0
    % Eliminate repeated positions
    [xval, xpos] = unique(x);
    Pdg = Pdg(xpos);
    Qdg = Qdg(xpos);

    for k = 1:length(xval)
        xk = round(xval(k));
        Pg(xk) = Pg(xk) + Pdg(k);

        Qg(xk) = Qg(xk) + Qdg(k);

    end
end

%Power flow
[V, del] = power_flow(Ybus, Busn, Btype, V, del, Pg, Qg, Pl, Ql);

% Losses calculation
[Lij, Si, I, Iij, Sij] = system_states(V, del, Ybus, Nl, Nr, Sb);

% plot voltage
system_report_plot(DGType, ndg, V)

% printed report 
system_report_print(DGType, ndg,...
    x, Pdg*Sb,Qdg*Sb,...
    V, del, Pl*Sb, Ql*Sb, Lij*Sb, Si*Sb, Sij*Sb)



end