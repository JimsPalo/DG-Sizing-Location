function system_report_plot(DGType, ndg, V)
%--------------------------------------------------------------------------
% Program for Bus Power Injections, Line & Power flows (p.u) report.
% Input:
%     V, 
%--------------------------------------------------------------------------
plot(V)

ylabel('Voltaje [pu]')
xlabel('Nodes')

label = string(ndg) +"DG, type " + string(DGType);
legend(label)