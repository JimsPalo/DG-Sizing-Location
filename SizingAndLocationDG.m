%% DG Location
clear, clc, close all
%% Required input data

% System data
Busdata = xlsread('PSdata.xlsx', 'Busdata');
Linedata = xlsread('PSdata.xlsx', 'Linedata');

%% Data retrieval from Linedata

Nbr = Linedata(:,1);    % Line number
Nl = Linedata(:,2);     % Nl, From bus
Nr = Linedata(:,3);     % Nr, To bus
R = Linedata(:,4);      % R(i), Line resistance
X = Linedata(:,5);      % X(i), Line reactance
Imax = Linedata(:,6);   % Maximum current

%% Data retrieval from Busdata

Busn = Busdata(:,1);    % Bus number
Btype = Busdata(:,2);   % Type of bus 1-Slack, 2-PV, 3-PQ
Pl = Busdata(:,4);      % Pl(i):Load of bus i
Ql = Busdata(:,5);      % Ql(i):Load of bus i


%% Y-Matrix

Ybus = ybus(length(Busn), Nl, Nr, X, R); %Computing Ybus

%% Initial values
arraysize = size(Busn);

Pg = zeros(arraysize);
Qg = zeros(arraysize);

V = ones(arraysize);
del = zeros(arraysize);

%% base Values

Vb = 12.66;     % kV
Sb = 100;       % MVA

Zb = Vb^2/Sb; 

%% Power flow initial case

[V,del] = power_flow(Ybus*Zb, Busn, Btype, V, del, Pg, Qg, Pl/Sb/1e3, Ql/Sb/1e3 );

plot(V)
hold on

[Lij] = system_states(V, del, Ybus*Zb, Nl, Nr, Sb);

Ploss0 = real(sum(Lij));

%% DG boundaries

% ndg = 1;        % Number of DG
% DGType = 2;     % 1 for P, 2 for P and Q

MinP = 0;
MaxP = 3715;    % kW

MinQ = 0;
MaxQ = 2300;    % kVA

%% Handle functions 

% Handle function to genererate several study cases
SizingDg = @(DGType, ndg) sizing_opt(...
                       ndg, DGType, MinP, MaxP, MinQ, MaxQ,...
                       Ybus*Zb, Busn, Btype, V, del, Pg, Qg, Pl/Sb/1e3, Ql/Sb/1e3 , ...
                       Nl, Nr, Sb,...
                       Ploss0);
                   
% Handle function to genererate several reports
SystemReport = @(dgtype, Ndg, xf, Pf, Qf) general_report(...
                    dgtype, Ndg,...
                    xf, Pf,Qf,...
                    Ybus*Zb, Busn, Btype, V, del, Pg, Qg, Pl/Sb/1e3, Ql/Sb/1e3 ,...
                    Nl, Nr, Sb);
%% Type I

% 1 DG with P generation
[xcal, Pgcal, Qgcal] = SizingDg(1, 1);

SystemReport(1, 1, round(xcal), Pgcal/Sb/1e3, Qgcal/Sb/1e3)

% 2 DG with P generation
[xcal, Pgcal, Qgcal] = SizingDg(1, 2);

SystemReport(1, 2, round(xcal), Pgcal/Sb/1e3, Qgcal/Sb/1e3)

% 2 DG with P generation
[xcal, Pgcal, Qgcal] = SizingDg(1, 3);

SystemReport(1, 3, round(xcal), Pgcal/Sb/1e3, Qgcal/Sb/1e3)

legend('0 DG','1 DG', '2 DG', '3 DG')

%% Type II
figure, hold on

% 0 DG with PQ generation
SystemReport(2, 0, 0, 0, 0)

% 1 DG with PQ generation
[xcal, Pgcal, Qgcal] = SizingDg(2, 1);

SystemReport(2, 1, round(xcal), Pgcal/Sb/1e3, Qgcal/Sb/1e3)

% 2 DG with PQ generation
[xcal, Pgcal, Qgcal] = SizingDg(2, 2);

SystemReport(2, 2, round(xcal), Pgcal/Sb/1e3, Qgcal/Sb/1e3)

% 2 DG with PQ generation
[xcal, Pgcal, Qgcal] = SizingDg(2, 3);

SystemReport(2, 3, round(xcal), Pgcal/Sb/1e3, Qgcal/Sb/1e3)

legend('0 DG','1 DG', '2 DG', '3 DG')