function system_report_print(DGType, ndg,...
    x, Pdg,Qdg,...
    V, del, Pl, Ql, Lij, Si, Sij)
%--------------------------------------------------------------------------
% Program for Bus Power Injections, Line & Power flows (p.u) report.
% Input:
%     V, 
%     del, 
%     Pl, 
%     Ql, 
%     Lij, 
%     Si, 
%     Sij
%--------------------------------------------------------------------------

nb = length(V);
nl = length(Lij);               % No. of Branches.

Lpij = real(Lij);               % Active power losses
Lqij = imag(Lij);               % Reactive power losses

% Power injected
Pi = real(Si);
Qi = imag(Si);

% Generated power
Pg = Pi+Pl;
Qg = Qi+Ql;

Del = 180/pi*del;               % Bus Voltage Angles in Degree.

%% Z-matrix 
Z=zeros(nb,10);

k=1;
l=1;
for m=1:nb
    for n=1:nb
        if ( n>m ) && ( Sij(m,n)~=0 )
            Z(k,1)=m;
            Z(k,2)=n;
            Z(k,3)= real(Sij(m,n));
            Z(k,4)= imag(Sij(m,n));
            k=k+1;
        else 
          if Sij(m,n)~=0
            Z(l,5)=m;
            Z(l,6)=n;
            Z(l,7)= real(Sij(m,n));
            Z(l,8)= imag(Sij(m,n));
            l=l+1;
          end
        end
    end
end

for m=1:nl
   Z(m,9)=Lpij(m,1);
   Z(m,10)=Lqij(m,1);
end 

%% Report

disp('#########################################################################################');
disp('-----------------------------------------------------------------------------------------');
disp('                              Newton Raphson Loadflow Analysis ');
disp('-----------------------------------------------------------------------------------------');
disp('| Bus |    V   |  Angle  |     Injection      |     Generation     |          Load      |');
disp('| No  |   pu   |  Degree |    MW   |   MVar   |    MW   |  Mvar    |     MW     |  MVar | ');
for m = 1:nb
    disp('-----------------------------------------------------------------------------------------');
    fprintf('%3g', m); fprintf('  %8.4f', V(m)); fprintf('   %8.4f', Del(m));
    fprintf('  %8.3f', Pi(m)); fprintf('   %8.3f', Qi(m)); 
    fprintf('  %8.3f', Pg(m)); fprintf('   %8.3f', Qg(m)); 
    fprintf('  %8.3f', Pl(m)); fprintf('   %8.3f', Ql(m)); fprintf('\n');
end
disp('-----------------------------------------------------------------------------------------');
fprintf(' Total                  ');
fprintf('  %8.3f', sum(Pi)); 
fprintf('   %8.3f', sum(Qi)); 
fprintf('  %8.3f', sum(Pi+Pl));
fprintf('   %8.3f', sum(Qi+Ql));
fprintf('  %8.3f', sum(Pl)); 
fprintf('   %8.3f', sum(Ql)); 
fprintf('\n');
disp('-----------------------------------------------------------------------------------------');
disp('#########################################################################################');

disp('-------------------------------------------------------------------------------------');
disp('                              Line FLow and Losses ');
disp('-------------------------------------------------------------------------------------');
disp('|From|To |    P    |    Q     | From| To |    P     |   Q     |      Line Loss      |');
disp('|Bus |Bus|   MW    |   MVar   | Bus | Bus|    MW    |  MVar   |     MW   |    MVar  |');
for m = 1:nl
    p = Z(m,1); q = Z(m,2); Pij=Z(m,3); Qij=Z(m,4); Pji= Z(m,7); Qji=Z(m,8); Lpij= Z(m,9); Lqij=Z(m,10);
    disp('-------------------------------------------------------------------------------------');
    fprintf('%4g', p); fprintf('%4g', q); fprintf('  %8.3f', Pij); fprintf('   %8.3f', Qij); 
    fprintf('   %4g', q); fprintf('%4g', p); fprintf('   %8.3f', Pji); fprintf('   %8.3f', Qji);
    fprintf('  %8.3f', Lpij); fprintf('   %8.3f', Lqij);
    fprintf('\n');
end
disp('-------------------------------------------------------------------------------------');
fprintf('   Total Loss                                                 ');
fprintf('  %8.3f', sum(Z(:,9))); fprintf('   %8.3f', sum(Z(:,10)));  fprintf('\n');
disp('-------------------------------------------------------------------------------------');
disp('#####################################################################################');

disp('-------------------------------------------------------------------------------------');
disp('                              DG Size and Locations ');
fprintf('DGs Type: %4g\n', DGType);
fprintf('Number of DG: %4g\n', ndg);
disp('-------------------------------------------------------------------------------------');
disp('|N   |    P    |    Q     | ');
disp('|Bus |   MW    |   MVar   | ');
for m = 1:length(x)
    xp = x(m); Pgp = Pdg(m); Qgp = Qdg(m);
    disp('-------------------------------------------------------------------------------------');
    fprintf('%4g', xp); fprintf('  %8.3f', Pgp); fprintf('   %8.3f', Qgp); 
    fprintf('\n');
end

end