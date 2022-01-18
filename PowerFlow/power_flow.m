function [V, del, J1,  J2, J3, J4] = power_flow(Y, bus, type, V, del, Pg, Qg, Pl, Ql )

% Program for Newton-Raphson Load Flow Analysis.
% Input:
%     Y, admitance matrix
%         type: array(bus x bus)
% 
%     bus, bus number
%         type: array(1 x bus )
% 
%     type, bus type
%         1, slack
%         2, Pv
%         3, PQ
%         type: array(1 x bus ) 
% 
%     V, initial voltage in pu
%         type: array(1 x bus ) 
% 
%     del, initial voltage angle in radians
%         type: array(1 x bus )
% 
%     Pg, active power generated in pu
%         type: array(1 x bus ) 
% 
%     Qg, reactive power generated in pu
%         type: array(1 x bus ) 
% 
%     Pl, active power suplied in pu
%         type: array(1 x bus ) 
% 
%     Ql, active power suplied in pu 
%         type: array(1 x bus )
% 
% Output:
%         V, Voltaje in pu
%         del, voltaje angle in radians
%         J1, J2, J3, J4,are the jacobian's components
%
%--------------------------------------------------------------------------

nbus = length(bus);         % bus lenght

Psp = Pg - Pl;              % P Specified PGi - PLi
Qsp = Qg - Ql;              % Q Specified QGi - QLi

G = real(Y);                % Conductance matrix..
B = imag(Y);                % Susceptance matrix..

pq = find(type == 3);               % PQ Buses..
npq = length(pq);                   % No. of PQ buses..

Tol = 1;  
Iter = 1;
while (Tol > 1e-10)   % Iteration starting..
    
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    % Calculate P and Q
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end
    
    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];     % Jacobian Matrix..
    
    X = J^-1*M;           % Correction Vector
    dTh = X(1:nbus-1);      % Change in Voltage Angle..
    dV = X(nbus:end);       % Change in Voltage Magnitude..
    
    % Updating State Vectors..
    del(2:nbus) = dTh + del(2:nbus);    % Voltage Angle..
    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);        % Voltage Magnitude..
            k = k+1;
        end
    end
    
    Iter = Iter + 1;
    Tol = max(abs(M));                  % Tolerance..
    
end