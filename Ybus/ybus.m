function Ybus = ybus (nbus, nl, nr, X, R)
% Input: 
%     nbus, buses number
%     nl, from bus 
%     nr, to bus 
%     X, reactance
%     R, resistance
% Output:
%     Ybus, Admittance matrix

%% Impedances
nbr = length(nl); % Branches number
Z = R + 1j*X;
y = ones(nbr,1)./Z; % Branch admittance

%% Initialize Ybus to zero
Ybus = zeros(nbus);

%% Formation of the off diagonal elements
for k = 1:nbr
    Ybus(nl(k),nr(k)) = Ybus(nl(k),nr(k))-y(k);
    Ybus(nr(k),nl(k)) = Ybus(nl(k),nr(k));
end

%% Formation of the diagonal elements
for n = 1:nbus
    for m = (n+1):nbus
        Ybus(n,n) = Ybus(n,n)-Ybus(n,m);
    end
    for m = 1:n-1
        Ybus(n,n) = Ybus(n,n)-Ybus(n,m);
    end
end