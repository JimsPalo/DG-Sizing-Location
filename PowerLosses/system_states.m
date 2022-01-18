function [Lij, Si, I, Iij, Sij] = system_states(V, del, Y, fb, tb, BMva)
% Program for Bus Power Injections, Line & Power flows (p.u).

nb = length(V);

[A, B] = pol2cart(del, V);           % Converting polar to rectangular.
Vm = A + 1j*B;

nl = length(fb);                     % No. of Branches..

Iij = zeros(nb,nb);
Sij = zeros(nb,nb);
Si = zeros(nb,1);

%% Bus Current Injections.
 I = Y*Vm;
 
%% Line Current Flows.
for m = 1:nl
    p = fb(m); q = tb(m);
    Iij(p,q) = -(Vm(p) - Vm(q))*Y(p,q); % Y(m,n) = -y(m,n)..
    Iij(q,p) = -Iij(p,q);
end
Iij = sparse(Iij);

%% Line Power Flows.
for m = 1:nb
    for n = 1:nb
        if m ~= n
            Sij(m,n) = Vm(m)*conj(Iij(m,n))*BMva;
        end
    end
end
Sij = sparse(Sij);
 
%% Line Losses.
Lij = zeros(nl,1);
for m = 1:nl
    p = fb(m); q = tb(m);
    Lij(m) = Sij(p,q) + Sij(q,p);
    
end

%% Bus Power Injections.
for i = 1:nb
    for k = 1:nb
        Si(i) = Si(i) + conj(Vm(i))* Vm(k)*Y(i,k)*BMva;
    end
end
Si = conj(Si);

end