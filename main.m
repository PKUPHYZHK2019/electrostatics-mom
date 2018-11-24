%% Electomagnetic Fields: Electrostatic problem (Integral and MoM)
% "a" sugarú huzalból "R" sugarú kört hajlítunk (a<<R).
% A kör síkjában, annak középpontjától "h" távolságban
% egy "Q" nagyságú ponttöltés áll, levegöben. A huzal össztöltése zérus.
% Írjanak Matlab függvényt, amely a momentum módszerrel kiszámítja a --
% megosztás miatt kialakuló -- vonalmenti töltéssürüséget a huzal
% mentén.
%% Configuration parameters
a = 5e-3;  % [m] Radius of wire
R = 1;     % [m] Radius of the circle 

Q0 = 1e-9; % [C] Charge
h  = 24;    % [m] Distance of the charge from the center of circle
eps0 = 8.854e-12;
K    = 1./(4*pi*eps0);

N = 128; 
dtheta = 2*pi/N; % Number and size of segments

%% Potentials
% Green function of potential (wire)
G = @(theta) 1./1./((sqrt((R-a).^2 + (R).^2 - 2*(R)*(R-a)*cos(theta))));

% Potential due to the 'Q' charge
Q_pot = @(theta) K*Q0./sqrt(h^2 + (R-a).^2 -2 * h *(R-a)*cos(pi - theta));

%% System matrix
A = zeros(N,N);
for j = 1:N
    t1 = dtheta*(j-1-1/2);
    t2 = dtheta*(j-1+1/2);
    A(1,j) = K * integral(G,t1,t2);
end
% Create the symmetric A matrix
for k = 2:N
    A(k,:) = [fliplr(A(1,2:k)), A(1, 1:(N-k+1))];
end

% New S system matrix
S = ones(N+1,N+1);
S(2:N+1,2:N+1) = -1 * A;
S(1, 1) = 0;

% Right side of the linear equation system
b = zeros(N,1);
for z = 1:N
    b(z,1) = Q_pot((z-1)*dtheta);
end

b = [0; b];
% Potential and line charge density 
q = S\b;
% Potential of wire
Phi = q(1);

% Integral of the line charge density on the circle
q = q(2:end)*1e12;
Q = sum(q)*dtheta;

disp(['Integral of charge = ', num2str(Q), ' pC']);
disp(['Potential of wire  = ', num2str(Phi), ' V']);
hold on
%% Plot results
stairs([0:dtheta:N*dtheta], [q;q(end)], 'Linewidth', 1.5)
grid on
axis([ 0,2*pi,min(q)*1.05,max(q)*1.05])
xlabel('\theta (rad)', 'Fontname', 'FixedWidth', 'FontSize', 14)
ylabel('q(\theta) (pC/rad)', 'Fontname', 'FixedWidth', 'FontSize', 14)
set(gca, 'Fontname', 'FixedWidth', 'FontSize', 14)