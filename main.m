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
h = 2;    % [m] Distance of the charge from the center of circle
Q_deg = pi;% [rad] Angle of charge in x-y plane
eps0 = 8.854e-12;
K = 1./(4*pi*eps0);

N = 128; dtheta = 2*pi/N; % Number and size of segments

%% Potentials
% Green function of potential (wire)
GreenFunc = @(theta) ...
    K * 1./((sqrt(a^2 + (R-a).^2 + (R).^2 - 2*(R)*(R-a)*cos(theta))));

% Potential due to the ´Q´ charge (+ image charges)
Q1 = -R/h*Q0;
Q2 = -Q1;
PhiQ  = @(theta) K*Q0./((sqrt(h^2 + (R-a).^2 -2*h*(R-a)*cos(pi - theta))));
PhiQ1 = @(theta) K*Q1./((sqrt((R^4)/(h^2) + (R-a).^2 -2*((R^2)/h)*(R-a)*cos(pi - theta))));
PhiQ2 = K*Q2./(R-a);

%% System matrix
S = zeros(N,N);
for j = 1:N
    t1 = dtheta*(j-1-1/2);
    t2 = dtheta*(j-1+1/2);
    S(1,j) = integral(GreenFunc,t1,t2);
end

% Create the symmetric S matrix
for k = 2:N
    S(k,:) = [fliplr(S(1,2:k)), S(1, 1:(N-k+1))];
end

% Right side of the linear equation system
b = zeros(1,N);
for z = 0:N-1
    b(1,z+1) = PhiQ(z*dtheta) + PhiQ1(z*dtheta) + PhiQ2;

end

% Line charge density
q = S\b'*1e12;

% Integral of the line charge density on the circle
Q = sum(q)*dtheta;


disp(['Az ossztoltes = ', num2str(Q), ' pC']);
hold on
%% Plot results
stairs([0:dtheta:N*dtheta], [q;q(end)], 'Linewidth', 1.5)
%plot(0:dtheta:(N-1)*dtheta,q);
grid on
%axis([ 0,2*pi,0,max(q)*1.05])
xlabel('\theta (rad)', 'Fontname', 'FixedWidth', 'FontSize', 14)
ylabel('q(\theta) (pC/rad)', 'Fontname', 'FixedWidth', 'FontSize', 14)
set(gca, 'Fontname', 'FixedWidth', 'FontSize', 14)