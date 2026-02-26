%% STATE-SPACE AEROELASTIC SOLVER (VARIABLE SURGE - IMPLICIT METHOD)
% Full 7-DOF System with Greenberg & Unsteady Drag
clear all; clc; close all;

%% 1. Parameters 
p.rho_inf = 1.225; 
p.b = 1; 
p.a = -0.5; 
p.r_alpha = sqrt(6/25); 
p.omega_alpha = 5*pi; 
p.omega_h = 3*pi; 
p.x_alpha = 0.2; 
p.mu = 20; 
p.mass = p.mu * pi * p.rho_inf * p.b^2; 

% Surge / Flight parameters
U_inf = 36; % Initial airspeed
p.m_sys = p.mass; % Total aircraft mass (arbitrary scale for surge inertia)
p.Cd0 = 0.015;         % Parasitic drag coefficient

% R.T.Jones Wagner Rational approximation
p.psi1 = 0.165; p.psi2 = 0.335;
p.eps1 = 0.0445; p.eps2 = 0.3;
p.Phi_0 = 1 - p.psi1 - p.psi2;

%% 2. Constant Matrices Definition
% Total Mass Matrix (Structural + Aero Added Mass)
M_str = [1, p.x_alpha; p.x_alpha, p.r_alpha^2];
M_aero = (1/p.mu) * [1, -p.a; -p.a, (1/8 + p.a^2)];
p.M_tot = M_str + M_aero;

% Structural Stiffness (Physical, independent of U)
p.K_s = p.omega_alpha^2 * [(p.omega_h/p.omega_alpha)^2, 0; 0, p.r_alpha^2];

% Aerodynamic Damping Base Matrix
p.D_NC_base = [0, 1; 0, (0.5 - p.a)];

% Lag States Base Matrix
p.Lag_base = [p.eps1*p.psi1, p.eps2*p.psi2; ...
             -p.eps1*p.psi1*(0.5+p.a), -p.eps2*p.psi2*(0.5+p.a)];

%% 3. Initial Conditions & Equilibrium Thrust
% Initial State vector y = [h/b, alpha, dh/dt, dalpha/dt, x5, x6, U]
x1_0 = 0.1;            % Small plunge perturbation
x2_0 = deg2rad(2);      % Small pitch perturbation (2 deg)
x7_0 = U_inf;           % Initial speed

y0 = [x1_0; x2_0; 0; 0; 0; 0; x7_0];

% Calculate initial Drag to set constant Thrust (equilibrium at t=0)
W_0 = x2_0; % Because rates and dot_U are 0 initially
Psi_0 = p.Phi_0 * W_0;
Cl_0 = 2*pi * Psi_0;
Cd_0 = p.Cd0 + x2_0 * Cl_0 - pi*(2*Psi_0)^2;
p.Thrust = 0.5 * p.rho_inf * (2*p.b) * x7_0^2 * Cd_0;

%% 4. Consistent Initial Conditions for ode15i

% ode15i requires f(t, y0, yp0) = 0 exactly. We use MATLAB's decic.
yp0_guess = zeros(7,1); % Guess initial derivatives are zero
fixed_y = [1 1 1 1 1 1 1]; % Fix all initial positions
fixed_yp = [0 0 0 0 0 0 0]; % Let MATLAB calculate all initial derivatives

disp('Calculating consistent initial conditions for implicit solver...');
[y0_mod, yp0_mod] = decic(@(t,y,yp) implicit_aeroelastic_system(t, y, yp, p), ...
                          0, y0, fixed_y, yp0_guess, fixed_yp);

%% 5. Time Integration (ode15i)

t_simulation = 20; % Seconds for the simulation to run
t_span = [0 t_simulation]; % Time for the simulation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

disp('Simulating Aeroelastic System with Surge & Greenberg...');
[t_out, Y_sol] = ode15i(@(t,y,yp) implicit_aeroelastic_system(t, y, yp, p), ...
                        t_span, y0_mod, yp0_mod, options);

%% -----------------------------------------------------------------------
%  SECTION 2: Flutter velocity and frequency using State-Space Stability 
%  (REAL-TIME FORMULATION)
%  -----------------------------------------------------------------------

% Flutter speed range to be evaluated (m/s)
U_vec = 1:0.5:150; 
stability_log = zeros(length(U_vec), 1); % Guarda la parte real máxima (Amortiguamiento)
freq_log = zeros(length(U_vec), 1);      % Guarda la frecuencia asociada en Hz

%% 1. Constant matrices (not dependent on U)
% Inertia Matrices
M_str = [1, p.x_alpha; p.x_alpha, p.r_alpha^2];
M_aero = (1/p.mu) * [1, -p.a; -p.a, (1/8 + p.a^2)];
M_tot = M_str + M_aero;
invM = inv(M_tot); 

% Structural stiffness
K_s = p.omega_alpha^2 * [(p.omega_h/p.omega_alpha)^2, 0; 0, p.r_alpha^2];

% Base aerodynamic matrices
K_circ_base = (2*p.Phi_0/p.mu) * [0, 1; 0, -(0.5+p.a)];
C_circ_base = (2*p.Phi_0/p.mu) * [1, (0.5-p.a); -(0.5+p.a), -(0.5+p.a)*(0.5-p.a)];
C_NC_base   = (1/p.mu) * [0, 1; 0, (0.5-p.a)];
K_lag_base  = (2/p.mu) * [p.eps1*p.psi1, p.eps2*p.psi2; -p.eps1*p.psi1*(0.5+p.a), -p.eps2*p.psi2*(0.5+p.a)];

% Lag States Matrices (Eq 5 y 6)
D12_base = [0, 1; 0, 1];
D34_lag  = [1, (0.5-p.a); 1, (0.5-p.a)]; % Esta NO depende de U
D56_base = [-p.eps1, 0; 0, -p.eps2];

found_flutter = false;
U_flutter = NaN;
freq_flutter = NaN;

%% 2. Loop to search Flutter point (Barrido de U)
for i = 1:length(U_vec)
    U = U_vec(i);
    
    % Matrices that depend on U are actulized
    K_circ = (U/p.b)^2 * K_circ_base;
    C_circ = (U/p.b)   * C_circ_base;
    C_NC   = (U/p.b)   * C_NC_base;
    K_lag  = (U/p.b)^2 * K_lag_base;
    
    % Total matrices
    K_tot = K_s + K_circ;
    C_tot = C_NC + C_circ;
    
    % Lag state matrices are actualized
    D12_lag = (U/p.b) * D12_base;
    D56_lag = (U/p.b) * D56_base;
    
    % 3. State system matrix (6 by 6)
    A_loop = zeros(6);
    A_loop(1:2, 3:4) = eye(2);
    A_loop(3:4, 1:2) = -invM * K_tot;     % F_elast y parte de F_circ
    A_loop(3:4, 3:4) = -invM * C_tot;     % F_amort y parte de F_circ
    A_loop(3:4, 5:6) = -invM * K_lag;     % F_lag
    A_loop(5:6, 1:2) = D12_lag;           % Input posicional a la estela
    A_loop(5:6, 3:4) = D34_lag;           % Input de velocidad a la estela
    A_loop(5:6, 5:6) = D56_lag;           % Dinámica interna de decaimiento
    
    % 4. Eigenvlaues
    lambdas = eig(A_loop);
    
    % Larger real part of the eigenvalues (seraching for the most unstable
    % mode)
    [max_real, idx_max] = max(real(lambdas));
    
    stability_log(i) = max_real;
    
    % Imaginary part for the frequency (Hz)
    freq_log(i) = abs(imag(lambdas(idx_max))) / (2*pi); 
    
    % 5. Cross the zero root (Stable -> Unstable)
    if i > 1 && stability_log(i) >= 0 && stability_log(i-1) < 0 && ~found_flutter
        % Linear interpolation for better accuracy on flutter velocity
        % calculation
        slope = (stability_log(i) - stability_log(i-1)) / (U_vec(i) - U_vec(i-1));
        U_flutter = U_vec(i-1) - stability_log(i-1) / slope;
        
        % Linear interpolation for the flutter frequency
        f_slope = (freq_log(i) - freq_log(i-1)) / (U_vec(i) - U_vec(i-1));
        freq_flutter = freq_log(i-1) + f_slope * (U_flutter - U_vec(i-1));
        
        found_flutter = true;
    end
end

%% 3. Results V-g / V-f
fprintf('--------------------------------------------------\n');
if found_flutter
    fprintf('FLUTTER CLÁSICO DETECTADO (Análisis Lineal):\n');
    fprintf('Velocidad Crítica (U_f):   %.2f m/s\n', U_flutter);
    fprintf('Frecuencia Crítica (f_f):  %.2f Hz\n', freq_flutter);
else
    fprintf('No se encontró flutter en el rango 0 - %.0f m/s\n', max(U_vec));
end
fprintf('--------------------------------------------------\n');

%% 6. Results & Plotting
h_b_t   = Y_sol(:, 1);
alpha_t = Y_sol(:, 2);
U_t     = Y_sol(:, 7);

figure('Name', 'Aeroelastic Response', 'Color', 'w', 'Position', [100 100 800 600]);

subplot(3,1,1)
plot(t_out, rad2deg(alpha_t), 'b', 'LineWidth', 1.5); grid on;
ylabel('Pitch (deg)'); title(['Response at Initial U = ' num2str(U_inf) ' m/s']);

subplot(3,1,2)
plot(t_out, h_b_t, 'r', 'LineWidth', 1.5); grid on;
ylabel('Plunge (h/b)'); 

subplot(3,1,3)
plot(t_out, U_t, 'k', 'LineWidth', 1.5); grid on;
ylabel('Airspeed U (m/s)'); xlabel('Time (s)');



%% =======================================================================
%  IMPLICIT ODE FUNCTION: f(t, y, y') = 0
%  =======================================================================
function res = implicit_aeroelastic_system(~, y, yp, p)
    % 1. Extract States
    x1 = y(1); x2 = y(2); x3 = y(3); x4 = y(4); x5 = y(5); x6 = y(6); x7 = y(7);
    dx1 = yp(1); dx2 = yp(2); dx3 = yp(3); dx4 = yp(4); dx5 = yp(5); dx6 = yp(6); dx7 = yp(7);
    
    if x7 < 0.1, x7 = 0.1; end
    
    % 2. Common Aerodynamic Inputs
    W = x2 + (p.b/x7)*x3 + (0.5 - p.a)*(p.b/x7)*x4 + (p.b*dx7 / x7^2)*x2;
    Psi = p.Phi_0 * W + p.eps1*p.psi1*x5 + p.eps2*p.psi2*x6;

    % 3. Forces for Pitch/Plunge (Eq 3 & 4)
    F_Elast = -p.K_s * [x1; x2];
    F_Amort = -(1/p.mu) * (x7/p.b) * p.D_NC_base * [x3; x4];
    F_Lag   = -(2/p.mu) * (x7/p.b)^2 * p.Lag_base * [x5; x6];
    F_Circ  = -(2*p.Phi_0/p.mu) * (x7/p.b)^2 * W * [1; -(0.5+p.a)];
    F_Greenberg = (1/p.mu) * (dx7/p.b) * x2 * [-1; p.a]; 
    
    LHS_34 = p.M_tot * [dx3; dx4];
    RHS_34 = F_Elast + F_Amort + F_Lag + F_Circ + F_Greenberg;

    % 4. Forces for Surge (Eq 7) 
    % Unsteady Lift 
    Cl_NC = pi * ( (p.b/x7)^2*dx3 + (p.b/x7)*x4 - p.a*(p.b/x7)^2*dx4 + (p.b*dx7/x7^2)*x2 );
    Cl_tot = Cl_NC + 2*pi*Psi;
    
    % Unsteady Drag (Garrick)
    % C_s = 2*pi * (Psi - 0.5 * (b/U) * alpha_dot)^2
    Cs = 2*pi * (Psi - 0.5*(p.b/x7)*x4)^2;
    Cd_tot = p.Cd0 + x2 * Cl_tot - Cs;
    
    D_total = 0.5 * p.rho_inf * (2*p.b) * x7^2 * Cd_tot;

    % 5. Assemble Residual Vector
    res = zeros(7,1);
    res(1) = dx1 - x3; 
    res(2) = dx2 - x4; 
    res(3:4) = LHS_34 - RHS_34; 
    res(5) = dx5 - (x7/p.b)*(W - p.eps1*x5);
    res(6) = dx6 - (x7/p.b)*(W - p.eps2*x6);
    res(7) = p.m_sys * dx7 - p.Thrust + D_total;
end
