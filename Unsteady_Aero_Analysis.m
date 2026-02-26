%% Unsteady Aerodynamic Forces Analysis: Theodorsen vs Wagner vs State-Space
close all;
clc;
clear all;

%% Parameters
a = 0; % Non dimensional distance between EA and geometric center
alpha = 1; % Amplitude of AoA perturbation
h = 1; % Amplitude of plunge perturbation
k = 0.5; % Reduced frequency
C = Theodorsen_function(k); 
LAG = deg2rad(-90); % Phase lag

%% Pitch and plunge motion of the airfoil
n = 4; % Number of cycles
tau_lim = 2*pi/k*n;
tau = linspace(0, tau_lim, 1000); % Time array 
dtau = tau(2) - tau(1);

% --- Motion Definition (Analytical) ---
alpha_t = real(alpha*exp(1i*k*tau)); 
h_t     = real(h*exp(1i*k*tau - 1i*LAG));

% Velocity 
alpha_dot_t = real(1i*k * alpha * exp(1i*k*tau));
h_dot_t     = real(1i*k * h * exp(1i*k*tau - 1i*LAG));

% Acceleration 
alpha_ddot_t = real(-k^2 * alpha * exp(1i*k*tau));
h_ddot_t     = real(-k^2 * h * exp(1i*k*tau - 1i*LAG));

%% 1. THEODORSEN (Frequency Domain Reference)
C_aero = Unsteady_aero(k,a);
% Time Domain reconstruction from Frequency Domain solution
for i = 1:1000
    Cl_unsteady_t(i) = real(-1*(C_aero(1,1)*h*exp(1i*k*tau(i) - 1i*LAG) + C_aero(1,2)*alpha*exp(1i*k*tau(i)))); 
    Cm_unsteady_t(i) = real(C_aero(2,1)*h*exp(1i*k*tau(i) - 1i*LAG) + C_aero(2,2)*alpha*exp(1i*k*tau(i))); 
end

% Steady Reference
C_l_steady = 4*pi*1;
C_m_steady = 4*pi*1*(0.5+a);
for i = 1:1000
    Cl_steady_t(i) = real(C_l_steady*exp(1i*k*tau(i))); 
    Cm_steady_t(i) = real(C_m_steady*exp(1i*k*tau(i))); 
end

% Drag Reference (Garrick)
for i = 1:1000
    Cd_Lift_t(i) = real(alpha*exp(1i*k*tau(i)))*real(Cl_unsteady_t(i)); 
    Cd_LES_t(i) = -pi*(real(2*C*(alpha*exp(1i*k*tau(i)) + 1i*k*h*exp(1i*k*tau(i) + - 1i*LAG) + 1i*k*(0.5-a)*alpha*exp(1i*k*tau(i)))-1i*k*alpha*exp(1i*k*tau(i))))^2; 
    Cd_t(i) = Cd_Lift_t(i) + Cd_LES_t(i);
end

%% 2. WAGNER (Integral Convolution Method)
% --- Movimiento General (aquí coincide con armónico) ---
sigma = 0; 
% Recalculamos movimiento para asegurar consistencia si cambiaras sigma
alpha_t = alpha * exp(sigma*tau) .* cos(k*tau);
h_t     = h     * exp(sigma*tau) .* cos(k*tau - LAG);

% Derivadas Numéricas (para consistencia general)
alpha_dot_t = gradient(alpha_t, dtau);
h_dot_t     = gradient(h_t, dtau);
alpha_ddot_t = gradient(alpha_dot_t, dtau);
h_ddot_t     = gradient(h_dot_t, dtau);

% Wagner function (R.T. Jones approximation)
phi_wagner = 1 - 0.165*exp(-0.0455*tau) - 0.335*exp(-0.3*tau);

% A) Non-Circulatory Terms (Apparent Mass)
Cl_nc = 2*pi * (h_ddot_t + alpha_dot_t - a*alpha_ddot_t);
Cm_nc = 2*pi * ( a*h_ddot_t - (0.5-a)*alpha_dot_t - (1/8+a^2)*alpha_ddot_t );

% B) Downwash 3/4 chord
w_34 = h_dot_t + alpha_t + (0.5 - a)*alpha_dot_t;
% Derivative for Duhamel
dw_34 = gradient(w_34, dtau); 

% Convolution Loop
Cl_circ_wagner = zeros(1, length(tau));
Cs_wagner = zeros(1, length(tau)); 
for i = 1:length(tau)
    t_curr = tau(i);
    if i == 1
        integral_val = 0;
    else
        taus_int = tau(1:i);
        dw_hist = dw_34(1:i);
        phi_kernel = 1 - 0.165*exp(-0.0455*(t_curr - taus_int)) - 0.335*exp(-0.3*(t_curr - taus_int));
        integral_val = trapz(taus_int, dw_hist .* phi_kernel);
    end
    
    Q_eff = integral_val + w_34(1) * phi_wagner(i); % Incluye transitorio inicial
    
    Cl_circ_wagner(i) = 4*pi * Q_eff;
    
    % Suction
    S_t = 2 * Q_eff - alpha_dot_t(i);
    Cs_wagner(i) = pi * (S_t)^2;
end

Cm_circ_wagner = Cl_circ_wagner * (0.5 + a); 

Cl_wagner_t = Cl_nc + Cl_circ_wagner;
Cm_wagner_t = Cm_nc + Cm_circ_wagner;
Cd_wagner_t = Cl_wagner_t .* alpha_t - Cs_wagner;

%% 3. STATE-SPACE LAG METHOD (Método de Estados)
% Este es el método que usarás en el solver aeroelástico.
% Evita la integral de convolución usando 2 estados internos (x1, x2).

% Constantes de Jones
psi1 = 0.165; eps1 = 0.0455;
psi2 = 0.335; eps2 = 0.3;

% Inicialización de estados de Lag
x1 = 0; 
x2 = 0;

% Arrays para guardar historia
Cl_state_t = zeros(1, length(tau));
Cm_state_t = zeros(1, length(tau));
Cd_state_t = zeros(1, length(tau));

for i = 1:length(tau)
    % 1. Downwash actual (Ya calculado arriba)
    w_curr = w_34(i);
    
    % 2. Ecuación de Salida (Output Equation)
    % Q_eff se reconstruye con los estados actuales y el input actual
    % Fórmula derivada de la transformada de Laplace de la aprox. de Jones
    Q_eff_state = w_curr*(1 - psi1 - psi2) + x1*(psi1*eps1) + x2*(psi2*eps2);
    
    % 3. Fuerzas Circulatorias (Slope 4pi)
    Cl_circ_state = 4*pi * Q_eff_state;
    Cm_circ_state = Cl_circ_state * (0.5 + a);
    
    % 4. Fuerzas Totales (Sumamos la NC que calculamos antes, es idéntica)
    Cl_state_t(i) = Cl_nc(i) + Cl_circ_state;
    Cm_state_t(i) = Cm_nc(i) + Cm_circ_state;
    
    % 5. Drag (Suction)
    S_state = 2 * Q_eff_state - alpha_dot_t(i);
    Cs_state = pi * S_state^2;
    Cd_state_t(i) = Cl_state_t(i) * alpha_t(i) - Cs_state;
    
    % 6. Actualizar Estados (Integración Temporal)
    % Resolvemos dx/dtau = -eps*x + w_34
    % Usamos discretización exacta (asumiendo w constante entre pasos dtau)
    % para máxima estabilidad, aunque Euler funcionaría para dtau pequeño.
    if i < length(tau)
        % Update para el siguiente paso
        x1 = x1*exp(-eps1*dtau) + w_curr * (1 - exp(-eps1*dtau))/eps1; % Más preciso que Euler
        x2 = x2*exp(-eps2*dtau) + w_curr * (1 - exp(-eps2*dtau))/eps2;
    end
end

%% Energetic Analysis
% Cálculo de potencia instantánea
P_Lift_inst   = -Cl_unsteady_t .* h_dot_t;
P_Moment_inst = Cm_unsteady_t .* alpha_dot_t;
P_Total_inst  = P_Lift_inst + P_Moment_inst;

%% Unsteady Aerodynamic functions

function Q_aero = Unsteady_aero(k,a)
    
    C = Theodorsen_function(k);
    Q_aero = zeros(2,2);
    Q_aero(1,1) = -2*pi*(-k^2 + 2*1i*k*C);
    Q_aero(1,2) = -2*pi*(a*k^2 + 1i*k + 2*C*(1 + (0.5-a)*1i*k));
    Q_aero(2,1) = 2*pi*(-a*k^2 + 2*1i*k*C*(0.5+a));
    Q_aero(2,2)= 2*pi*(-1i*k*(0.5-a) + k^2*(1/8+a^2) + 2*C*(0.5+a)*(1 + (0.5-a)*1i*k));

end

function C = Theodorsen_function(k)
    
    C = besselh(1,2,k)/(besselh(1,2,k) + 1i*besselh(0,2,k));

end

%% Figures Comparison
% Figure 1: Lift Coefficient
figure(1)
plot(tau, Cl_unsteady_t/C_l_steady,'r','LineWidth',1); hold on 
plot(tau, Cl_wagner_t/C_l_steady, 'g--', 'LineWidth', 2);
plot(tau, Cl_state_t/C_l_steady, 'b-.', 'LineWidth', 2); % NUEVO: State Space
plot(tau, Cl_steady_t/C_l_steady,'k:','LineWidth',1);
legend('Theodorsen (Freq)', 'Wagner (Integral)', 'State-Space (Lag)', 'Steady', ...
    'Interpreter', 'latex','FontSize', 10, 'Location', 'best');
xlabel('$\tau$','interpreter','latex'); ylabel('$C_l / C_{l,steady}$','interpreter','latex')
title('Lift Coefficient Comparison'); grid on; xlim([0 tau_lim]);
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);

% Figure 2: Moment Coefficient
figure(2)
plot(tau, Cm_unsteady_t/C_m_steady,'r','LineWidth',1); hold on 
plot(tau, Cm_wagner_t/C_m_steady, 'g--', 'LineWidth', 2);
plot(tau, Cm_state_t/C_m_steady, 'b-.', 'LineWidth', 2); % NUEVO: State Space
plot(tau, Cm_steady_t/C_m_steady,'k:','LineWidth',1);
legend('Theodorsen', 'Wagner', 'State-Space', 'Steady', ...
    'Interpreter', 'latex','FontSize', 10, 'Location', 'best');
xlabel('$\tau$','interpreter','latex'); ylabel('$C_m / C_{m,steady}$','interpreter','latex')
title('Moment Coefficient Comparison'); grid on; xlim([0 tau_lim]);
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);

% Figure 3: Drag Components (Reference)
figure(3)
plot(tau, Cd_LES_t/C_l_steady,'b','LineWidth',1); hold on
plot(tau, -Cd_Lift_t/C_l_steady,'r','LineWidth',1);
legend('$C_{s}$', '$\alpha C_{l}$', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex'); ylabel('Drag Coeff Components','interpreter','latex')
title('Unsteady Drag Components (Theodorsen)'); grid on; xlim([0 tau_lim]);
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);

% Figure 4: Total Drag
figure(4)
plot(tau, Cd_t/C_l_steady,'k','LineWidth',1); hold on
plot(tau, Cd_wagner_t/C_l_steady,'g--', 'LineWidth', 2);
plot(tau, Cd_state_t/C_l_steady, 'b-.', 'LineWidth', 2); % NUEVO: State Space
legend('Theodorsen', 'Wagner', 'State-Space', 'Interpreter', 'latex','FontSize', 10, 'Location', 'best');
xlabel('$\tau$','interpreter','latex'); ylabel('$C_d / C_{l,steady}$','interpreter','latex')
title('Total Unsteady Drag Comparison'); grid on; xlim([0 tau_lim]);
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);

% Figure 5: Power
figure(5)
plot(tau, P_Lift_inst, 'b', 'LineWidth', 1); hold on
plot(tau, P_Moment_inst, 'g', 'LineWidth', 1); 
plot(tau, P_Total_inst, 'r', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.5);
legend('$P_{Lift}$', '$P_{Moment}$', '$P_{Total}$', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex'); ylabel('Power ($C_P$)','interpreter','latex');
title('Aerodynamic Power Transfer'); grid on; xlim([0 tau_lim]);
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
