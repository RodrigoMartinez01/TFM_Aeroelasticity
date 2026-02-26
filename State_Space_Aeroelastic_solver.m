%% STATE-SPACE AEROELASTIC SOLVER

clear all;
clc;
close all;

%% Parameters 

rho_inf = 1.225; % Airflow density, depends on the altitude studied
b = 1; % Semi-chord
a = -0.5; % Dimensionless distance from the Elastic Axis to the geometrical half point of the airfoil
r_alpha = sqrt(6/25); % Dimensionless radius of the airfoil
omega_alpha = 5*pi; % Uncoupled natural frequency of torsion mode
omega_h = 3*pi; % Uncoupled natural frequency of bending mode
x_alpha = 0.2; % Dimensionless distance from the CG to the geometrical half point of the airfoil: x_alpha = e - a (in Hodges biblio.)
g_alpha = 0.0; % Damping coefficient of the structure
mu = 20; % Masic parameter mu = M/(pi*rho_inf*b^2)
mass = mu*pi*rho_inf*b^2; % Mass of the airfoil (per unit of semi-span)
U_inf = 33; % Airspeed at which the aeroelastic simulations is carried out

% R.T.Jones Wagner Rational approximation values

psi1 = 0.165; psi2 = 0.335;
eps1 = 0.0445; eps2 = 0.3;
Phi_0 = 1 - psi1 - psi2;

%% Matrices definition

% Structural matrices

M_str = [1 x_alpha; x_alpha r_alpha^2];
k_alpha = b*omega_alpha/U_inf;
K_str = k_alpha^2*[(omega_h/omega_alpha)^2 0; 0 r_alpha^2];

% Aerodynamic matrices

M_aero = 1/mu*[1 -a; -a (1/8 + a^2)];
K_aero = 2*Phi_0/mu*[0 1; 0 -(0.5+a)];
C_aero_NC = 1/mu*[0 1; 0 (0.5-a)];
C_aero_C = 2*Phi_0/mu*[1 (0.5-a); -(0.5+a) -(1/4-a^2)];
C_aero = C_aero_NC + C_aero_C;
E_lag = 2/mu*[-eps1*psi1 -eps2*psi2; eps1*psi1*(0.5+a) eps2*psi2*(0.5+a)];

% Lag states matrices (equations 5 and 6)

D12_lag = [0 1; 0 1];
D34_lag = [1 (0.5-a); 1 (0.5-a)];
D56_lag = [-eps1 0; 0 -eps2];

% Assemble the overall aeroelastic system matrix

M_tot = M_str + M_aero;
K_tot = K_str + K_aero;

A = zeros(6);
A(1:2,3:4) = eye(2,2);
A(3:4,1:2) = -M_tot\K_tot;
A(3:4,3:4) = -M_tot\C_aero;
A(3:4,5:6) = M_tot\E_lag;
A(5:6,1:2) = D12_lag;
A(5:6,3:4) = D34_lag;
A(5:6,5:6) = D56_lag;

%% 3. Time Integration (ODE45)
% Non dimensional time span for the simulation
tau_max = 200; 
tau_span = [0 tau_max];

% Initial conditions (perturbation of plunge or pitch)
alpha_0 = deg2rad(5); 
h_b_0 = 0.2;
X0 = zeros(6,1);
X0(1) = h_b_0;
X0(2) = alpha_0; 

% Aeroelastic global system function definition
sys_model = @(t, x) A * x;

% Resolution of the Sstate-Space ODE system
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[tau, X_sol] = ode45(sys_model, tau_span, X0, options);

% Results (plunge and pitch evolution in time)
h_b_t   = X_sol(:, 1); % Estado x1: h/b
alpha_t = X_sol(:, 2); % Estado x2: alpha

%% Eigenvalue Analysis
[Phi, Lambda] = eig(A);
eigenvalues = diag(Lambda);

% Filtrar frecuencias espurias o muy altas (opcional)
% Nos interesan los modos estructurales (partes imaginarias cercanas a las frecuencias naturales)

% Interpretación rápida
if any(real(eigenvalues) > 0)
    fprintf('\nALERTA: El sistema es INESTABLE (Flutter) a %.2f m/s\n', U_inf);
else
    fprintf('\nEl sistema es ESTABLE a %.2f m/s\n', U_inf);
end

%% -----------------------------------------------------------------------
%  SECTION 2: Flutter velocity and frequency using State-Space Stability 
%  -----------------------------------------------------------------------

% Flutter speed range to be evaluated (m/s)
U_vec = 1:0.5:150; 
stability_log = zeros(length(U_vec), 1); % Guarda la parte real máxima
freq_log = zeros(length(U_vec), 1);      % Guarda la frecuencia asociada

% Matrices constantes (No dependen de U)
% M_str, M_aero, K_aero, C_aero, E_lag son constantes en este modelo adimensionalizado
% La ÚNICA matriz que cambia con U es K_str debido a k_alpha^2

M_tot_loop = M_str + M_aero;
invM = inv(M_tot_loop); % Invertimos una sola vez para ahorrar cómputo (opcional)

found_flutter = false;
U_flutter = NaN;
freq_flutter = NaN;

for i = 1:length(U_vec)
    U_val = U_vec(i);
    
    % Actualizamos Rigidez Estructural (Depende de U)
    % k_alpha = b * omega_alpha / U
    % K_str = (b*wa/U)^2 * [...]
    k_a_loop = b * omega_alpha / U_val;
    K_str_loop = k_a_loop^2 * [(omega_h/omega_alpha)^2 0; 0 r_alpha^2];
    
    % K Total Actualizada
    K_tot_loop = K_str_loop + K_aero;
    
    % Re-ensamblamos la matriz A para esta velocidad
    A_loop = zeros(6);
    A_loop(1:2,3:4) = eye(2,2);
    A_loop(3:4,1:2) = -M_tot\K_tot_loop;
    A_loop(3:4,3:4) = -M_tot\C_aero;
    A_loop(3:4,5:6) =  M_tot\E_lag;
    A_loop(5:6,1:2) = D12_lag;
    A_loop(5:6,3:4) = D34_lag;
    A_loop(5:6,5:6) = D56_lag;
    
    % Autovalores
    lambdas = eig(A_loop);
    
    % Buscamos el autovalor más inestable (mayor parte real)
    [max_real, idx_max] = max(real(lambdas));
    
    stability_log(i) = max_real;
    freq_log(i) = abs(imag(lambdas(idx_max))); % Frecuencia adimensional (k)
    
    % Detección de cruce por cero (Estable -> Inestable)
    if i > 1 && stability_log(i) > 0 && stability_log(i-1) < 0
        % Interpolación lineal para mayor precisión
        slope = (stability_log(i) - stability_log(i-1)) / (U_vec(i) - U_vec(i-1));
        U_flutter = U_vec(i-1) + (0 - stability_log(i-1)) / slope;
        
        % Interpolación de frecuencia
        f_slope = (freq_log(i) - freq_log(i-1)) / (U_vec(i) - U_vec(i-1));
        k_flutter = freq_log(i-1) + f_slope * (U_flutter - U_vec(i-1));
        
        % Convertir k adimensional a Hz físicos: f = k * U / (2*pi*b)
        freq_flutter = k_flutter * U_flutter / (2*pi*b);
        
        found_flutter = true;
        break; % Paramos al encontrar el primer modo inestable
    end
end

fprintf('--------------------------------------------------\n');
if found_flutter
    fprintf('FLUTTER DETECTADO:\n');
    fprintf('Velocidad Crítica (U_f):   %.2f m/s\n', U_flutter);
    fprintf('Frecuencia Crítica (f_f):  %.2f Hz\n', freq_flutter);
else
    fprintf('No se encontró flutter en el rango 0 - %.0f m/s\n', max(U_vec));
end
fprintf('--------------------------------------------------\n');


%% -----------------------------------------------------------------------
%  SECCIÓN 2: ANIMACIÓN VISUAL DEL PERFIL (CORREGIDA Y ROBUSTA)
%  -----------------------------------------------------------------------
% Aseguramos que la figura no se solape con anteriores
if ishandle(100), close(100); end
figure(100); set(gcf, 'Name', 'Aeroelastic Animation', 'Color', 'w');

% 1. RE-CÁLCULO DEL SISTEMA PARA U_INF (Aislamiento del bucle anterior)
% Es vital recalcular K_str y A para U_inf=20, pues el bucle anterior
% dejó la matriz A configurada para U=150 (inestable).

k_alpha_anim = b * omega_alpha / U_inf;
K_str_anim = k_alpha_anim^2 * [(omega_h/omega_alpha)^2 0; 0 r_alpha^2];
K_tot_anim = K_str_anim + K_aero; % K_aero ya estaba definida

% Ensamblaje local de A para la simulación
A_anim = zeros(6);
A_anim(1:2,3:4) = eye(2,2);
A_anim(3:4,1:2) = -M_tot\K_tot_anim;
A_anim(3:4,3:4) = -M_tot\C_aero;
A_anim(3:4,5:6) =  M_tot\E_lag;
A_anim(5:6,1:2) = D12_lag;
A_anim(5:6,3:4) = D34_lag;
A_anim(5:6,5:6) = D56_lag;

% 2. SOLUCIÓN TEMPORAL
% Simulamos 2 segundos reales (ajustar tau_max según U_inf)
% t_real = tau * b / U. Queremos ver unos cuantos ciclos.
tau_max_anim = 200; 
tspan_anim = [0 tau_max_anim];

% Condición Inicial: Pequeño desplazamiento en h y alpha
X0_anim = zeros(6,1);
X0_anim(1) = 0.25;          % h/b = 0.1
X0_anim(2) = deg2rad(5);   % alpha = 5 deg

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[tau_sol, X_sol] = ode45(@(t,x) A_anim*x, tspan_anim, X0_anim, options);

x1 = X_sol(:,1); % h/b
x2 = X_sol(:,2); % alpha (rad)

% 3. GEOMETRÍA DEL PERFIL (NACA 0012 Simplificado)
th = linspace(0, 2*pi, 80);
x_naca = 0.5 * (cos(th) + 1); 
y_naca = 0.12/0.2 * (0.2969*sqrt(x_naca) - 0.1260*x_naca - 0.3516*x_naca.^2 + 0.2843*x_naca.^3 - 0.1015*x_naca.^4);
y_naca(41:end) = -y_naca(41:end); 
% Escalar a dimensiones físicas: de -b a +b
X_foil_local = (x_naca - 0.5) * 2 * b; 
Y_foil_local = y_naca * 2 * b;

% 4. PREPARACIÓN DE LA INTERPOLACIÓN
% Creamos un vector de tiempo constante para que el vídeo sea fluido
dt_frame = 1; % Velocidad de la animación
tau_frames = 0:dt_frame:tau_max_anim;
h_interp   = interp1(tau_sol, x1, tau_frames);
alp_interp = interp1(tau_sol, x2, tau_frames);

% 5. CONFIGURACIÓN GRÁFICA
ax = gca;
axis equal;
grid on;
% Límites fijos para que no salte la cámara
limit = 2.5 * b; 
xlim([-limit, limit]);
ylim([-limit, limit]);
xlabel('X (m)'); ylabel('Y (m) (Positivo arriba)');
hold on;

% Dibujos estáticos
yline(0, 'k--', 'Streamline', 'Color', [0.7 0.7 0.7]);
plot(a*b, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % EA Pivote
text(a*b, -0.3*b, 'EA', 'HorizontalAlignment', 'center', 'FontSize', 8);

% Objeto dinámico (el perfil)
h_patch = patch(X_foil_local, Y_foil_local, 'c', 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.8);
htitle = title('Iniciando...');

fprintf('Reproduciendo animación para U = %.1f m/s...\n', U_inf);

% 6. BUCLE DE ANIMACIÓN
for i = 1:length(tau_frames)
    % Recuperar estado interpolado
    hh = h_interp(i) * b;   % h dimensional (metros)
    aa = alp_interp(i);     % alpha (radianes)

    % --- MATRIZ DE TRANSFORMACIÓN ---
    % 1. Trasladar al origen (EA) -> 2. Rotar -> 3. Trasladar de vuelta -> 4. Plunge

    % Coordenadas relativas al Eje Elástico
    X_rel = X_foil_local - (a*b);
    Y_rel = Y_foil_local;

    % Rotación (Pitch positivo = Nariz Arriba)
    % Nota: En ejes de plot estándar, Y+ es arriba. 
    % Alpha positivo levanta la nariz (x negativo sube, x positivo baja)
    c_a = cos(-aa); 
    s_a = sin(-aa);

    X_rot = X_rel * c_a - Y_rel * s_a;
    Y_rot = X_rel * s_a + Y_rel * c_a;

    % Traslación final: Volver a sumar (a*b) y aplicar h (h>0 es bajar)
    X_final = X_rot + (a*b);
    Y_final = Y_rot - hh; 

    % Actualizar gráficos
    set(h_patch, 'XData', X_final, 'YData', Y_final);
    set(htitle, 'String', sprintf('\\tau: %.1f | U: %.1f m/s | h/b: %.3f | \\alpha: %.1f^o', ...
        tau_frames(i), U_inf, h_interp(i), rad2deg(aa)));

    drawnow; 

    % Control de velocidad (si va muy rápido en PCs potentes)
    pause(0.02);
end
fprintf('Animación completada.\n');

%% 4. Figures 

% --- Figura 1: Respuesta en Pitch (Alpha) ---
figure(1)
plot(tau, rad2deg(alpha_t), 'b', 'LineWidth', 1.5);
hold on
yline(0, 'k--', 'LineWidth', 1); % Referencia cero
legend('$\alpha(\tau)$ (Pitch)', 'Equilibrium', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$ (Dimensionless Time)','interpreter','latex')
ylabel('Pitch Angle $\alpha$ (deg)','interpreter','latex')
title(['Aeroelastic Response at U = ' num2str(U_inf) ' m/s'],'interpreter','latex');
xlim([0 tau_max])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

% --- Figura 2: Respuesta en Plunge (h/b) ---
figure(2)
plot(tau, h_b_t, 'r', 'LineWidth', 1.5);
hold on
yline(0, 'k--', 'LineWidth', 1);
legend('$h/b(\tau)$ (Plunge)', 'Equilibrium', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$ (Dimensionless Time)','interpreter','latex')
ylabel('Plunge Displacement $h/b$','interpreter','latex')
title(['Aeroelastic Response at U = ' num2str(U_inf) ' m/s'],'interpreter','latex');
xlim([0 tau_max])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

% --- Figura 3: Plano de Fase (Alpha vs h/b) - Opcional pero útil ---
figure(3)
plot(h_b_t, rad2deg(alpha_t), 'k', 'LineWidth', 1);
xlabel('$h/b$','interpreter','latex')
ylabel('$\alpha$ (deg)','interpreter','latex')
title('Phase Plane: Coupling','interpreter','latex');
grid on; box on;
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

