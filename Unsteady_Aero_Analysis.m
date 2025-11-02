%% Unsteady Aerodynamic Forces Analysis

close all;
clc;
clear all;

% Parameters

a = -0.5; % Non dimensional distance between the rotation axis (Elastic Axis) and the geometrical center (x = 0) 
alpha = deg2rad(1); % Amplitude of AoA perturbation. 
k =  0.5; % Reduced frequency. Here is a parameter, in flutter solvers it's part of the solution itself
C = Theodorsen_function(k); % Calls Theodorsen function, and introduces the value of k (reduced freq.)

% Complex Unsteady Lift Amplitude from Theodorsen Model

C_l = 2*pi*(a*k^2 + 1i*k + 2*C*(1 + (1/2 - a)*1i*k))*alpha; 
C_l_amplitude = abs(C_l); % Amplitude of the Lift Coefficient 
C_l_phase = angle(C_l); % Phase angle with respect to the AoA motion

% Complex Steady Lift Amplitude (assuming instantaneous response, and no
% trailing edge unsteady vortex shedding)

C_l_steady = 4*pi*alpha;
C_l_steady_amplitude  = abs(C_l_steady);
C_l_steady_phase = angle(C_l_steady);

% Unsteady Drag, from Garrick's Model

C_d = 2*pi*(a*k^2 + 1i*k + 2*C*(1 + (1/2 - a)*1i*k))*alpha^2 - ...
    pi*(2*C*(1 + 1i*k*(1/2 - a)) - 1i*k)^2*alpha^2; % C_l*alpha - C_s
C_d_amplitude = abs(C_d); % Amplitude of the Lift Coefficient 
C_d_phase = angle(C_d); % Phase angle with respect to the AoA motion

% Circulatory and Non Circulatory Lift contributions. Just to study C_l
% behavior in more detail.

C_l_C = 2*pi*(2*C*(1 + (1/2 - a)*1i*k))*alpha; % Circulatory contribution
C_l_C_amplitude = abs(C_l_C); % Amplitude of the Lift Coefficient 
C_l_C_phase = angle(C_l_C); % Phase angle with respect to the AoA motion

C_l_NC = 2*pi*(a*k^2 + 1i*k)*alpha; % Non circulatory contribution
C_l_NC_amplitude = abs(C_l_NC); % Amplitude of the Lift Coefficient 
C_l_NC_phase = angle(C_l_NC); % Phase angle with respect to the AoA motion

% Aerodynamic forces in real Time Domain

n = 4; % Number of cycles to plot the response
tau_lim = 2*pi/k*n;
tau = linspace(0, tau_lim, 1000); % Time array 
alpha_t = alpha*cos(k*tau);
Cl_t = C_l_amplitude*cos(k*tau + C_l_phase);
Cl_t_steady = C_l_steady_amplitude*cos(k*tau + C_l_steady_phase);
Cd_t = C_d_amplitude*cos(2*k*tau + C_d_phase); % The 2 comes from the alpha^2;
Cl_C_t = C_l_C_amplitude*cos(k*tau + C_l_C_phase);
Cl_NC_t = C_l_NC_amplitude*cos(k*tau + C_l_NC_phase);

% Time Plots


figure(1)

plot(tau, alpha_t,'g','LineWidth',1);
hold on
plot(tau, Cl_t,'b','LineWidth',1);
hold on
plot(tau, Cl_t_steady,'black--','LineWidth',1);
hold on 
plot(tau, Cd_t,'r','LineWidth',1);
hold on 
legend('$\alpha(t)$','$C_{l,unsteady}$', '$C_{l,steady}$', '$C_d$','Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex')
ylabel('Aerodynamic Coefficients','interpreter','latex')
xlim([0 tau_lim])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

figure(2)

plot(tau, alpha_t,'g','LineWidth',1);
hold on
plot(tau, Cl_t,'b','LineWidth',1);
hold on
plot(tau, Cl_t_steady,'black--','LineWidth',1);
hold on 
plot(tau, Cl_C_t,'m','LineWidth',1);
hold on 
plot(tau, Cl_NC_t,'r','LineWidth',1);
legend('$\alpha(t)$','$C_{l,unsteady}$', '$C_{l,steady}$', '$C_{l,C}$','$C_{l,NC}$','Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex')
ylabel('Aerodynamic Coefficients','interpreter','latex')
xlim([0 tau_lim])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')


%% Theodorsen's Function Analysis

k_range = linspace(0,5,1000);
C_plot = zeros(1000,1);

for i = 1:1000
    C_plot(i) = Theodorsen_function(k_range((i)));
end

figure(3)
plot(real(C_plot), imag(C_plot),'b','LineWidth',1);

figure(4)
plot(k_range, imag(C_plot),'b','LineWidth',1);

figure(5)
plot(k_range, real(C_plot),'b','LineWidth',1);

function C = Theodorsen_function(k)
    
    C = besselh(1,2,k)./(besselh(1,2,k) + 1i*besselh(0,2,k));

end
