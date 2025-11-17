%% Unsteady Aerodynamic Forces Analysis

close all;
clc;
clear all;

%% Parameters

a = 0; % Non dimensional distance between the rotation axis (Elastic Axis) and the geometrical center (x = 0). a = 0 in Patil's (EA at 50% chord)
alpha = 1; % Amplitude of AoA perturbation. 1 rad in Patil's, 1 deg in Simpson's
k =  0.5; % Reduced frequency. Here is a parameter, in flutter solvers it's part of the solution itself. k = 0.5 for Patil's paper validation
C = Theodorsen_function(k); % Calls Theodorsen's function, and introduces the value of k (reduced freq.)
F = real(C); % This split is done for applying Simpson's mathematical path
G = imag(C);

%% Complex Unsteady Lift Amplitude from Theodorsen's Model

C_l = 2*pi*(a*k^2 + 1i*k + 2*C*(1 + (1/2 - a)*1i*k))*alpha; 
C_l_amplitude = abs(C_l); % Amplitude of the Lift Coefficient 
C_l_phase = angle(C_l); % Phase angle with respect to the AoA motion

%% Complex Steady Lift Amplitude (assuming instantaneous response, and no trailing edge unsteady vortex shedding)

C_l_steady = 4*pi*alpha;
C_l_steady_amplitude  = abs(C_l_steady);
C_l_steady_phase = angle(C_l_steady);

%% Time domain calculations

n = 4; % Number of cycles to plot the response
tau_lim = 2*pi/k*n;
tau = linspace(0, tau_lim, 1000); % Time array 
alpha_t = alpha*real(1i*k*tau); % Real for cos(ks), Imag for sin(ks)

for i = 1:1000
    Cl_t(i) = real(C_l*exp(1i*k*tau(i))); % Change to Imag for alpha * sin(ks) definition as in Simpson's
end

Cl_t_steady = C_l_steady_amplitude*cos(k*tau + C_l_steady_phase); % Change to cos() or to sin() depending on alpha definition

%% Simpsons model (Imaginary as alpha = \bar alpha * sin(ks)); 2 subscript, is for Simpson's paper.

Cl2_t = 2*pi*alpha*(k*cos(k*tau) + a*k^2*sin(k*tau) + 2*F*(sin(k*tau) + (1/2-a)*k*cos(k*tau)) + 2*G*(cos(k*tau) - (1/2-a)*k*sin(k*tau)));

Gamma_1 = 2*(F - k*G*(1/2-a));
Gamma_2 = 2*(G + k*F*(1/2-a)) - k;

for i = 1:1000
    Cd_Lift2_t(i) = Cl2_t(i)*imag(alpha*exp(1i*k*tau(i))); % Imag as alpha = alpha sin(ks) in Simpson's paper
    Cd_LES2_t(i) = -pi*alpha^2*(Gamma_1*sin(k*tau(i)) + Gamma_2*cos(k*tau(i)))^2;
    Cd2_t(i) = Cd_Lift2_t(i) + Cd_LES2_t(i); % Total drag coefficient
end

%% This way is the introduced by my self, using complex exponential expressions instead of working with sin() and cos()

for i = 1:1000
    Cd_Lift_t(i) = real(C_l*exp(1i*k*tau(i)))*real(alpha*exp(1i*k*tau(i)));
    Cd_LES_t(i) = -pi*(real((2*C*(1 + 1i*k*(1/2 - a)) - 1i*k)*exp(1i*k*tau(i))))^2*alpha^2; % (real S)^2 instead of (real(S))^2 --> Garrick's
    Cd_t(i) = Cd_Lift_t(i) + Cd_LES_t(i);
end

% Figures 2, 3 and 4, are just to validate this way of computation of the
% Cd against Simpson's paper formulas. Figure 4, is used to validate the
% unsteady drag in time domain against Patil's paper. Note that for 2
% different alphas sin vs cos (90 deg of phase lag), the phase lag induced
% in unsteady drag is 180 deg, due to Cd acting at twice alpha's frequency.
% To validate against Simpson's, you need to use the imaginary parts, as
% Simpson's formulas are given for alpha = sin(ks). On the other hand, 
% for validating against Patil's, you need to set all the complex to real
% as alpha = cos(ks).

%% Theodorsen's Function Definition

function C = Theodorsen_function(k)
    
    C = besselh(1,2,k)./(besselh(1,2,k) + 1i*besselh(0,2,k));

end

%% Figures

% Figures 1 and 3, can be compared to Patil's paper figures 2c and 2d for
% results validation. Figure 2 plots Cl*alpha and Cs separately as Cd =
% Cl*alpha - Cs; IMPORTANT NOTE: The same as in Patil's, the plots show the
% normalized unsteady lift and drag. Normalized means that they are both
% divided by the steady drag amplitude (2pi*alpha alpha or 4pi*alpha,
% depending if we are adimensionalizing the aerodynamic forces with b = c/2
% ---> 4pi or with c --> 2pi;

figure(1)

plot(tau, Cl_t/C_l_steady_amplitude,'r','LineWidth',1);
hold on
% plot(tau, Cl2_t/C_l_steady_amplitude,'r','LineWidth',1);
% hold on 
plot(tau, Cl_t_steady/C_l_steady_amplitude,'black--','LineWidth',1);
hold on 
legend('$C_{l,unsteady}$', '$C_{l,steady}$','Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex')
ylabel('Aerodynamic Coefficients','interpreter','latex')
xlim([0 tau_lim])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')


% figure(2)
% 
% plot(tau, Cd_LES_t/C_l_steady_amplitude,'b','LineWidth',1);
% hold on
% plot(tau , Cd_LES2_t/C_l_steady_amplitude, 'r', 'LineWidth',1);
% legend('$C_{d,unsteady}$', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
% xlabel('$\tau$','interpreter','latex')
% ylabel('Aerodynamic Coefficients','interpreter','latex')
% xlim([0 tau_lim])
% grid on
% grid minor
% box on
% set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
% set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

figure(2)

plot(tau, Cd_LES_t/C_l_steady_amplitude,'b','LineWidth',1);
hold on
plot(tau, -Cd_Lift_t/C_l_steady_amplitude,'r','LineWidth',1);
% hold on
% plot(tau , Cd_Lift2_t/C_l_steady_amplitude, 'r', 'LineWidth',1);
legend('$C_{s}$', '$\alpha C_{l}$', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex')
ylabel('Aerodynamic Coefficients','interpreter','latex')
xlim([0 tau_lim])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

figure(3)

plot(tau, Cd_t/C_l_steady_amplitude,'black','LineWidth',1);
% hold on
% plot(tau , Cd2_t/C_l_steady_amplitude, 'r', 'LineWidth',1);
legend('$C_{d,unsteady}$', 'Interpreter', 'latex','FontSize', 11, 'Location', 'best');
xlabel('$\tau$','interpreter','latex')
ylabel('Aerodynamic Coefficients','interpreter','latex')
xlim([0 tau_lim])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

