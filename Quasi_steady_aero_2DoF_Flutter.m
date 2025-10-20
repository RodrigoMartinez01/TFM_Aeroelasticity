%% Case 3: Quasi-steady aero, with inertia coupling between torsion and bending, and no structural damping

clear all;
clc;
close all;

%% Parameters definition

rho_inf = 1.225;
mu = 20;
b = 1;
M = mu*pi*rho_inf*b^2;
a = -0.25;
r_alpha = 0.5;
omega_alpha = 8*pi;
omega_h = 4*pi;
x_alpha = 0.1;

%% Aeroelastic system solution in the Laplace domain "p"

N_iter = 1000; 
Damping_11 = zeros(N_iter,1);
Damping_12 = zeros(N_iter,1);
Damping_21 = zeros(N_iter,1);
Damping_22 = zeros(N_iter,1);
Freq_11 = zeros(N_iter,1);
Freq_12 = zeros(N_iter,1);
Freq_21 = zeros(N_iter,1);
Freq_22 = zeros(N_iter,1);
p_11 = zeros(N_iter,1);
p_12 = zeros(N_iter,1);
p_21 = zeros(N_iter,1);
p_22 = zeros(N_iter,1);
Lag_11 = zeros(N_iter,1);
Lag_12 = zeros(N_iter,1);
Lag_21 = zeros(N_iter,1);
Lag_22 = zeros(N_iter,1);
Ratio_11 = zeros(N_iter,1);
Ratio_12 = zeros(N_iter,1);
Ratio_21 = zeros(N_iter,1);
Ratio_22 = zeros(N_iter,1);
U_min = 0.1;
U_max = 100;
U_range = linspace(U_min,U_max,N_iter);

for i = 1:N_iter

    U_inf = U_range(i);
    A4 = r_alpha^2 - x_alpha^2;
    A3 = 2*pi*rho_inf*U_inf*b/M * (r_alpha^2 + x_alpha*(0.5+a));
    A2 = omega_alpha^2*r_alpha^2 + omega_h^2*r_alpha^2 - 2*pi*rho_inf*U_inf^2*(0.5+a)/M - 2*pi*x_alpha*rho_inf*U_inf^2/M;
    A1 = 2*pi*rho_inf*U_inf*b/M*(omega_alpha^2*r_alpha^2 - 2*pi*rho_inf*U_inf^2*(0.5+a)/M) + (2*pi*rho_inf/M)^2*b*U_inf^3*(0.5+a);
    A0 = omega_h^2*omega_alpha^2*r_alpha^2 - 2*pi*omega_h^2*rho_inf*U_inf^2*(0.5+a)/M;
    Characteristic_Polynomial = [A4 A3 A2 A1 A0];
    p = roots(Characteristic_Polynomial);
    p_11(i) = p(1);
    p_12(i) = p(2);
    p_21(i) = p(3);
    p_22(i) = p(4);
    Damping_11(i) = real(p_11(i));
    Damping_12(i) = real(p_12(i));
    Damping_21(i) = real(p_21(i));
    Damping_22(i) = real(p_22(i)); 
    Freq_11(i) = imag(p_11(i));
    Freq_12(i) = imag(p_12(i));
    Freq_21(i) = imag(p_21(i));
    Freq_22(i) = imag(p_22(i));
    % Lag_11(i+1) = angle(-1*((p_11(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_11(i+1))^2*x_alpha));
    % Lag_12(i+1) = angle(-1*((p_12(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_12(i+1))^2*x_alpha));
    % Lag_21(i+1) = angle(-1*((p_21(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_21(i+1))^2*x_alpha));
    % Lag_22(i+1) = angle(-1*((p_22(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_22(i+1))^2*x_alpha));
    % Ratio_11(i+1) = abs(-1*((p_11(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_11(i+1))^2*x_alpha));
    % Ratio_12(i+1) = abs(-1*((p_12(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_12(i+1))^2*x_alpha));
    % Ratio_21(i+1) = abs(-1*((p_21(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_21(i+1))^2*x_alpha));
    % Ratio_22(i+1) = abs(-1*((p_22(i+1))^2*r_alpha^2 + w_alpha^2*r_alpha^2 - pi*rho_inf*v_inf^2/M*(0.5+a))/((p_22(i+1))^2*x_alpha));
    
end

%% Flutter point 

idx_mode1 = find(Damping_11(2:end-1).*Damping_11(3:end) < 0);
idx_mode2 = find(Damping_21(2:end-1).*Damping_21(3:end) < 0);

if isempty(idx_mode1) && isempty(idx_mode2)

   disp('No flutter point in the reduced frequency range studied. Change k_range or problem parameters');

elseif isempty(idx_mode2)

    idx_flutter_mode1 = idx_mode1(length(idx_mode1));
    idx_flutter = idx_flutter_mode1;
    U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
        (Damping_11(idx_flutter+1) - Damping_11(idx_flutter)))*(-Damping_11(idx_flutter));
    Freq_flutter = (Freq_11(idx_flutter) + ((Freq_11(idx_flutter+1) - Freq_11(idx_flutter))/...
        (Damping_11(idx_flutter+1) - Damping_11(idx_flutter)))*(-Damping_11(idx_flutter)))/(2*pi);

elseif isempty(idx_mode1)

    idx_flutter_mode2 = idx_mode2(length(idx_mode2));
    idx_flutter = idx_flutter_mode2;
    U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
        (Damping_21(idx_flutter+1) - Damping_21(idx_flutter)))*(-Damping_21(idx_flutter));
    Freq_flutter = (Freq_21(idx_flutter) + ((Freq_21(idx_flutter+1) - Freq_21(idx_flutter))/...
        (Damping_21(idx_flutter+1) - Damping_21(idx_flutter)))*(-Damping_21(idx_flutter)))/(2*pi);
else

    [idx_flutter, mode] = min([idx_mode1 ;idx_mode2]);
    
    if mode == 1
        U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
            (Damping_11(idx_flutter+1) - Damping_11(idx_flutter)))*(-Damping_11(idx_flutter));
        Freq_flutter = (Freq_11(idx_flutter) + ((Freq_11(idx_flutter+1) - Freq_11(idx_flutter))/...
            (Damping_11(idx_flutter+1) - Damping_11(idx_flutter)))*(-Damping_11(idx_flutter)))/(2*pi);
    else
        U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
            (Damping_21(idx_flutter+1) - Damping_21(idx_flutter)))*(-Damping_21(idx_flutter));
        Freq_flutter = (omega_mode(idx_flutter) + ((omega_mode(idx_flutter+1) - omega_mode(idx_flutter))/...
            (Damping_21(idx_flutter+1) - Damping_21(idx_flutter)))*(-Damping_21(idx_flutter)))/(2*pi);
    end

end



%% Vg and Vf plots

zero_damping = 0:0.1:100;

figure(1)

plot(U_range,Freq_11/(2*pi),'Blue')
hold on
plot(U_range, Freq_21/(2*pi),'Red')
hold on
xline(U_flutter, 'k--', 'LineWidth',1.2);
hold on
yline(Freq_flutter, 'k--', 'LineWidth',1.2);
hold on 
plot(U_flutter, Freq_flutter, 'md', 'MarkerSize',7, 'MarkerFaceColor','m')
xlabel('$U_{\infty}\mathrm{[m/s]}$','interpreter','latex')
ylabel('$f\mathrm{[Hz]}$','interpreter','latex')
xlim([0 U_flutter+2])
%ylim([0.4 2])
grid on
grid minor 
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

figure(2)

plot(U_range, Damping_11,'Blue')
hold on
plot(U_range, Damping_21,'Red')
hold on
yline(0, 'k-', 'LineWidth',0.5);
hold on
xline(U_flutter, 'k--', 'LineWidth',1.2);
hold on
plot(U_flutter, 0, 'md', 'MarkerSize',7, 'MarkerFaceColor','m')
xlabel('$U_{\infty}\mathrm{[m/s]}$','interpreter','latex')
ylabel('$\sigma\mathrm{[-]}$','interpreter','latex')
xlim([0 U_flutter+2])
%ylim([-0.2 0.2])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

% figure(3)
% 
% plot(V,Lag_11)
% hold on
% plot(V, Lag_21)
% hold off;
% 
% figure(4)
% 
% plot(V, Ratio_11)
% hold on
% plot(V, Ratio_21)
% hold off;