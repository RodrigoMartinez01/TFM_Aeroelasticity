%% "k" method (also call Vg-Vf) for Classical flutter

clear all;
clc;
close all;

%% Parameters of the problem definition

rho_inf = 1.225; % Airflow density, depends on the altitude studied
b = 1; % Semi-chord
a = -0.2; % Dimensionless distance from the Elastic Axis to the geometrical half point of the airfoil
r_alpha = sqrt(6/25); % Dimensionless radius of the airfoil
omega_alpha = 5*pi; % Uncoupled natural frequency of torsion mode
omega_h = 2*pi; % Uncoupled natural frequency of bending mode
x_alpha = 0.1; % Dimensionless distance from the CG to the geometrical half point of the airfoil: x_alpha = e - a (in Hodges biblio.)
g_alpha = 0.0; % Damping coefficient of the structure
mu = 20; % Masic parameter mu = M/(pi*rho_inf*b^2)
mass = mu*pi*rho_inf*b^2; % Mass of the airfoil (per unit of semi-span)
n_DoF = 2; % Number of degrees of freedom of the system, and thus number of aeroelastic modes considered

%% Definition of non dimensionalized mass and stiffness matrices

M = [1 x_alpha; x_alpha r_alpha^2];
K = [(omega_h/omega_alpha)^2 0; 0 r_alpha^2];

%% Implementation of the "k" (also called V-g) method

N = 1000; % Number of iterations. 
k_min = 0.01;
k_max = 10;
k_range = linspace(k_min,k_max,N);
omega_aeroelastic = zeros(N,n_DoF); % Matrix that stores all the aeroelastic frequencies for all modes and "k" iterations
omega_mode = zeros(N,n_DoF);
g_modal = zeros(N,n_DoF);
U_inf = zeros(N,n_DoF);
g_global = zeros(N,n_DoF);

for i = 1:N

    k = k_range(i); % Flutter problem solved from k = 0.01 to k = 1;
    Q_aero = Unsteady_aero(k,a);
    A_matrix = K\(M + 1/(2*pi*mu*k^2)*Q_aero);
    Eval = eigs(A_matrix);

    for j = 1:n_DoF

        omega_mode(i,j) = omega_alpha/sqrt(real(Eval(j)));
        omega_aeroelastic(i,j) = omega_mode(i,j)/(2*pi);
        g_modal(i,j) = imag(Eval(j))/real(Eval(j));
        U_inf(i,j) = omega_mode(i,j)*b/k;
        g_global(i,j) = g_modal(i,j) - g_alpha;
       
    end
end

% Important: In the k method, as k is studied from k_min to k_max on
% interest, as k_min (starting point of the computation of aeroelastic
% modes) is k_min = omega*b/U_inf ---> U_inf is large, and as k increases
% the solution tends to give smaller U_inf (is not exactly a linear
% progression since it is coupled with the omega solution in a non linear 
% aeroelastic behaviour). That's the reason why the aeroelastic modes 
% solutionsn (omegas and dampings) are ordered in inverse fashion, i.e, 
% the first aeroelastic damping corresponds to the largest U_inf obtained
% in the solutions, and the last damping term, to the smallest U_inf.

%% Flutter point calculation

idx_mode1 = find(g_global(1:end-1,1).*g_global(2:end,1) < 0);
idx_mode2 = find(g_global(1:end-1,2).*g_global(2:end,2) < 0);

if isempty(idx_mode1) && isempty(idx_mode2)

   disp('No flutter point in the reduced frequency range studied. Change k_range or problem parameters');

elseif isempty(idx_mode1)

    idx_flutter_mode2 = idx_mode2(length(idx_mode2));
    idx_flutter = idx_flutter_mode2;
    U_flutter = U_inf(idx_flutter,2) + ((U_inf(idx_flutter+1,2) - U_inf(idx_flutter,2))/...
        (g_global(idx_flutter+1,2) - g_global(idx_flutter,2)))*(-g_global(idx_flutter,2));
    Freq_flutter = (omega_mode(idx_flutter,2) + ((omega_mode(idx_flutter+1,2) - omega_mode(idx_flutter,2))/...
        (g_global(idx_flutter+1,2) - g_global(idx_flutter,2)))*(-g_global(idx_flutter,2)));

elseif isempty(idx_mode2)

    idx_flutter_mode1 = idx_mode1(length(idx_mode1));
    idx_flutter = idx_flutter_mode1;
    U_flutter = U_inf(idx_flutter,1) + ((U_inf(idx_flutter+1,1) - U_inf(idx_flutter,1))/...
        (g_global(idx_flutter+1,1) - g_global(idx_flutter,1)))*(-g_global(idx_flutter,1));
    Freq_flutter = (omega_mode(idx_flutter,1) + ((omega_mode(idx_flutter+1,1) - omega_mode(idx_flutter,1))/...
        (g_global(idx_flutter+1,1) - g_global(idx_flutter,1)))*(-g_global(idx_flutter,1)));
else 

    [idx_flutter, mode] = min([idx_mode1(length(idx_mode1)) idx_mode2(length(idx_mode2))]);
    U_flutter = U_inf(idx_flutter, mode) + ((U_inf(idx_flutter+1, mode) - U_inf(idx_flutter, mode))/...
        (g_global(idx_flutter+1,mode) - g_global(idx_flutter,mode)))*(-g_global(idx_flutter,mode));
    Freq_flutter = (omega_mode(idx_flutter,mode) + ((omega_mode(idx_flutter+1,mode) - omega_mode(idx_flutter,mode))/...
        (g_global(idx_flutter+1,mode) - g_global(idx_flutter,mode)))*(-g_global(idx_flutter,mode)));

end

% Dimensionless flutter point

U_F_nondim = U_flutter/(omega_alpha*b);
omega_F_nondim = Freq_flutter/omega_alpha;

%% Vg-Vf plots

U_max = 25;

figure(1)
    
plot(U_inf(:,1)/(b*omega_alpha), g_global(:,1),'Blue','LineWidth',1);
hold on
plot(U_inf(:,2)/(b*omega_alpha), g_global(:,2),'Red', 'LineWidth',1);
hold on 
yline(0, 'k-', 'LineWidth',0.5);
hold on
xline(U_flutter/(b*omega_alpha), 'k--', 'LineWidth',1.2);
hold on
plot(U_flutter/(b*omega_alpha), 0, 'md', 'MarkerSize',7, 'MarkerFaceColor','m')
xlabel('$U_{\infty}/\omega_\alpha b\mathrm{[-]}$','interpreter','latex')
ylabel('$(g-g_{\alpha})\mathrm{[-]}$','interpreter','latex')
xlim([0 U_flutter/(b*omega_alpha)])
%ylim([-0.5 0.2])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

figure(2)

plot(U_inf(:,1)/(b*omega_alpha), omega_mode(:,1)/omega_alpha,'Blue', 'LineWidth',1);
hold on
plot(U_inf(:,2)/(b*omega_alpha), omega_mode(:,2)/omega_alpha,'Red', 'LineWidth',1);
hold on
xline(U_flutter/(b*omega_alpha), 'k--', 'LineWidth',1.2);
hold on
yline(Freq_flutter/omega_alpha, 'k--', 'LineWidth',1.2);
hold on 
plot(U_flutter/(b*omega_alpha), Freq_flutter/omega_alpha, 'md', 'MarkerSize',7, 'MarkerFaceColor','m')
xlabel('$U_{\infty}/\omega_\alpha b\mathrm{[-]}$','interpreter','latex')
ylabel('$\omega/\omega_\alpha \mathrm{[-]}$','interpreter','latex')
xlim([0 U_flutter/(b*omega_alpha)])
%ylim([0.4 2])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

%% Unsteady aerodynamic functions

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
