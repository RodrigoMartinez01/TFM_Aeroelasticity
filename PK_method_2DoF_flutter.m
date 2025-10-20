%% "PK" method for Classical flutter

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
mu = 20; % Mass parameter mu = M/(pi*rho_inf*b^2)
mass = mu*pi*rho_inf*b^2; % Mass of the airfoil (per unit of semi-span)
n_DoF = 2; % Number of degrees of freedom of the system, and thus number of aeroelastic modes considered

%% Definition of non dimensionalize mass and stiffness matrices and natural modes frequencies

M = [1 x_alpha; x_alpha r_alpha^2];
K = [omega_h^2 0; 0 omega_alpha^2*r_alpha^2];
A_nat = K\M;
lambda_nat = eigs(A_nat);
omega_nat = zeros(length(lambda_nat),1);

for i = 1:length(omega_nat)
    omega_nat(i) = 1/sqrt(lambda_nat(i)); % Natural frequencies of the structural modes (with inertia coupling)
end

%% Pk method iterations

N_iter = 2000; 
U_max = 100; % IMPORTANT: Adjust the value depending on the imput parameters.
U_min = 0.01;
U_range = linspace(U_min,U_max,N_iter);
tol = 0.00001;
G = zeros(2*n_DoF, 2*n_DoF);
omega_aero_pole = zeros(N_iter,n_DoF);
damping_aero_pole = zeros(N_iter,n_DoF);
exit = false;

for i = 1:N_iter        

    U_inf = U_range(i);

    for j = 1:n_DoF

        k_0 = b*omega_nat(j)/U_inf; 
        k_iter = k_0;
        epsilon = tol + 0.1;

        while epsilon >= tol  
           
           if k_iter == 0

               disp('The p-k Method fails due to large damping coefficients which moves away the Laplace variable "p" from the imaginary axis at')
               U_inf
               exit = true;   
               break;

           else

           [B0_aero, B1_aero] = Unsteady_aero(k_iter,a);
           G(1:n_DoF, 1:n_DoF) = zeros(n_DoF,n_DoF);
           G(1:n_DoF, (n_DoF+1):2*n_DoF) = eye(n_DoF,n_DoF);
           G((n_DoF+1):2*n_DoF, 1:n_DoF) = -inv(M)*(K - 0.5*rho_inf*U_inf^2/mass * B0_aero);
           G((n_DoF+1):2*n_DoF, (n_DoF+1):2*n_DoF) = 0.5*rho_inf*U_inf*b/mass*inv(M)*B1_aero;
           lambda_aerolastic = eigs(G);
           omega_G = imag(lambda_aerolastic);
           k_G = b/U_inf*omega_G;
           d = zeros(n_DoF,1);

           for n = 1:2*n_DoF
               d(n) = abs(k_G(n) - k_iter);
           end

           [d_min, r] = min(d);
           omega_aero_pole(i,j) = omega_G(r);
           damping_aero_pole(i,j) = real(lambda_aerolastic(r));
           epsilon = abs(k_iter - k_G(r));
           k_iter = k_G(r);

           end
        end
        if exit, break; end
    end
    if exit, break; end
end

% Rationale: For each airspeed U_inf, there is another loop for each
% aeroelastic mode (which are the structural modes modified by the unsteady
% aerodynamics forces). For the first iteration, the reduced frequency is
% assumed to be equal to the corresponding structural one. Then the unsteady
% terms are computed and matrix G and its eigenvalues are calculated. We
% pick the closest eigenvalue to the one used in the previous iteration and
% we calculate the error until it is below the imposed tolerance.

%% Flutter point calculation

idx_mode1 = find(damping_aero_pole(1:end-1,1).*damping_aero_pole(2:end,1) < 0);
idx_mode2 = find(damping_aero_pole(1:end-1,2).*damping_aero_pole(2:end,2) < 0);

if isempty(idx_mode1) && isempty(idx_mode2)

    disp('Atention!: No flutter point in the airspeed range studied, adjust U_max (increase) or change the problem parameters.')
    return

elseif isempty(idx_mode1)

    idx_flutter_mode2 = idx_mode2(1);
    idx_flutter = idx_flutter_mode2;
    U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
        (damping_aero_pole(idx_flutter+1,2) - damping_aero_pole(idx_flutter,2)))*(-damping_aero_pole(idx_flutter,2));
    Freq_flutter = (omega_aero_pole(idx_flutter,2) + ((omega_aero_pole(idx_flutter+1,2) - omega_aero_pole(idx_flutter,2))/...
        (damping_aero_pole(idx_flutter+1,2) - damping_aero_pole(idx_flutter,2)))*(-damping_aero_pole(idx_flutter,2)));

elseif isempty(idx_mode2)

    idx_flutter_mode1 = idx_mode1(1);
    idx_flutter = idx_flutter_mode1;
    U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
        (damping_aero_pole(idx_flutter+1,1) - damping_aero_pole(idx_flutter,1)))*(-damping_aero_pole(idx_flutter,1));
    Freq_flutter = (omega_aero_pole(idx_flutter,1) + ((omega_aero_pole(idx_flutter+1,1) - omega_aero_pole(idx_flutter,1))/...
        (damping_aero_pole(idx_flutter+1,1) - damping_aero_pole(idx_flutter,1)))*(-damping_aero_pole(idx_flutter,1)));
else 

    [idx_flutter, mode] = min([idx_mode1(1) idx_mode2(1)]);
    U_flutter = U_range(idx_flutter) + ((U_range(idx_flutter+1) - U_range(idx_flutter))/...
        (damping_aero_pole(idx_flutter+1,mode) - damping_aero_pole(idx_flutter,mode)))*(-damping_aero_pole(idx_flutter,mode));
    Freq_flutter = (omega_aero_pole(idx_flutter,mode) + ((omega_aero_pole(idx_flutter+1,mode) - omega_aero_pole(idx_flutter,mode))/...
        (damping_aero_pole(idx_flutter+1,mode) - damping_aero_pole(idx_flutter,mode)))*(-damping_aero_pole(idx_flutter,mode)));

end

% Rationale: This section finds the flutter point by linear interpolation
% of the two air speeds U_inf in which the damping of one of the
% aeroelastic modes changes from negative to positive. The interpolation is
% done between these two points and damping = 0, which corresponds to the
% exact flutter point. 

% Dimensionless flutter point

U_F_nondim = U_flutter/(omega_alpha*b);
omega_F_nondim = Freq_flutter/omega_alpha;


%% Vg-Vf plots

figure(1)
    
plot(U_range/(b*omega_alpha), damping_aero_pole(:,1),'Blue','LineWidth',1);
hold on
plot(U_range/(b*omega_alpha), damping_aero_pole(:,2),'Red','LineWidth',1);
hold on 
yline(0, 'k-', 'LineWidth',0.5);
hold on
xline(U_flutter/(b*omega_alpha), 'k--', 'LineWidth',1.2);
hold on
plot(U_flutter/(b*omega_alpha), 0, 'md', 'MarkerSize',7, 'MarkerFaceColor','m')
xlabel('$U_{\infty}/ \omega_\alpha b\mathrm{[-]}$','interpreter','latex')
ylabel('$\sigma_{mode}$','interpreter','latex')
xlim([0 U_inf/(b*omega_alpha)])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

figure(2)

plot(U_range/(b*omega_alpha), omega_aero_pole(:,1)/omega_alpha,'Blue','LineWidth',1);
hold on
plot(U_range/(b*omega_alpha), omega_aero_pole(:,2)/omega_alpha,'Red','LineWidth',1);
hold on
xline(U_flutter/(b*omega_alpha), 'k--', 'LineWidth',1.2);
hold on
yline(Freq_flutter/omega_alpha, 'k--', 'LineWidth',1.2);
hold on 
plot(U_flutter/(b*omega_alpha), Freq_flutter/omega_alpha, 'md', 'MarkerSize',7, 'MarkerFaceColor','m')
xlabel('$U_{\infty}/ \omega_\alpha b\mathrm{[-]}$','interpreter','latex')
ylabel('$\omega/\omega_\alpha\mathrm{[-]}$','interpreter','latex')
xlim([0 U_inf/(b*omega_alpha)])
grid on
grid minor
box on
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 4]);
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman')

%% Unsteady Aerodynamic functions

function [B0_aero, B1_aero] = Unsteady_aero(k,a)
    
    C = Theodorsen_function(k);

    B0_aero = zeros(2,2);
    B1_aero = zeros(2,2);

    B0_aero(1,1) = real(-2*pi*(-k^2 + 2*1i*k*C));
    B0_aero(1,2) = real(-2*pi*(a*k^2 + 1i*k + 2*C*(1 + (0.5-a)*1i*k)));
    B0_aero(2,1) = real(2*pi*(-a*k^2 + 2*1i*k*C*(0.5+a)));
    B0_aero(2,2) = real(2*pi*(-1i*k*(0.5-a) + k^2*(1/8+a^2) + 2*C*(0.5+a)*(1 + (0.5-a)*1i*k)));

    B1_aero(1,1) = 1/k*imag(-2*pi*(-k^2 + 2*1i*k*C));
    B1_aero(1,2) = 1/k*imag(-2*pi*(a*k^2 + 1i*k + 2*C*(1 + (0.5-a)*1i*k)));
    B1_aero(2,1) = 1/k*imag(2*pi*(-a*k^2 + 2*1i*k*C*(0.5+a)));
    B1_aero(2,2) = 1/k*imag(2*pi*(-1i*k*(0.5-a) + k^2*(1/8+a^2) + 2*C*(0.5+a)*(1 + (0.5-a)*1i*k)));

end

function C = Theodorsen_function(k)
    
    C = besselh(1,2,k)/(besselh(1,2,k) + 1i*besselh(0,2,k));

end
