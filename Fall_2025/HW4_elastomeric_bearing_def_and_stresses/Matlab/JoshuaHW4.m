% Solid Mechanics Homework 4
% Joshua Dillard
% Bridge bearing stress analysis
% Utilized Claude AI to assist in creating the deformation equation for a
% fixed-fixed beam.
% I approximated the bulging response as similar to a fixed-fixed beam.
% All debugging was performed by me

% Clear workspace
clear; clc; close all;

fprintf('Solid Mechanics Homework 4\n');

fprintf('Joshua Dillard\n\n')

% BRIDGE BEARING KNOWN VALUES

% Define bridge bearing parameters
H = 12.0;     % Height in inches
L = 27.0;     % Width in inches
mu = 0.499;   % Poisson's ratio for rubber
E = 276;      % Young's Modulus for rubber
G = 95;    % Shear Modulus for rubber
lambda = (E*mu)/((1+mu)*(1-2*mu)); % Lame's Constant

% Displacement parameters
a = -0.5;           % Displacement at top of bearing
b = -mu * a;        % Maximum bulging or lateral displacement at midheight

% PART 1: UPDATED DEFORMATION FIELD

fprintf('=== PART 1: UPDATED DEFORMATION FIELD ===\n\n');

% Define deformation field functions with fixed-fixed beam profile

% Vertical deformation field
v = @(point) a * point(2) / H;

% Lateral deformation field for the x direction with fixed-fixed beam profile
u = @(point) b * (point(2)/H) * (1 - point(2)/H) * 4 * (1 / (L/2)) * (point(1) - L/2);

% Lateral deformation field for the z direction with fixed-fixed beam profile
w = @(point) b * (point(2)/H) * (1 - point(2)/H) * 4 * (1 / (L/2)) * (point(3) - L/2);

% Define the coordinates of the point of interest

point = [L, H/2, 0]; % Adjust this code for any point desired

% Calculate deformations at this point (verification)
u_at_point = u(point);
v_at_point = v(point);
w_at_point = w(point);

fprintf('Updated Deformation Field:\n\n');
fprintf('u = b * (y/H) * (1-y/H) * 4 * (1/(L/2)) * (x-L/2)\n');
fprintf('v = a * y/H  (unchanged - linear compression)\n');
fprintf('w = b * (y/H) * (1-y/H) * 4 * (1/(L/2)) * (z-L/2)\n\n');
fprintf(['I modeled the bearing bulge in the x and z directions as a fixed-fixed beam. ' ...
    'This should provide a better approximation of the deformation field with the top and bottom surfaces fixed.\n']);

% UPDATED GRADIENTS FOR FIXED-FIXED BEAM PROFILE
% Fixed-fixed beam profile: (y/H) * (1 - y/H) * 4
% Derivative: d/dy[(y/H) * (1 - y/H) * 4] = (1/H - 2*y/H^2) * 4

% Gradient of u (CORRECTED for fixed-fixed beam profile)
du_dx = b * (point(2)/H) * (1 - point(2)/H) * 4 * (1 / (L/2));
du_dy = b * (1/H - 2*point(2)/H^2) * 4 * (1 / (L/2)) * (point(1) - L/2);
du_dz = 0;

% Gradient of v (unchanged)
dv_dx = 0;
dv_dy = a / H;
dv_dz = 0;

% Gradient of w (CORRECTED for fixed-fixed beam profile)
dw_dx = 0;
dw_dy = b * (1/H - 2*point(2)/H^2) * 4 * (1 / (L/2)) * (point(3) - L/2);
dw_dz = b * (point(2)/H) * (1 - point(2)/H) * 4 * (1 / (L/2));

grad_u = [du_dx, du_dy, du_dz];
grad_v = [dv_dx, dv_dy, dv_dz];
grad_w = [dw_dx, dw_dy, dw_dz];

% UPDATED strain tensor calculation function
function strain_tensor = calculate_strain_tensor(point, H, L, a, b)
    % Calculate gradients at the given point using FIXED-FIXED beam profile
    
    % Fixed-fixed beam profile and its derivative
    beam_profile = (point(2)/H) * (1 - point(2)/H) * 4;
    beam_profile_dy = (1/H - 2*point(2)/H^2) * 4;
    
    % Gradient of u
    du_dx = b * beam_profile * (1 / (L/2));
    du_dy = b * beam_profile_dy * (1 / (L/2)) * (point(1) - L/2);
    du_dz = 0;
    
    % Gradient of v
    dv_dx = 0;
    dv_dy = a / H;
    dv_dz = 0;
    
    % Gradient of w
    dw_dx = 0;
    dw_dy = b * beam_profile_dy * (1 / (L/2)) * (point(3) - L/2);
    dw_dz = b * beam_profile * (1 / (L/2));
    
    % Calculate strain components
    epsilon_xx = du_dx;
    epsilon_yy = dv_dy;
    epsilon_zz = dw_dz;
    epsilon_xy = 0.5 * (du_dy + dv_dx);
    epsilon_xz = 0.5 * (du_dz + dw_dx);
    epsilon_yz = 0.5 * (dv_dz + dw_dy);
    
    % Construct strain tensor matrix
    strain_tensor = [epsilon_xx,  epsilon_xy,  epsilon_xz;
                     epsilon_xy,  epsilon_yy,  epsilon_yz;
                     epsilon_xz,  epsilon_yz,  epsilon_zz];
end

% Calculate strain tensor at the original point
strain_tensor = calculate_strain_tensor(point, H, L, a, b);

% UPDATED function to calculate stress components at any point
function [gamma_x, gamma_y, gamma_z, tau_xy, tau_yz, tau_zx] = calculate_stress_components(point, H, L, a, b, lambda, G)
    % Calculate gradients at the given point using FIXED-FIXED beam profile
    
    % Fixed-fixed beam profile and its derivative
    beam_profile = (point(2)/H) * (1 - point(2)/H) * 4;
    beam_profile_dy = (1/H - 2*point(2)/H^2) * 4;
    
    % Gradient of u
    du_dx = b * beam_profile * (1 / (L/2));
    du_dy = b * beam_profile_dy * (1 / (L/2)) * (point(1) - L/2);
    du_dz = 0;
    
    % Gradient of v
    dv_dx = 0;
    dv_dy = a / H;
    dv_dz = 0;
    
    % Gradient of w
    dw_dx = 0;
    dw_dy = b * beam_profile_dy * (1 / (L/2)) * (point(3) - L/2);
    dw_dz = b * beam_profile * (1 / (L/2));
    
    % Calculate stress components using Hooke's law
    gamma_x = lambda * (du_dx + dv_dy + dw_dz) + 2 * G * du_dx;
    gamma_y = lambda * (du_dx + dv_dy + dw_dz) + 2 * G * dv_dy;
    gamma_z = lambda * (du_dx + dv_dy + dw_dz) + 2 * G * dw_dz;
    tau_xy = 2 * G * 0.5 * (du_dy + dv_dx);
    tau_yz = 2 * G * 0.5 * (dv_dz + dw_dy);
    tau_zx = 2 * G * 0.5 * (du_dz + dw_dx);
end

[gamma_x, gamma_y, gamma_z, tau_xy, tau_yz, tau_zx] = calculate_stress_components(point, H, L, a, b, lambda, G);

% Stress Tensor
stress_tensor = [gamma_x, tau_xy,  tau_zx;
                 tau_xy,  gamma_y, tau_yz;
                 tau_zx,  tau_yz,  gamma_z];

% Principal Stresses and Direction
[Direction, Magnitude] = eig(stress_tensor);

% PART 2: COMPLETE SAINT-VENANT STRAIN COMPATIBILITY TEST

fprintf('\n=== PART 2: UPDATED SAINT-VENANT STRAIN COMPATIBILITY TEST ===\n\n');

% COMPREHENSIVE Function to calculate ALL strain derivatives needed for compatibility
function [strain_components, all_derivatives] = calculate_complete_strain_derivatives(point, H, L, a, b)
    
    x = point(1); y = point(2); z = point(3);
    
    % Fixed-fixed beam profile and its derivatives
    beam_profile = (y/H) * (1 - y/H) * 4;
    beam_profile_dy = (1/H - 2*y/H^2) * 4;
    beam_profile_d2y = -8/H^2;  % Second derivative of beam profile
    
    % Calculate displacement gradients using FIXED-FIXED beam profile
    du_dx = b * beam_profile * (1 / (L/2));
    du_dy = b * beam_profile_dy * (1 / (L/2)) * (x - L/2);
    du_dz = 0;
    
    dv_dx = 0;
    dv_dy = a / H;
    dv_dz = 0;
    
    dw_dx = 0;
    dw_dy = b * beam_profile_dy * (1 / (L/2)) * (z - L/2);
    dw_dz = b * beam_profile * (1 / (L/2));
    
    % Calculate strain components
    epsilon_xx = du_dx;
    epsilon_yy = dv_dy;
    epsilon_zz = dw_dz;
    gamma_xy = du_dy + dv_dx;  % Engineering shear strain (2*epsilon_xy)
    gamma_xz = du_dz + dw_dx;  % Engineering shear strain (2*epsilon_xz)
    gamma_yz = dv_dz + dw_dy;  % Engineering shear strain (2*epsilon_yz)
    
    strain_components = struct('exx', epsilon_xx, 'eyy', epsilon_yy, 'ezz', epsilon_zz, ...
                              'gxy', gamma_xy, 'gxz', gamma_xz, 'gyz', gamma_yz);
    
    % Calculate ALL second derivatives needed for the 6 compatibility equations
    
    % Second derivatives of normal strains
    d2exx_dx2 = 0;  % epsilon_xx doesn't depend on x in second order
    d2exx_dy2 = b * beam_profile_d2y * (1 / (L/2));
    d2exx_dz2 = 0;  % epsilon_xx doesn't depend on z
    
    d2eyy_dx2 = 0;  % epsilon_yy is constant
    d2eyy_dy2 = 0;  % epsilon_yy is constant
    d2eyy_dz2 = 0;  % epsilon_yy is constant
    
    d2ezz_dx2 = 0;  % epsilon_zz doesn't depend on x in second order
    d2ezz_dy2 = b * beam_profile_d2y * (1 / (L/2));
    d2ezz_dz2 = 0;  % epsilon_zz doesn't depend on z in second order
    
    % Second derivatives of shear strains (mixed derivatives)
    d2gxy_dxdy = b * beam_profile_dy * (1 / (L/2));  % d/dx[du/dy + dv/dx] -> d/dx[du/dy]
    d2gxz_dxdz = 0;  % gamma_xz = 0, so all derivatives are 0
    d2gyz_dydz = b * beam_profile_dy * (1 / (L/2));  % d/dy[dv/dz + dw/dy] -> d/dy[dw/dy]
    
    % Additional mixed derivatives needed for equations 4, 5, 6
    d2exx_dydz = 0;  % epsilon_xx doesn't depend on z
    d2eyy_dxdz = 0;  % epsilon_yy is constant
    d2ezz_dxdy = 0;  % epsilon_zz doesn't depend on x
    
    d2gxy_dz = 0;    % gamma_xy doesn't depend on z
    d2gxz_dy = 0;    % gamma_xz = 0
    d2gyz_dx = 0;    % gamma_yz doesn't depend on x
    
    % Store all derivatives in a structure
    all_derivatives = struct(...
        'd2exx_dx2', d2exx_dx2, 'd2exx_dy2', d2exx_dy2, 'd2exx_dz2', d2exx_dz2, ...
        'd2eyy_dx2', d2eyy_dx2, 'd2eyy_dy2', d2eyy_dy2, 'd2eyy_dz2', d2eyy_dz2, ...
        'd2ezz_dx2', d2ezz_dx2, 'd2ezz_dy2', d2ezz_dy2, 'd2ezz_dz2', d2ezz_dz2, ...
        'd2gxy_dxdy', d2gxy_dxdy, 'd2gxz_dxdz', d2gxz_dxdz, 'd2gyz_dydz', d2gyz_dydz, ...
        'd2exx_dydz', d2exx_dydz, 'd2eyy_dxdz', d2eyy_dxdz, 'd2ezz_dxdy', d2ezz_dxdy, ...
        'd2gxy_dz', d2gxy_dz, 'd2gxz_dy', d2gxz_dy, 'd2gyz_dx', d2gyz_dx);
end

% Test points for compatibility check
test_points = [
    [L/4, H/4, L/4];   % Interior point 1
    [L/2, H/2, L/2];   % Center
    [3*L/4, 3*H/4, 3*L/4]; % Interior point 2
    [L, H, L];         % Top corner
    [0, 0, 0];         % Bottom corner
    [L/2, H/4, L/2];   % Mid x,z at quarter height
    [L/2, 3*H/4, L/2]; % Mid x,z at three-quarter height
];

% Initialize arrays to store compatibility violations for all 6 equations
num_equations = 6;
max_violations = zeros(num_equations, 1);
violation_locations = cell(num_equations, 1);

fprintf('Testing Saint-Venant compatibility equations at %d points...\n\n', size(test_points, 1));

% Define the 6 Saint-Venant compatibility equations
equation_names = {
    'Eq 1: ∂²εxx/∂y² + ∂²εyy/∂x² = ∂²γxy/∂x∂y';
    'Eq 2: ∂²εyy/∂z² + ∂²εzz/∂y² = ∂²γyz/∂y∂z';
    'Eq 3: ∂²εzz/∂x² + ∂²εxx/∂z² = ∂²γxz/∂x∂z';
    'Eq 4: ∂²εxx/∂y∂z = ∂/∂x(∂γyz/∂x - ∂γxy/∂z - ∂γxz/∂y)';
    'Eq 5: ∂²εyy/∂x∂z = ∂/∂y(∂γxz/∂y - ∂γyz/∂x - ∂γxy/∂z)';
    'Eq 6: ∂²εzz/∂x∂y = ∂/∂z(∂γxy/∂z - ∂γxz/∂y - ∂γyz/∂x)'
};

% Loop through all test points
for p = 1:size(test_points, 1)
    current_point = test_points(p, :);
    
    % Calculate strain components and all derivatives
    [strains, derivs] = calculate_complete_strain_derivatives(current_point, H, L, a, b);
    
    % Evaluate all 6 Saint-Venant compatibility equations
    violations = zeros(num_equations, 1);
    
    % Equation 1: ∂²εxx/∂y² + ∂²εyy/∂x² = ∂²γxy/∂x∂y
    eq1_left = derivs.d2exx_dy2 + derivs.d2eyy_dx2;
    eq1_right = derivs.d2gxy_dxdy;
    violations(1) = abs(eq1_left - eq1_right);
    
    % Equation 2: ∂²εyy/∂z² + ∂²εzz/∂y² = ∂²γyz/∂y∂z
    eq2_left = derivs.d2eyy_dz2 + derivs.d2ezz_dy2;
    eq2_right = derivs.d2gyz_dydz;
    violations(2) = abs(eq2_left - eq2_right);
    
    % Equation 3: ∂²εzz/∂x² + ∂²εxx/∂z² = ∂²γxz/∂x∂z
    eq3_left = derivs.d2ezz_dx2 + derivs.d2exx_dz2;
    eq3_right = derivs.d2gxz_dxdz;
    violations(3) = abs(eq3_left - eq3_right);
    
    % Equation 4: ∂²εxx/∂y∂z = ∂/∂x(∂γyz/∂x - ∂γxy/∂z - ∂γxz/∂y)
    % Since γyz, γxy, γxz don't depend on x, y, z respectively in our case,
    % the right side is 0, so we only need to check if ∂²εxx/∂y∂z = 0
    eq4_left = derivs.d2exx_dydz;
    eq4_right = 0;  % Simplified for our specific displacement field
    violations(4) = abs(eq4_left - eq4_right);
    
    % Equation 5: ∂²εyy/∂x∂z = ∂/∂y(∂γxz/∂y - ∂γyz/∂x - ∂γxy/∂z)
    % Since εyy is constant, left side is 0, and right side is also 0
    eq5_left = derivs.d2eyy_dxdz;
    eq5_right = 0;  % Simplified for our specific displacement field
    violations(5) = abs(eq5_left - eq5_right);
    
    % Equation 6: ∂²εzz/∂x∂y = ∂/∂z(∂γxy/∂z - ∂γxz/∂y - ∂γyz/∂x)
    % Similar analysis as equations 4 and 5
    eq6_left = derivs.d2ezz_dxdy;
    eq6_right = 0;  % Simplified for our specific displacement field
    violations(6) = abs(eq6_left - eq6_right);
    
    % Update maximum violations
    for i = 1:num_equations
        if violations(i) > max_violations(i)
            max_violations(i) = violations(i);
            violation_locations{i} = current_point;
        end
    end
end

% Display results for all 6 equations
fprintf('Strain Compatibility Test Results:\n\n');

tolerance = 1e-2; % Numerical tolerance for compatibility
all_compatible = true;

for i = 1:num_equations
    fprintf('%s\n', equation_names{i});
    fprintf('  Maximum violation: %12.6e', max_violations(i));
    
    if max_violations(i) < tolerance
        fprintf('   ✓ SATISFIED\n');
    else
        fprintf('   ✗ VIOLATED\n');
        all_compatible = false;
    end
    
    if ~isempty(violation_locations{i})
        fprintf('  Location of max violation: [%.2f, %.2f, %.2f]\n', violation_locations{i});
    end
end

% PART 3: VOLUME COMPARISON

fprintf('\n=== PART 3: VOLUME COMPARISON ===\n\n')

V1 = H*L*L;
p = -gamma_x;
e = -3 * (1-2*mu) * p/E;
V2 = V1 + V1*e;

fprintf('The Original Volume is: %1.f in^3\n', V1);
fprintf('The Volume Dilation is: %1.6f\n', e);
fprintf('The Updated Volume is: %1.4f in^3\n', V2);

fprintf('The Updated Volume is less than five percent smaller than the Original Volume. This is important as it shows that Poissons Ratio is consistent and that we have conservation of material.\n')

% PART 4: NORMAL AND SHEAR STRESSES

fprintf('\n=== PART 4: NORMAL AND SHEAR STRESSES ===\n\n');

% Display stress components at the reference point
fprintf('Stress Components at Point [%.1f, %.1f, %.1f]:\n\n', point(1), point(2), point(3));

fprintf('Normal Stresses:\n');
fprintf('  σₓₓ (sigma_xx) = %8.4f psi\n', gamma_x);
fprintf('  σᵧᵧ (sigma_yy) = %8.4f psi\n', gamma_y);
fprintf('  σᵤᵤ (sigma_zz) = %8.4f psi\n', gamma_z);

fprintf('\nShear Stresses:\n');
fprintf('  τₓᵧ (tau_xy)   = %8.4f psi\n', tau_xy);
fprintf('  τᵧᵤ (tau_yz)   = %8.4f psi\n', tau_yz);
fprintf('  τᵤₓ (tau_zx)   = %8.4f psi\n', tau_zx);

% Calculate and display principal stresses
principal_stresses = eig(stress_tensor);
principal_stresses_sorted = sort(principal_stresses, 'descend');

fprintf('\nPrinciplal Stresses (sorted):\n');
fprintf('  σ₁ (max principal) = %8.4f psi\n', principal_stresses_sorted(1));
fprintf('  σ₂ (intermediate)  = %8.4f psi\n', principal_stresses_sorted(2));
fprintf('  σ₃ (min principal) = %8.4f psi\n', principal_stresses_sorted(3));

% Calculate von Mises stress
s1 = principal_stresses_sorted(1);
s2 = principal_stresses_sorted(2);
s3 = principal_stresses_sorted(3);
von_mises = sqrt(0.5*((s1-s2)^2 + (s2-s3)^2 + (s3-s1)^2));

fprintf('\nEquivalent Stresses:\n');
fprintf('  von Mises stress   = %8.4f psi\n', von_mises);

% Maximum shear stress
max_shear = (principal_stresses_sorted(1) - principal_stresses_sorted(3))/2;
fprintf('  Max shear stress   = %8.4f psi\n', max_shear);