% Solid Mechanics Homework 1 - Enhanced with 3D Visualization
% Joshua Dillard
% Bridge bearing deformation and strain analysis
% Utilized Claude AI to convert base julia code to MATLAB script
% All equations and code for the z direction (w) are my own in Parts 1, 2,
% and 3
% Utilized Claude AI as the basis for developing the code for Parts 4 and 5
% All debugging was performed by me

% Clear workspace
clear; clc; close all;

% Define bridge bearing parameters
H = 12.0;   % Height in inches
L = 27.0;   % Width in inches
mu = 0.5;   % Poisson's ratio for natural rubber

% Displacement parameters
a = -2.0;           % Displacement at top of bearing, plane A-A-A-A
b = -mu * a;        % Maximum bulging or lateral displacement at midheight

% Part 1 Displacement Fields

% Define deformation field functions

% Vertical deformation field
v = @(point) a * point(2) / H;

% Lateral deformation field for the x direction  
u = @(point) b * sin(pi * point(2) / H) * (1 / (L/2)) * (point(1) - L/2);

% Lateral deformation field for the z direction  
w = @(point) b * sin(pi * point(2) / H) * (1 / (L/2)) * (point(3) - L/2);

fprintf('=== PART 1: === \n');

fprintf('\nDisplacement Field Equation for w: \n\n');

fprintf('w = b * sin(pi * y / H) * (1 / (L/2)) * (w - L/2) \n');

% Define the coordinates of the point of interest
point = [L, H/2, L]; % Adjust this code for any point desired

% Calculate deformations at this point (verification)
u_at_point = u(point);
v_at_point = v(point);
w_at_point = w(point);

fprintf('\nDeformation at point [%.1f, %.1f, %.1f]:\n', point(1), point(2), point(3));
fprintf('u = %.6f\n', u_at_point);
fprintf('v = %.6f\n', v_at_point);
fprintf('w = %.6f\n', w_at_point);

% Part 2 Gradients

fprintf('\n=== PART2: === \n');

fprintf('\nGeneral Gradients:\n\n');

% Gradient of u
du_dx = b * sin(pi * point(2) / H) * (1 / (L/2));
du_dy = b * cos(pi * point(2) / H) * (pi / H) * (1 / (L/2)) * (point(1) - L/2);
du_dz = 0;

% Gradient of v
dv_dx = 0;
dv_dy = a / H;
dv_dz = 0;

% Gradient of w
dw_dx = 0;
dw_dy = b * cos(pi * point(2) / H) * (pi / H) * (1 / (L/2)) * (point(3) - L/2);
dw_dz = b * sin(pi * point(2) / H) * (1 / (L/2));

grad_u = [du_dx, du_dy, du_dz];
grad_v = [dv_dx, dv_dy, dv_dz];
grad_w = [dw_dx, dw_dy, dw_dz];

fprintf('∇u = [∂u/∂x, ∂u/∂y, ∂u/∂z] = [b*sin(pi*y/H)*(1/(L/2)), b*cos(pi*y/H)*(pi/H)*(x-L/2)/(L/2), 0]\n');
fprintf('∇v = [∂v/∂x, ∂v/∂y, ∂v/∂z] = [0, a/H, 0]\n');
fprintf('∇w = [∂w/∂x, ∂w/∂y, ∂w/∂z] = [0, b*cos(pi*y/H)*(pi/H)*(z-L/2)/(L/2), b*sin(pi*y/H)/(L/2)]\n');

fprintf('\nSpecific Gradients:\n\n')

fprintf('∇u = [%.8f, %.8f, %.8f]\n', grad_u(1), grad_u(2), grad_u(3));
fprintf('∇v = [%.8f, %.8f, %.8f]\n', grad_v(1), grad_v(2), grad_v(3));
fprintf('∇w = [%.8f, %.8f, %.8f]\n', grad_w(1), grad_w(2), grad_w(3));

% Part 3 Strain Tensor

fprintf('\n=== PART 3: ===\n');

fprintf('\nStrain Tensor Components:\n\n');

fprintf('εxx = ∂u/∂x = b*sin(pi*y/H)/(L/2)\n');
fprintf('εyy = ∂v/∂y = a/H\n');
fprintf('εzz = ∂w/∂z = b*sin(pi*y/H)/(L/2)\n');
fprintf('γxy = γyx = (1/2)(∂u/∂y + ∂v/∂x) = (1/2)*b*cos(pi*y/H)*(pi/H)*(x-L/2)/(L/2)\n');
fprintf('γxz = γzx = (1/2)(∂u/∂z + ∂w/∂x) = 0\n');
fprintf('γyz = γzy = (1/2)(∂v/∂z + ∂w/∂y) = (1/2)*b*cos(pi*y/H)*(pi/H)*(z-L/2)/(L/2)\n');

% Alternative compact notation
fprintf('\nCompact matrix notation:\n\n');
fprintf('εij = \n');
fprintf('⎡ ε₁₁  γ₁₂  γ₁₃ ⎤   ⎡ εxx  γxy  γxz ⎤\n');
fprintf('⎢ γ₂₁  ε₂₂  γ₂₃ ⎥ = ⎢ γyx  εyy  γyz ⎥\n');
fprintf('⎣ γ₃₁  γ₃₂  ε₃₃ ⎦   ⎣ γzx  γzy  εzz ⎦\n');

% Original strain tensor calculation function
function strain_tensor = calculate_strain_tensor(point, H, L, a, b)
    % Calculate gradients at the given point
    
    % Gradient of u
    du_dx = b * sin(pi * point(2) / H) * (1 / (L/2));
    du_dy = b * cos(pi * point(2) / H) * (pi / H) * (1 / (L/2)) * (point(1) - L/2);
    du_dz = 0;
    
    % Gradient of v
    dv_dx = 0;
    dv_dy = a / H;
    dv_dz = 0;
    
    % Gradient of w
    dw_dx = 0;
    dw_dy = b * cos(pi * point(2) / H) * (pi / H) * (1 / (L/2)) * (point(3) - L/2);
    dw_dz = b * sin(pi * point(2) / H) * (1 / (L/2));
    
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

fprintf('\nStrain Tensor at Point [%.1f, %.1f, %.1f]:\n\n', point(1), point(2), point(3));

fprintf('[%12.8f  %12.8f  %12.8f]\n', strain_tensor(1,:));
fprintf('[%12.8f  %12.8f  %12.8f]\n', strain_tensor(2,:));
fprintf('[%12.8f  %12.8f  %12.8f]\n', strain_tensor(3,:));

% Part 4 and 5 Visual Representations

fprintf('\n=== PART 4 and 5: ===\n');

% Create 3D grid for visualization (rectangular prism)
nx = 20; ny = 15; nz = 20;  % Grid resolution
[X, Y, Z] = meshgrid(linspace(0, L, nx), linspace(0, H, ny), linspace(0, L, nz));

% Calculate deformed coordinates
X_def = zeros(size(X));
Y_def = zeros(size(Y));
Z_def = zeros(size(Z));

% Apply deformation field to each grid point
for i = 1:numel(X)
    point_orig = [X(i), Y(i), Z(i)];
    X_def(i) = X(i) + u(point_orig);
    Y_def(i) = Y(i) + v(point_orig);
    Z_def(i) = Z(i) + w(point_orig);
end

% Helper function to plot rectangular prism wireframe
function plot_prism_wireframe(X, Y, Z, color, linestyle, alpha)
    % Get the corner points of the rectangular prism
    corners_x = [0, L, L, 0, 0, L, L, 0];
    corners_y = [0, 0, 0, 0, H, H, H, H];
    corners_z = [0, 0, L, L, 0, 0, L, L];
    
    % Plot edges of the rectangular prism
    % Bottom face
    plot3([0,L,L,0,0], [0,0,0,0,0], [0,0,L,L,0], color, 'LineWidth', 2, 'LineStyle', linestyle);
    hold on;
    % Top face
    plot3([0,L,L,0,0], [H,H,H,H,H], [0,0,L,L,0], color, 'LineWidth', 2, 'LineStyle', linestyle);
    % Vertical edges
    plot3([0,0], [0,H], [0,0], color, 'LineWidth', 2, 'LineStyle', linestyle);
    plot3([L,L], [0,H], [0,0], color, 'LineWidth', 2, 'LineStyle', linestyle);
    plot3([L,L], [0,H], [L,L], color, 'LineWidth', 2, 'LineStyle', linestyle);
    plot3([0,0], [0,H], [L,L], color, 'LineWidth', 2, 'LineStyle', linestyle);
end

% 1. Plot 3D Deformed Shape with proper rectangular prism visualization
figure('Position', [100, 100, 1400, 600]);

% Original shape
subplot(1,2,1);
% Plot the rectangular grid as a wireframe to show the prism structure
for i = 1:nx
    for j = 1:ny
        plot3(squeeze(X(j,i,:)), squeeze(Y(j,i,:)), squeeze(Z(j,i,:)), 'b-', 'LineWidth', 0.5);
        hold on;
    end
end
for i = 1:nx
    for k = 1:nz
        plot3(squeeze(X(:,i,k)), squeeze(Y(:,i,k)), squeeze(Z(:,i,k)), 'b-', 'LineWidth', 0.5);
    end
end
for j = 1:ny
    for k = 1:nz
        plot3(squeeze(X(j,:,k)), squeeze(Y(j,:,k)), squeeze(Z(j,:,k)), 'b-', 'LineWidth', 0.5);
    end
end

title('Original Bridge Bearing (Rectangular Prism)');
xlabel('X (inches)'); ylabel('Y (inches)'); zlabel('Z (inches)');
axis equal; grid on;
view(-30, 30); % Adjusted view to show Y as vertical
xlim([-2 L+2]); ylim([-1 H+1]); zlim([-2 L+2]);

% Deformed shape
subplot(1,2,2);
% Plot the deformed grid
for i = 1:nx
    for j = 1:ny
        plot3(squeeze(X_def(j,i,:)), squeeze(Y_def(j,i,:)), squeeze(Z_def(j,i,:)), 'r-', 'LineWidth', 0.5);
        hold on;
    end
end
for i = 1:nx
    for k = 1:nz
        plot3(squeeze(X_def(:,i,k)), squeeze(Y_def(:,i,k)), squeeze(Z_def(:,i,k)), 'r-', 'LineWidth', 0.5);
    end
end
for j = 1:ny
    for k = 1:nz
        plot3(squeeze(X_def(j,:,k)), squeeze(Y_def(j,:,k)), squeeze(Z_def(j,:,k)), 'r-', 'LineWidth', 0.5);
    end
end

title('Deformed Bridge Bearing (2" compression, bulging)');
xlabel('X (inches)'); ylabel('Y (inches)'); zlabel('Z (inches)');
axis equal; grid on;
view(-30, 30); % Adjusted view to show Y as vertical
xlim([-2 L+2]); ylim([-3 H+1]); zlim([-2 L+2]);

% Overlay comparison plot showing bulging behavior
figure('Position', [100, 750, 1000, 700]);

% Plot original prism outline
plot3([0,L,L,0,0], [0,0,0,0,0], [0,0,L,L,0], 'b-', 'LineWidth', 3);
hold on;
plot3([0,L,L,0,0], [H,H,H,H,H], [0,0,L,L,0], 'b-', 'LineWidth', 3);
plot3([0,0], [0,H], [0,0], 'b-', 'LineWidth', 3);
plot3([L,L], [0,H], [0,0], 'b-', 'LineWidth', 3);
plot3([L,L], [0,H], [L,L], 'b-', 'LineWidth', 3);
plot3([0,0], [0,H], [L,L], 'b-', 'LineWidth', 3);

% Plot selected deformed grid lines to show bulging
step = 4; % Show every 4th line to avoid clutter
for i = 1:step:nx
    for j = 1:step:ny
        plot3(squeeze(X_def(j,i,:)), squeeze(Y_def(j,i,:)), squeeze(Z_def(j,i,:)), 'r-', 'LineWidth', 1);
    end
end
for j = 1:step:ny
    for k = 1:step:nz
        plot3(squeeze(X_def(j,:,k)), squeeze(Y_def(j,:,k)), squeeze(Z_def(j,:,k)), 'r-', 'LineWidth', 1);
    end
end

% Highlight the bulging at mid-height (y = H/2)
mid_y_idx = round(ny/2);
plot3(squeeze(X_def(mid_y_idx,:,:)), squeeze(Y_def(mid_y_idx,:,:)), squeeze(Z_def(mid_y_idx,:,:)), 'g-', 'LineWidth', 2);

% Show displacement vectors at key points
y_levels = [1, round(ny/4), round(ny/2), round(3*ny/4), ny];
colors = ['k', 'm', 'g', 'm', 'k'];
for level_idx = 1:length(y_levels)
    j = y_levels(level_idx);
    for i = 1:3:nx
        for k = 1:3:nz
            if norm([X_def(j,i,k)-X(j,i,k), Y_def(j,i,k)-Y(j,i,k), Z_def(j,i,k)-Z(j,i,k)]) > 0.01
                quiver3(X(j,i,k), Y(j,i,k), Z(j,i,k), ...
                       X_def(j,i,k)-X(j,i,k), Y_def(j,i,k)-Y(j,i,k), Z_def(j,i,k)-Z(j,i,k), ...
                       0, colors(level_idx), 'LineWidth', 1.5);
            end
        end
    end
end

legend('Original Prism', '', '', '', '', '', 'Deformed Grid', '', 'Mid-height Bulging', 'Location', 'best');
title('Bridge Bearing: Compression and Lateral Bulging');
xlabel('X (inches)'); ylabel('Y (inches)'); zlabel('Z (inches)');
axis equal; grid on;
view(-30, 30); % Adjusted view to show Y as vertical
xlim([-2 L+2]); ylim([-3 H+1]); zlim([-2 L+2]);

%% 2. Plot 3D Strain Field Contours
fprintf('\nCalculating strain fields...\n');

% Create a finer grid for strain visualization
nx_strain = 25; ny_strain = 20; nz_strain = 25;
[Xs, Ys, Zs] = meshgrid(linspace(0, L, nx_strain), linspace(0, H, ny_strain), linspace(0, L, nz_strain));

% Calculate strain components at each grid point
epsilon_xx_field = zeros(size(Xs));
epsilon_yy_field = zeros(size(Ys));
epsilon_zz_field = zeros(size(Zs));
epsilon_xy_field = zeros(size(Xs));
epsilon_yz_field = zeros(size(Ys));
von_mises_field = zeros(size(Xs));

for i = 1:numel(Xs)
    point_strain = [Xs(i), Ys(i), Zs(i)];
    strain_tensor_local = calculate_strain_tensor(point_strain, H, L, a, b);
    
    epsilon_xx_field(i) = strain_tensor_local(1,1);
    epsilon_yy_field(i) = strain_tensor_local(2,2);
    epsilon_zz_field(i) = strain_tensor_local(3,3);
    epsilon_xy_field(i) = strain_tensor_local(1,2);
    epsilon_yz_field(i) = strain_tensor_local(2,3);
    
    % Calculate von Mises equivalent strain
    exx = strain_tensor_local(1,1);
    eyy = strain_tensor_local(2,2);
    ezz = strain_tensor_local(3,3);
    exy = strain_tensor_local(1,2);
    exz = strain_tensor_local(1,3);
    eyz = strain_tensor_local(2,3);
    
    von_mises_field(i) = sqrt(2/3) * sqrt((exx-eyy)^2 + (eyy-ezz)^2 + (ezz-exx)^2 + ...
                                         6*(exy^2 + exz^2 + eyz^2));
end

% Plot strain field contours with proper Y-vertical orientation
figure('Position', [200, 200, 1400, 1000]);

% Normal strain εxx (lateral bulging strain)
subplot(3,3,1);
slice(Xs, Ys, Zs, epsilon_xx_field, [L/4, 3*L/4], [H/4, H/2, 3*H/4], L/2);
colorbar; colormap(jet);
title('Normal Strain ε_{xx} (X-direction bulging)');
xlabel('X'); ylabel('Y'); zlabel('Z');
shading interp;
view(-30, 30); % Y-axis vertical view

% Normal strain εyy (compression strain)
subplot(3,3,2);
slice(Xs, Ys, Zs, epsilon_yy_field, L/2, [H/4, H/2, 3*H/4], L/2);
colorbar; colormap(jet);
title('Normal Strain ε_{yy} (Y-direction compression)');
xlabel('X'); ylabel('Y - VERTICAL'); zlabel('Z');
shading interp;
view(-30, 30); % Y-axis vertical view

% Normal strain εzz (lateral bulging strain)
subplot(3,3,3);
slice(Xs, Ys, Zs, epsilon_zz_field, L/2, [H/4, H/2, 3*H/4], [L/4, 3*L/4]);
colorbar; colormap(jet);
title('Normal Strain ε_{zz} (Z-direction bulging)');
xlabel('X'); ylabel('Y'); zlabel('Z');
shading interp;
view(-30, 30); % Y-axis vertical view

% Shear strain εxy
subplot(3,3,4);
slice(Xs, Ys, Zs, epsilon_xy_field, [L/4, 3*L/4], H/2, L/2);
colorbar; colormap(jet);
title('Shear Strain ε_{xy}');
xlabel('X'); ylabel('Y'); zlabel('Z');
shading interp;
view(-30, 30); % Y-axis vertical view

% Shear strain εyz
subplot(3,3,5);
slice(Xs, Ys, Zs, epsilon_yz_field, L/2, H/2, [L/4, 3*L/4]);
colorbar; colormap(jet);
title('Shear Strain ε_{yz}');
xlabel('X'); ylabel('Y'); zlabel('Z');
shading interp;
view(-30, 30); % Y-axis vertical view

% von Mises equivalent strain
subplot(3,3,6);
slice(Xs, Ys, Zs, von_mises_field, L/2, [H/4, H/2, 3*H/4], L/2);
colorbar; colormap(jet);
title('von Mises Equivalent Strain');
xlabel('X'); ylabel('Y'); zlabel('Z');
shading interp;
view(-30, 30); % Y-axis vertical view

% Cross-sections showing strain variation
% XY plane at mid-Z
subplot(3,3,7);
contourf(squeeze(Xs(:,:,round(nz_strain/2))), squeeze(Ys(:,:,round(nz_strain/2))), ...
         squeeze(von_mises_field(:,:,round(nz_strain/2))), 20);
colorbar; colormap(jet);
title('von Mises Strain - XY plane (mid-Z)');
xlabel('X'); ylabel('Y');
axis equal;

% YZ plane at mid-X
subplot(3,3,8);
contourf(squeeze(Ys(round(ny_strain/2),:,:))', squeeze(Zs(round(ny_strain/2),:,:))', ...
         squeeze(von_mises_field(round(ny_strain/2),:,:))', 20);
colorbar; colormap(jet);
title('von Mises Strain - YZ plane (mid-X)');
xlabel('Y'); ylabel('Z');
axis equal;

% XZ plane at mid-Y (maximum bulging)
subplot(3,3,9);
contourf(squeeze(Xs(round(ny_strain/2),:,:)), squeeze(Zs(round(ny_strain/2),:,:)), ...
         squeeze(von_mises_field(round(ny_strain/2),:,:)), 20);
colorbar; colormap(jet);
title('von Mises Strain - XZ plane (mid-Y, max bulging)');
xlabel('X'); ylabel('Z');
axis equal;

% Additional 3D strain visualization showing bulging pattern
figure('Position', [300, 300, 1000, 800]);

% Create isosurfaces for different strain levels
max_strain = max(von_mises_field(:));
isovalues = [0.15, 0.3, 0.5, 0.7] * max_strain;
colors = {'green', 'blue', 'yellow', 'red'};
alphas = [0.2, 0.3, 0.4, 0.5];

for i = 1:length(isovalues)
    if any(von_mises_field(:) >= isovalues(i))
        p = patch(isosurface(Xs, Ys, Zs, von_mises_field, isovalues(i)));
        isonormals(Xs, Ys, Zs, von_mises_field, p);
        p.FaceColor = colors{i};
        p.EdgeColor = 'none';
        p.FaceAlpha = alphas(i);
        hold on;
    end
end

% Add wireframe of original prism for reference
plot3([0,L,L,0,0], [0,0,0,0,0], [0,0,L,L,0], 'k--', 'LineWidth', 2);
plot3([0,L,L,0,0], [H,H,H,H,H], [0,0,L,L,0], 'k--', 'LineWidth', 2);
plot3([0,0], [0,H], [0,0], 'k--', 'LineWidth', 2);
plot3([L,L], [0,H], [0,0], 'k--', 'LineWidth', 2);
plot3([L,L], [0,H], [L,L], 'k--', 'LineWidth', 2);
plot3([0,0], [0,H], [L,L], 'k--', 'LineWidth', 2);

title('3D von Mises Strain Distribution - Bridge Bearing Bulging');
xlabel('X (inches)'); ylabel('Y (inches)'); zlabel('Z (inches)');
lighting gouraud;
camlight;
axis equal;
grid on;
view(-30, 30); % Y-axis vertical view
xlim([0 L]); ylim([0 H]); zlim([0 L]);

% Add text annotations
text(L/2, H+1, L/2, 'Maximum bulging at mid-height (Y = H/2)', 'HorizontalAlignment', 'center');
text(L+2, H/2, L/2, 'Lateral bulging', 'HorizontalAlignment', 'center');

fprintf('\nVisualization complete!\n');
fprintf('- Figure 1: Original vs Deformed rectangular prism comparison\n');
fprintf('- Figure 2: Deformation with lateral bulging and displacement vectors\n');
fprintf('- Figure 3: Complete strain field analysis with cross-sections\n');
fprintf('- Figure 4: 3D strain isosurfaces showing bulging pattern\n');
fprintf('\nKey observations:\n');
fprintf('- Maximum lateral bulging occurs at mid-height (Y = H/2 = 6 inches)\n');
fprintf('- Compression strain (εyy) is constant = %.4f\n', a/H);
fprintf('- Lateral bulging strains (εxx, εzz) vary sinusoidally with height\n');
fprintf('- Maximum bulging strain = %.4f at mid-height edges\n', b/(L/2));