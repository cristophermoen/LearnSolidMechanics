% demo_q4_prof.m
% Q4 – 3D deformed shape (Professor-based, corner origin, y vertical)

clear; close all; clc;

% -------- Parameters --------
params.Lx = 0.10;    % x-length (m)
params.H  = 0.20;    % thickness, vertical y (m)
params.Lz = 0.10;    % z-length (m)

inch = 0.0254;                 % 1 inch in meters
params.Delta = 2.0 * inch;     % vertical compression at y=H

params.A = 0.02;               % bulging amplitude (m)

% -------- Plot deformed shape --------
plot_deformed_bearing(params, 'N', 61, 'Scale', 1.0);

% plot_deformed_bearing_prof(params, 'N', 61, 'Scale', 3.0);
