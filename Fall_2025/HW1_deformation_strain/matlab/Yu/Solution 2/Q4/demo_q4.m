clear; close all; clc;

% -------- Parameters --------
params.Lx = 0.10;           % m
params.Ly = 0.10;           % m
params.H  = 0.10;           % m

inch = 0.0254;              % 1 inch in meters
params.Delta = 2.0*inch;    % 2.0 inches vertical compression

% Drum amplitude
params.A  = 0.02;          % m  vertical bulging amplitude
params.B  = 0.01;          % m  lateral amplitude

% -------- Plot deformed shape --------
% N: grid density (default 41); Scale: displacement magnification factor
plot_deformed_bearing(params, 'N', 61, 'Scale', 1.0);
