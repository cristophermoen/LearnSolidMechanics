function plot_deformed_bearing(params, varargin)
% plot_deformed_bearing: 3D deformed shape of the bearing
%
% Inputs:
%   params.Lx, params.Ly, params.H : dimensions (m)
%   params.Delta : vertical compression at z=H (m)
%   params.A     : vertical bulging amplitude (m)
%   params.B     : lateral amplitude (m)
%
% Name-Value:
%   'N'      : grid nodes per edge
%   'Scale'  : displacement scale factor
%
% Convention:
%   The origin is at the center of the bottom surface; the Z axis is upward; z∈[0,H]
%
% Dependencies: uv_model.m, w_and_grad.m

% ---- parse options ----
p = inputParser;
addParameter(p, 'N', 41, @(x)isnumeric(x)&&x>=5);
addParameter(p, 'Scale', 1.0, @(x)isnumeric(x)&&x>0);
parse(p, varargin{:});
N = p.Results.N;
SCL = p.Results.Scale;

Lx = params.Lx; Ly = params.Ly; H = params.H;

% ---- build a hexahedral grid ----
xv = linspace(-Lx/2, Lx/2, N);
yv = linspace(-Ly/2, Ly/2, N);
zv = linspace(0, H, N);
[X, Y, Z] = meshgrid(xv, yv, zv);  % full grid

% surfaces only:
faces = { ...
    struct('X', X(:,:,1),   'Y', Y(:,:,1),   'Z', Z(:,:,1),   'name','Bottom'), ...
    struct('X', X(:,:,end), 'Y', Y(:,:,end), 'Z', Z(:,:,end), 'name','Top'), ...
    struct('X', squeeze(X(1,:,:))',   'Y', squeeze(Y(1,:,:))',   'Z', squeeze(Z(1,:,:))',   'name','Y- Left'), ...
    struct('X', squeeze(X(end,:,:))', 'Y', squeeze(Y(end,:,:))', 'Z', squeeze(Z(end,:,:))', 'name','Y- Right'), ...
    struct('X', squeeze(X(:,1,:))',   'Y', squeeze(Y(:,1,:))',   'Z', squeeze(Z(:,1,:))',   'name','X- Front'), ...
    struct('X', squeeze(X(:,end,:))', 'Y', squeeze(Y(:,end,:))', 'Z', squeeze(Z(:,end,:))', 'name','X- Back') ...
    };

% ---- compute displacements on each face ----
figure('Color','w'); hold on;
for k = 1:numel(faces)
    Xf = faces{k}.X; Yf = faces{k}.Y; Zf = faces{k}.Z;
    % lateral displacements (u,v) 
    [u,v,~,~,~,~,~,~] = uv_model(Xf, Yf, Zf, params);
    % vertical displacement (w) 
    [w,~,~,~] = w_and_grad(Xf, Yf, Zf, params);

    % apply scale factor 
    Xd = Xf + SCL*u;  Yd = Yf + SCL*v;  Zd = Zf + SCL*w;

    % draw surface 
    surf(Xd, Yd, Zd, Zf, 'EdgeColor',[0.3 0.3 0.3], 'FaceAlpha',0.85);
end

% ---- plot settings ----
axis equal; box on; view(40,25);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Deformed shape (Scale = %.2f)', SCL));

% draw the undeformed box
plot_undeformed_box(Lx, Ly, H);

colormap(parula); colorbar('eastoutside'); % Color by original z to indicate height
end

function plot_undeformed_box(Lx,Ly,H)
% Auxiliary: Draw the undeformed frame (dashed line)
cx = [-Lx/2 Lx/2];
cy = [-Ly/2 Ly/2];
cz = [0 H];
P = [cx(1) cy(1) cz(1);
     cx(2) cy(1) cz(1);
     cx(2) cy(2) cz(1);
     cx(1) cy(2) cz(1);
     cx(1) cy(1) cz(2);
     cx(2) cy(1) cz(2);
     cx(2) cy(2) cz(2);
     cx(1) cy(2) cz(2)];
E = [1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
for i=1:size(E,1)
    plot3(P(E(i,1),1),P(E(i,1),2),P(E(i,1),3), 'k--','LineWidth',1.0);
    plot3(P(E(i,2),1),P(E(i,2),2),P(E(i,2),3), 'k--','LineWidth',1.0);
end
end
