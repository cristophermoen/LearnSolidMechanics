function plot_deformed_bearing_prof(params, varargin)
% plot_deformed_bearing_prof: 3D deformed shape (Professor-based field)
%

% ---- parse options ----
p = inputParser;
addParameter(p, 'N', 41, @(x)isnumeric(x)&&x>=5);
addParameter(p, 'Scale', 1.0, @(x)isnumeric(x)&&x>0);
parse(p, varargin{:});
N   = p.Results.N;
SCL = p.Results.Scale;

Lx = params.Lx; H = params.H; Lz = params.Lz;

% ---- build a hexahedral grid ----
xv = linspace(0, Lx, N);
yv = linspace(0, H,  N);
zv = linspace(0, Lz, N);
[X, Y, Z] = meshgrid(xv, yv, zv);   

% Draw only exterior surfaces (6 faces)
faces = { ...
    struct('X', X(:,:,1),   'Y', Y(:,:,1),   'Z', Z(:,:,1),   'name','z = 0 (front)'), ...
    struct('X', X(:,:,end), 'Y', Y(:,:,end), 'Z', Z(:,:,end), 'name','z = Lz (back)'), ...
    struct('X', squeeze(X(1,:,:))',   'Y', squeeze(Y(1,:,:))',   'Z', squeeze(Z(1,:,:))',   'name','x = 0 (left)'), ...
    struct('X', squeeze(X(end,:,:))', 'Y', squeeze(Y(end,:,:))', 'Z', squeeze(Z(end,:,:))', 'name','x = Lx (right)'), ...
    struct('X', squeeze(X(:,1,:))',   'Y', squeeze(Y(:,1,:))',   'Z', squeeze(Z(:,1,:))',   'name','y = 0 (bottom)'), ...
    struct('X', squeeze(X(:,end,:))', 'Y', squeeze(Y(:,end,:))', 'Z', squeeze(Z(:,end,:))', 'name','y = H (top)') ...
    };

% ---- compute displacements on each face ----
figure('Color','w'); hold on;
for k = 1:numel(faces)
    Xf = faces{k}.X; Yf = faces{k}.Y; Zf = faces{k}.Z;

    
    F = uvw_model(Xf, Yf, Zf, params);
    u = F.u; v = F.v; w = F.w;

    
    Xd = Xf + SCL*u;      
    Yd = Zf + SCL*w;      
    Zd = Yf + SCL*v;      

    
    C  = Yf;

   
    surf(Xd, Yd, Zd, C, 'EdgeColor',[0.3 0.3 0.3], 'FaceAlpha',0.9);
end

% ---- plot settings ----
axis equal; box on; view(45,22);
xlabel('x (m)'); ylabel('z (m)'); zlabel('y (m)'); 
title(sprintf('Deformed shape (Professor-based, Scale = %.2f)', SCL));

plot_undeformed_box_xzy(Lx, Lz, H);

colormap(parula);
cb = colorbar('eastoutside'); cb.Label.String = 'original y (m)';
end

% --------- helper: undeformed box in (x,z,y) order ---------
function plot_undeformed_box_xzy(Lx, Lz, H)
cx = [0 Lx]; cz = [0 Lz]; cy = [0 H];
P = [cx(1) cz(1) cy(1);
     cx(2) cz(1) cy(1);
     cx(2) cz(2) cy(1);
     cx(1) cz(2) cy(1);
     cx(1) cz(1) cy(2);
     cx(2) cz(1) cy(2);
     cx(2) cz(2) cy(2);
     cx(1) cz(2) cy(2)];
E = [1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
for i=1:size(E,1)
    plot3(P(E(i,1),1),P(E(i,1),2),P(E(i,1),3), 'k--','LineWidth',1.0);
    plot3(P(E(i,2),1),P(E(i,2),2),P(E(i,2),3), 'k--','LineWidth',1.0);
end
end
