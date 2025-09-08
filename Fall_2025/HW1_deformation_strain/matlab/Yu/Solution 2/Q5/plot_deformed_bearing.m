function plot_deformed_bearing_prof(params, varargin)
% plot_deformed_bearing_prof: 3D deformed shape (Professor-based field)
% 绘制三维变形外形（方法二：教授笔记版；原点=角点；竖直=y）
%
% Inputs:
%   params.Lx, params.H, params.Lz : 尺寸 / dimensions (m)
%   params.Delta : 顶面竖向压缩量 Δ (m)
%   params.A     : 鼓曲幅值 (m)
%
% Name-Value:
%   'N'     : 网格密度（默认 41）
%   'Scale' : 位移放大系数（默认 1.0）
%
% 坐标约定：x∈[0,Lx], y∈[0,H](竖直), z∈[0,Lz]
% 绘图时用 (x+u, z+w, y+v) 作为 (X,Y,Z)，使 y 在图上为竖直轴

p = inputParser;
addParameter(p, 'N', 41, @(x)isnumeric(x)&&x>=5);
addParameter(p, 'Scale', 1.0, @(x)isnumeric(x)&&x>0);
parse(p, varargin{:});
N   = p.Results.N;
SCL = p.Results.Scale;

Lx = params.Lx; H = params.H; Lz = params.Lz;

xv = linspace(0, Lx, N);
yv = linspace(0, H,  N);
zv = linspace(0, Lz, N);
[X, Y, Z] = meshgrid(xv, yv, zv);   % X=x, Y=y(vert), Z=z

faces = { ...
    struct('X', X(:,:,1),   'Y', Y(:,:,1),   'Z', Z(:,:,1)), ...
    struct('X', X(:,:,end), 'Y', Y(:,:,end), 'Z', Z(:,:,end)), ...
    struct('X', squeeze(X(1,:,:))',   'Y', squeeze(Y(1,:,:))',   'Z', squeeze(Z(1,:,:))'), ...
    struct('X', squeeze(X(end,:,:))', 'Y', squeeze(Y(end,:,:))', 'Z', squeeze(Z(end,:,:))'), ...
    struct('X', squeeze(X(:,1,:))',   'Y', squeeze(Y(:,1,:))',   'Z', squeeze(Z(:,1,:))'), ...
    struct('X', squeeze(X(:,end,:))', 'Y', squeeze(Y(:,end,:))', 'Z', squeeze(Z(:,end,:))') ...
    };

figure('Color','w'); hold on;
for k = 1:numel(faces)
    Xf = faces{k}.X; Yf = faces{k}.Y; Zf = faces{k}.Z;
    F = uvw_prof_model(Xf, Yf, Zf, params);
    Xd = Xf + SCL*F.u;      % x + u
    Yd = Zf + SCL*F.w;      % 图上第二轴用 z + w
    Zd = Yf + SCL*F.v;      % 图上第三轴用 y + v (竖直)
    surf(Xd, Yd, Zd, Yf, 'EdgeColor',[0.3 0.3 0.3], 'FaceAlpha',0.9);
end

axis equal; box on; view(45,22);
xlabel('x (m)'); ylabel('z (m)'); zlabel('y (m)');
title(sprintf('Deformed shape (Professor-based, Scale = %.2f)', SCL));
plot_undeformed_box_xzy(Lx, Lz, H);
colorbar('eastoutside'); colormap(parula);
end

function plot_undeformed_box_xzy(Lx, Lz, H)
cx=[0 Lx]; cz=[0 Lz]; cy=[0 H];
P=[cx(1) cz(1) cy(1); cx(2) cz(1) cy(1); cx(2) cz(2) cy(1); cx(1) cz(2) cy(1); ...
   cx(1) cz(1) cy(2); cx(2) cz(1) cy(2); cx(2) cz(2) cy(2); cx(1) cz(2) cy(2)];
E=[1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
for i=1:size(E,1)
    plot3(P(E(i,1),1),P(E(i,1),2),P(E(i,1),3),'k--','LineWidth',1.0);
    plot3(P(E(i,2),1),P(E(i,2),2),P(E(i,2),3),'k--','LineWidth',1.0);
end
end
