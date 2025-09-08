function plot_surface_strain_prof(params, varargin)
% plot_surface_strain: strain "cloud" on 6 deformed outer faces

% ---------- parse options ----------
p = inputParser;
addParameter(p,'N',41,@(x)isnumeric(x)&&x>=5);
addParameter(p,'Scale',1.0,@(x)isnumeric(x)&&x>0);
addParameter(p,'Field','eq',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
N   = p.Results.N;
SCL = p.Results.Scale;
field = lower(string(p.Results.Field));

Lx = params.Lx;  H = params.H;  Lz = params.Lz;

% ---------- grid on a box (for 6 faces only) ----------
xv = linspace(0,Lx,N);
yv = linspace(0,H, N);
zv = linspace(0,Lz,N);
[X, Y, Z] = meshgrid(xv, yv, zv);

faces = { ...
    struct('X', X(:,:,1),   'Y', Y(:,:,1),   'Z', Z(:,:,1),   'name','z=0'), ...
    struct('X', X(:,:,end), 'Y', Y(:,:,end), 'Z', Z(:,:,end), 'name','z=Lz'), ...
    struct('X', squeeze(X(1,:,:))',   'Y', squeeze(Y(1,:,:))',   'Z', squeeze(Z(1,:,:))',   'name','x=0'), ...
    struct('X', squeeze(X(end,:,:))', 'Y', squeeze(Y(end,:,:))', 'Z', squeeze(Z(end,:,:))', 'name','x=Lx'), ...
    struct('X', squeeze(X(:,1,:))',   'Y', squeeze(Y(:,1,:))',   'Z', squeeze(Z(:,1,:))',   'name','y=0'), ...
    struct('X', squeeze(X(:,end,:))', 'Y', squeeze(Y(:,end,:))', 'Z', squeeze(Z(:,end,:))', 'name','y=H') ...
    };

% ---------- first pass: compute scalar range for unified colorbar ----------
vals = cell(1,numel(faces));
vmin = inf; vmax = -inf;

for k=1:numel(faces)
    Xf = faces{k}.X; Yf = faces{k}.Y; Zf = faces{k}.Z;
    S  = strain_tensor(Xf, Yf, Zf, params);
    vals{k} = pick_scalar(S, field);
    vmin = min(vmin, min(vals{k}(:)));
    vmax = max(vmax, max(vals{k}(:)));
end
if ~isfinite(vmax - vmin) || (vmax - vmin) <= 1e-12
    cmin = vmin - 1e-12; cmax = vmax + 1e-12;
else
    cmin = vmin; cmax = vmax;
end

% ---------- draw 6 deformed faces colored by strain ----------
figure('Color','w'); hold on;
for k=1:numel(faces)
    Xf = faces{k}.X; Yf = faces{k}.Y; Zf = faces{k}.Z;

    F = uvw_model(Xf, Yf, Zf, params);

    Xd = Xf + SCL*F.u;    
    Yd = Zf + SCL*F.w;    
    Zd = Yf + SCL*F.v;    

    C  = vals{k};         
    surf(Xd, Yd, Zd, C, 'EdgeColor','none', 'FaceAlpha', 0.95);
end

axis equal; box on; view(45,22);
xlabel('x (m)'); ylabel('z (m)'); zlabel('y (m)');

ttl = "Surface strain cloud — "+ field_label(field) + ...
      sprintf("  (Scale = %.2f)", SCL);
title(ttl, 'Interpreter','tex');

colormap(parula);
cb = colorbar('eastoutside');
cb.Label.String = field_label(field);
caxis([cmin cmax]);

plot_undeformed_box_xzy(Lx, Lz, H);
end

% ---------- helpers ----------
function val = pick_scalar(S, field)
switch field
    case "eq"   % Frobenius norm of small strain tensor
        val = sqrt(S.eps_xx.^2 + S.eps_yy.^2 + S.eps_zz.^2 + ...
                   2*(S.eps_xy.^2 + S.eps_xz.^2 + S.eps_yz.^2));
    case "exx", val = S.eps_xx;
    case "eyy", val = S.eps_yy;
    case "ezz", val = S.eps_zz;
    case "exy", val = S.eps_xy;
    case "exz", val = S.eps_xz;
    case "eyz", val = S.eps_yz;
    otherwise, error('Unknown Field: %s', field);
end
end

function s = field_label(field)
switch field
    case "eq",  s = '||\epsilon||_F';
    case "exx", s = '\epsilon_{xx}';
    case "eyy", s = '\epsilon_{yy}';
    case "ezz", s = '\epsilon_{zz}';
    case "exy", s = '\epsilon_{xy}';
    case "exz", s = '\epsilon_{xz}';
    case "eyz", s = '\epsilon_{yz}';
end
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
