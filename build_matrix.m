function [MATRIX, wn_k, lat, lat_deg, z, ny, dt, second] = build_matrix(wavenumber_factor, opt)

if nargin < 2 || isempty(opt); opt = struct(); end
opt = fill_defaults(opt);

%% build matrix for the EVP with given parameters
%current parameters: wavenumber_factor, dlat, hyperdiff_scale_step

%define units
meter = 1.;
second = 1.;
kg = 1.;
gram = 1e-3 * kg;
kelvin = 1.;
pa = kg / meter / second^2;
joule = kg * meter^2 / second^2;
%define physical constants
ggr=9.81 * meter / second^2;
Cp=1004 * joule /kg/ kelvin;
RE=6370e3 * meter; % radius of the Earth
beta=4*pi/(86400.*second)/RE;
wn_k=wavenumber_factor/RE;
%define model hyperparameters
nT=26; % 26 layers for T
nQ=14; % 14 layers for Q
dt=900.*second; % fixed time step

%% 1. initialize the reference profiles and grid etc.
%name of the netcdf file
nc='IDEAL2_2048x1024x28_2048_f0_rep3.nc';
% 1.1 read out the basic states
t_wavebg=ncread(nc,'T_WAVEBG');
q_wavebg=ncread(nc,'Q_WAVEBG');
t_wavebg=t_wavebg(:,1) * kelvin;
q_wavebg=q_wavebg(:,1) * gram / kg;
p=ncread(nc,'p')';
pres=ncread(nc,'p')';pres=pres'*(100*pa);
z=ncread(nc,'z') * meter;
nzm=numel(z);nz=nzm+1;

% 1.2 density at interface levels in kg/m3
for k=2:nzm
  rhow(k) =  (pres(k-1)-pres(k))/(z(k)-z(k-1))/ggr;
end
rhow(1) = 2*rhow(2) - rhow(3);
rhow(nz)= 2*rhow(nzm) - rhow(nzm-1);
rhow=rhow'; %make it a column vector
%find interface height
dz = 0.5*(z(1)+z(2));
for k=2:nzm
   adzw(k) = (z(k)-z(k-1))/dz;
end
adzw(1) = 1.;
adzw(nz) = adzw(nzm);
adz(1) = 1.;
for k=2:nzm-1
   adz(k) = 0.5*(z(k+1)-z(k-1))/dz;
end
adz(nzm) = adzw(nzm);
zi(1) = 0.;
for k=2:nz
   zi(k) = zi(k-1) + adz(k-1)*dz;
end
%interpolate layer density
rho=interp1(zi,rhow,z);

% 1.3 read mass weight in each layer
tmp=load('IDEAL_dm.mat');
dm=double(tmp.dm);
%weight variance by mass. Works best for vertically coherent signal in
%that subdividing a layer into multiple layers that are coherent does
%not change the answer.
mw_vec = sqrt([dm(1:26); (2.5^2)*dm(1:14)]);
mass_weight = spdiags(mw_vec, 0, numel(mw_vec), numel(mw_vec));

% 1.4 meridional grids
% test_bc: Neumann for all variables
% boundary value is now on the face centers, lat are cell center
max_lat = opt.max_lat; % degrees of boundary latitude for slat
dlat = opt.dlat; % degrees of grid spacing
slat_deg = (-max_lat:dlat:max_lat)';                 % faces (ny+1)
lat_deg  = 0.5*(slat_deg(1:end-1) + slat_deg(2:end)); % centers (ny)
ny = numel(lat_deg);
% interior faces for staggered v (Dirichlet v=0 at boundaries)
slat_int_deg = slat_deg(2:end-1);                    % (ny-1)
Nv = numel(slat_int_deg);                            % = ny-1

lat = lat_deg*pi/180*RE;      % meters (centers)
slat = slat_deg*pi/180*RE;    % meters (faces)
slat_int = slat(2:end-1);     % meters (interior faces)
dy = slat(2) - slat(1);        % uniform grid spacing in meters

%% 2. load state space model parameters
tmp=load('sys_randx2_40_120_free_prediction_estx0_nomean_tabs','sys');
sys=tmp.sys;
A=sys.A;B=sys.B*mass_weight;C=mass_weight\sys.C;
for i=1:26
    B(:,i) = B(:,i)/kelvin;
    C(i,:) = C(i,:)*kelvin;
end
for i=27:40
    B(:,i) = B(:,i)/(gram/kg);
    C(i,:) = C(i,:)*(gram/kg);
end

%% 3. set optional damping and diffusion
% 3.1 damp U, V, T everywhere
damping = opt.damping.enable;
epsdt   = (dt/(opt.damping.tau_day*86400)) * double(damping);

% 3.2 damp U, V at higher latitudes
damping_bound      = opt.sponge.enable;
damping_bound_lat  = opt.sponge.lat0_deg;           % degrees: sponge starts at |lat| >= this
damping_bound_time = opt.sponge.tau_hr*3600*second; % seconds: e-folding time at the boundary (e.g. 7200 like test_bc)
ramp_power         = opt.sponge.ramp_power;         % 1=linear, 2=quadratic

% 3.3 add diffusion to U, V in zonal and meridional direction
diffusion=opt.diffusion.enable;
nu_visc = opt.diffusion.nu_visc * meter^2 / second;

hyperdiff_scale_step = opt.hyperdiff.scale_step;
hyperdiff = (hyperdiff_scale_step < 1.0);
if hyperdiff
    nu4 = -log(hyperdiff_scale_step)/dt * (dy/pi)^4;   % m^4/s
else
    nu4 = 0.0;
end

%% 4. build the matrix
% vector        = [ x,   d(rhou)dz, d(rhov)dz(staggered)];
% vector length = ny*120 + ny*26 + (ny-1)*26
%          | x->x,         d(rhou)dz->x,         d(rhov)dz->x         |
% matrix = | x->d(rhou)dz, d(rhou)dz->d(rhou)dz, d(rhou)dz->d(rhov)dz |
%          | x->d(rhov)dz, d(rhou)dz->d(rhov)dz, d(rhov)dz->d(rhov)dz |
% assuming zero signal in higher latitudes than max_lat
disp(['Start building matrix for ' num2str(wavenumber_factor)])

% 4.1 auxiliary matrices (sparse)
I26 = speye(26);
I120 = speye(120);
% 1D meridional operators (sparse, BC-update)
% lat-grid: ny centers; slat-grid: Nv=ny-1 interior faces
e_ny = ones(ny,1);
e_nv = ones(Nv,1);

% second derivative on centers (Neumann-like at boundaries)
Lyy_lat_single = spdiags([e_ny -2*e_ny e_ny], -1:1, ny, ny) / (dy^2);
%ghost points has same value, so update BC for Lyy_lat_single
Lyy_lat_single(1,1) = -1 / (dy^2);
Lyy_lat_single(ny, ny) = -1 / (dy^2);

% second derivative on interior faces
Lyy_slat_single = spdiags([e_nv -2*e_nv e_nv], -1:1, Nv, Nv) / (dy^2);
%ghost points has same value, so update BC for Lyy_slat_single
Lyy_slat_single(1,1) = -1 / (dy^2);
Lyy_slat_single(Nv, Nv) = -1 / (dy^2);

% first derivatives
Dy_s2l_single = spdiags([-e_nv  e_nv], [-1 0], ny, Nv) / dy;
% if want to enforce zero signal at boundaries, then update BC for Dy_s2l_single
Dy_s2l_single(1,1) = 0.;
Dy_s2l_single(ny, Nv) = 0.;

Dy_l2s_single = spdiags([-e_nv  e_nv], [0 1], Nv, ny) / dy;

% interpolation between grids (for Coriolis coupling)
slat_to_lat_single = spdiags([0.5*e_nv, 0.5*e_nv], [-1 0], ny, Nv);
% if want to enforce zero gradient at boundaries, then update BC for slat_to_lat_single
slat_to_lat_single(1,1) = 1.;
slat_to_lat_single(ny, Nv) = 1.;

lat_to_slat_single = spdiags([0.5*e_nv 0.5*e_nv], [0 1], Nv, ny);

% lift to full (kron with I26)
partial2_y_lat_to_lat   = kron(Lyy_lat_single,  I26);
partial2_y_slat_to_slat = kron(Lyy_slat_single, I26);
partial_y_slat_to_lat   = kron(Dy_s2l_single,   I26);
partial_y_lat_to_slat   = kron(Dy_l2s_single,   I26);

slat_to_lat = kron(slat_to_lat_single, I26);
lat_to_slat = kron(lat_to_slat_single, I26);

% nu4 hyperdiffusion (sparse): (d^2/dy^2)^2
L4_lat_single  = Lyy_lat_single  * Lyy_lat_single;
L4_slat_single = Lyy_slat_single * Lyy_slat_single;
hdiff_u = hyperdiff * (nu4*dt) * kron(L4_lat_single,  I26);   % (ny*26)x(ny*26)
hdiff_v = hyperdiff * (nu4*dt) * kron(L4_slat_single, I26);   % (Nv*26)x(Nv*26)

%compute w from dudz, dvdz (divergence on lat grid)
m_to_int_u = -1j*wn_k*speye(26*ny);
m_to_int_v = -partial_y_slat_to_lat;

% build vertical laplacian (sparse) and integration operator
adz_v = zeros(nT+1,1);
adz_v(2:nT) = z(2:nT) - z(1:nT-1);
adz_v(1)    = 2*z(1);
adz_v(nT+1) = z(nT+1) - z(nT);
aa = zeros(nT,1); bb = zeros(nT,1); cc = zeros(nT,1);
for k=1:nT
  aa(k)=adz_v(k+1)/(adz_v(k)+adz_v(k+1));
  bb(k)=-1.;
  cc(k)=adz_v(k)/(adz_v(k)+adz_v(k+1));
end
% symmetric lower BC
aa(1)=0.;
bb(1)=-(2*adz_v(2)+adz_v(1))/(adz_v(1)+adz_v(2));
% symmetric upper BC
bb(nT)=-(2*adz_v(nT)+adz_v(nT+1))/(adz_v(nT)+adz_v(nT+1));
cc(nT)=0.;
% note that spdiags uses different convention for sub/super-diagonal than diag
sub = [aa(2:nT); 0];        % so -1 diagonal uses aa(2:nT)
sup = [0; cc(1:nT-1)];      % so +1 diagonal uses cc(1:nT-1)
tridiag_L = spdiags([sub bb sup], [-1 0 1], nT, nT);
rhs_diag = spdiags(adz_v(1:nT).*adz_v(2:nT+1)/2, 0, nT, nT);
Linv_rhs = tridiag_L \ rhs_diag;   % avoid inv(...)
Linv_rhs = sparse(Linv_rhs);
m_int = kron(speye(ny), Linv_rhs);
W_u = m_int*m_to_int_u;
W_v = m_int*m_to_int_v;

% 4.2 main parts for matrix
%x->x (memory)
Lx_x = kron(speye(ny), sparse(A));
damp_x = damping*kron(speye(ny), sparse(B*C)*epsdt);
pos = max(0, (abs(lat_deg(:)) - damping_bound_lat) / (max_lat - damping_bound_lat));
epsdt_lat = (pos.^ramp_power) * (dt / damping_bound_time);   % ny x 1
damp_x_bound = kron(spdiags(epsdt_lat, 0, ny, ny), sparse(B*C));   % (ny*120)x(ny*120), sparse
damp_x = damp_x + damp_x_bound;
%d(rhou)dz,d(rhov)dz->x (input)
coef_t=t_wavebg(1:26)*0.;
coef_q=q_wavebg(1:14)*0.;
coef_t(1)=(t_wavebg(2)-t_wavebg(1))/(z(2)-z(1));
coef_t(2:26)=(t_wavebg(3:27)-t_wavebg(1:25))./(z(3:27)-z(1:25));
coef_t=-(coef_t+ggr/Cp)./rho(1:26);
coef_q(1)=(q_wavebg(2)-q_wavebg(1))/(z(2)-z(1));
coef_q(2:14)=(q_wavebg(3:15)-q_wavebg(1:13))./(z(3:15)-z(1:13));
coef_q=-coef_q./rho(1:14);
coef_matrix = spdiags([coef_t; coef_q], 0, 40, 40);
trans_w = [speye(26); [speye(14) sparse(14,12)]];
matrix_1_5 = sparse(B*coef_matrix*trans_w*dt);
Lx_u = kron(speye(ny), matrix_1_5)*W_u;
Lx_v = kron(speye(ny), matrix_1_5)*W_v;
%d(rhov)dz->d(rhou)dz and d(rhou)dz->d(rhov)dz (Coriolis, sparse; BC-update)
% f = beta * y, evaluated on centers and interior faces
%f_lat  = beta * lat(:);        % (ny)
%f_slat = beta * slat_int(:);   % (Nv)
%Fu = kron(spdiags(f_lat,  0, ny, ny), I26) * dt;
%Fv = kron(spdiags(f_slat, 0, Nv, Nv), I26) * dt;
%Lu_v = Fu * slat_to_lat;
%Lv_u = -Fv * lat_to_slat;

% make Coriolis terms skew-symmetric to avoid spurious growth from Coriolis coupling
f_slat = beta * slat_int(:);                 % interior faces
F_big  = kron(spdiags(f_slat,0,Nv,Nv), I26) * dt;

Lu_v = slat_to_lat * F_big;                  % (ny*26) x (Nv*26)
Lv_u = -F_big * lat_to_slat;                 % (Nv*26) x (ny*26)


%x->d(rhou)dz,d(rhov)dz
tmp = spdiags(ggr*rho(1:26)./t_wavebg(1:26), 0, 26, 26) * C(1:26,:);
tmp = sparse(tmp);
Lu_x = -1j*wn_k*kron(speye(ny), tmp)*dt;
Lv_x = -partial_y_lat_to_slat*kron(speye(ny), tmp)*dt;
%d(rhou)dz->d(rhou)dz and %d(rhov)dz->d(rhov)dz
damp_u = damping*speye(ny*26)*epsdt;
damp_v = damping*speye(Nv*26)*epsdt;

% --- high-lat damping for u(lat) and v(slat interior) only (test_bc style, ramp^2) ---
% u lives on lat_deg (ny centers), v lives on slat_int_deg (Nv = ny-1 interior faces)
pos_u = max(0, (abs(lat_deg(:)) - damping_bound_lat) / (max_lat - damping_bound_lat));
damp_bound_u_single_lev = pos_u.^ramp_power / damping_bound_time * dt;    % (ny,1)
pos_v = max(0, (abs(slat_int_deg(:)) - damping_bound_lat) / (max_lat - damping_bound_lat));
damp_bound_v_single_lev = pos_v.^ramp_power / damping_bound_time * dt;    % (Nv,1)

% add to existing damp_u/damp_v (both sparse)
damp_u = damp_u + damping_bound * kron(spdiags(damp_bound_u_single_lev, 0, ny, ny), I26);
damp_v = damp_v + damping_bound * kron(spdiags(damp_bound_v_single_lev, 0, Nv, Nv), I26);

% add diffusion
%partial2_y_x = kron(Lyy_lat_single, I120);   % (ny*120) x (ny*120), sparse
%diff_x = diffusion * nu_visc * dt * ...
%         (-wn_k^2*speye(120*ny) + partial2_y_x);
Dlat = (-wn_k^2*speye(ny) + Lyy_lat_single);
diff_x = diffusion * nu_visc * dt * kron(Dlat, sparse(B*C));
diff_u = diffusion * nu_visc * dt * ...
         (-wn_k^2*speye(26*ny) + partial2_y_lat_to_lat);
diff_v = diffusion * nu_visc * dt * ...
         (-wn_k^2*speye(26*Nv) + partial2_y_slat_to_slat);

% 4.3 build the matrix from components
% option 1: fully explicit (forward)
if opt.stepper == "forward"
    MATRIX = [...
        Lx_x + diff_x - damp_x, Lx_u, Lx_v;
        Lu_x, eye(ny*26) + diff_u - damp_u - hdiff_u, Lu_v;
        Lv_x, Lv_u, eye(Nv*26) + diff_v - damp_v - hdiff_v];
% option 2: make it half implicit time-stepping,
%           forward for state vector, backward for u and v
elseif opt.stepper == "semi-implicit"
    nx = ny*120; nu_dim = ny*26; nv_dim = Nv*26;
    Zxu = sparse(nx, nu_dim);  Zxv = sparse(nx, nv_dim);
    MATRIX_LHS = [ ...
        speye(nx), Zxu, Zxv; ...
        -Lu_x, speye(nu_dim)+damp_u-diff_u+hdiff_u, -Lu_v; ...
        -Lv_x, -Lv_u, speye(nv_dim)+damp_v-diff_v+hdiff_v];
    MATRIX_RHS = [ ...
        Lx_x - damp_x+diff_x, Lx_u, Lx_v; ...
        sparse(nu_dim+nv_dim, nx), speye(nu_dim+nv_dim)];
    MATRIX = MATRIX_LHS \ MATRIX_RHS;
% option 3: --- CN/Cayley-style step mapping (test_bc style), keeping your block names ---
elseif opt.stepper == "CN"
    nx = ny*120; nu_dim = ny*26; nv_dim = Nv*26;

    Zxu = sparse(nx, nu_dim);  Zxv = sparse(nx, nv_dim);
    Zux = sparse(nu_dim, nx);  Zuv = sparse(nu_dim, nv_dim);
    Zvx = sparse(nv_dim, nx);  Zvu = sparse(nv_dim, nu_dim);

    I  = speye(nx + nu_dim + nv_dim);
    Ix = speye(nx);

    % Wave-coupling operator (off-diagonal only)
    R_wave = [ ...
        sparse(nx,nx),  Lx_u,  Lx_v; ...
        Lu_x,           sparse(nu_dim,nu_dim),  Lu_v; ...
        Lv_x,           Lv_u,  sparse(nv_dim,nv_dim)];

    % Dissipation (placed on RHS, like test_bc)
    % NOTE: your LHS-form used +damp -diff +hdiff. Here we add tendencies:
    %   -damp  + diff  - hdiff   (hdiff_u/v are positive semidefinite as you defined them)
    DISS = blkdiag( ...
        -damp_x + diff_x, ...
        (-damp_u + diff_u - hdiff_u), ...
        (-damp_v + diff_v - hdiff_v) );

    % Put the discrete x-step Lx_x into RHS (since it's already a one-step map)
    SSM_step = blkdiag(Lx_x - Ix, sparse(nu_dim+nv_dim, nu_dim+nv_dim));

    LHS = I - 0.5 * R_wave;
    RHS = I + 0.5 * R_wave + DISS + SSM_step;

    MATRIX = LHS \ RHS;
else
    error("Unknown stepper type: " + opt.stepper);

end


function opt = fill_defaults(opt)
opt = def(opt, "dlat", 2);
opt = def(opt, "max_lat", 60);

opt = def(opt, "stepper", "CN");

if ~isfield(opt,"bc"); opt.bc = struct(); end
opt.bc = def(opt.bc, "upper", "rigidlid");
opt.bc = def(opt.bc, "merid", "Neumann");

if ~isfield(opt,"damping"); opt.damping = struct(); end
opt.damping = def(opt.damping, "enable", true);
opt.damping = def(opt.damping, "tau_day", 1);

if ~isfield(opt,"sponge"); opt.sponge = struct(); end
opt.sponge = def(opt.sponge, "enable", true);
opt.sponge = def(opt.sponge, "lat0_deg", 30);
opt.sponge = def(opt.sponge, "tau_hr", 2);
opt.sponge = def(opt.sponge, "ramp_power", 2);

if ~isfield(opt,"diffusion"); opt.diffusion = struct(); end
opt.diffusion = def(opt.diffusion, "enable", true);
opt.diffusion = def(opt.diffusion, "nu_visc", 4e4);

if ~isfield(opt,"hyperdiff"); opt.hyperdiff = struct(); end
opt.hyperdiff = def(opt.hyperdiff, "scale_step", 1.0);

opt = def(opt, "out_root", "linear_wave_ssm_results");
opt = def(opt, "out_case", "test_bc");
end

function s = def(s, field, val)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = val;
end
end