%% build matrix for the EVP with given parameters
%current parameters: wavenumber_factor, dlat

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
mass_weight=diag(sqrt([dm(1:26);2.5^2*dm(1:14)]));

% 1.4 meridional grids
max_lat = 60; % degrees
lat_deg=-max_lat:dlat:max_lat;
% staggered grids
slat_deg=lat_deg-(lat_deg(2)-lat_deg(1))/2;
slat_deg(end+1)=lat_deg(end)+(lat_deg(2)-lat_deg(1))/2;
ny=size(lat_deg,2);
lat=lat_deg*pi/180*RE; % meters
slat=slat_deg*pi/180*RE;
dy=lat(2)-lat(1);

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
damping=true;
eps=0./(86400*second);
epsdt=eps*dt;

% 3.2 damp U, V at higher latitudes
damping_bound=true;
damping_bound_lat=40.;
damping_bound_time=14400.*second;

% 3.3 add diffusion to U, V in zonal and meridional direction
diffusion=true;
nu = 0e5 * meter^2 / second;

%% 4. build the matrix
% vector        = [ x,   d(rhou)dz, d(rhov)dz(staggered)];
% vector length = ny*120 + ny*26 + (ny+1)*26
%          | x->x,         d(rhou)dz->x,         d(rhov)dz->x         |
% matrix = | x->d(rhou)dz, d(rhou)dz->d(rhou)dz, d(rhou)dz->d(rhov)dz |
%          | x->d(rhov)dz, d(rhou)dz->d(rhov)dz, d(rhov)dz->d(rhov)dz |
% assuming zero signal in higher latitudes than max_lat
disp(['Start building matrix for ' num2str(wavenumber_factor)])

% 4.1 auxiliary matrices
%partial2_y_lat_to_lat
partial2_y_lat_to_lat_single_lev=zeros(ny, ny);
for i=2:ny-1
    partial2_y_lat_to_lat_single_lev(i,i-1:i+1)=[1, -2, 1]/(dy*dy);
end
partial2_y_lat_to_lat_single_lev(1,1)=-2/(dy*dy);
partial2_y_lat_to_lat_single_lev(1,2)=1/(dy*dy);
partial2_y_lat_to_lat_single_lev(ny,ny-1)=1/(dy*dy);
partial2_y_lat_to_lat_single_lev(ny,ny)=-2/(dy*dy);
partial2_y_lat_to_lat=kron(partial2_y_lat_to_lat_single_lev,eye(26));

%partial2_y_slat_to_slat
partial2_y_slat_to_slat_single_lev=zeros(ny+1, ny+1);
for i=2:ny
    partial2_y_slat_to_slat_single_lev(i,i-1:i+1)=[1, -2, 1]/(dy*dy);
end
partial2_y_slat_to_slat_single_lev(1,1)=-2/(dy*dy);
partial2_y_slat_to_slat_single_lev(1,2)=1/(dy*dy);
partial2_y_slat_to_slat_single_lev(ny+1,ny)=1/(dy*dy);
partial2_y_slat_to_slat_single_lev(ny+1,ny+1)=-2/(dy*dy);
partial2_y_slat_to_slat=kron(partial2_y_slat_to_slat_single_lev,eye(26));

%partial_y_slat_to_lat
partial_y_slat_to_lat_single_lev=zeros(ny, ny+1);
for i=1:ny
    partial_y_slat_to_lat_single_lev(i,i:i+1)=[-1, 1]/dy;
end
partial_y_slat_to_lat=kron(partial_y_slat_to_lat_single_lev,eye(26));

%partial_y_lat_to_slat
partial_y_lat_to_slat_single_lev=zeros(ny+1, ny);
for i=2:ny
    partial_y_lat_to_slat_single_lev(i,i-1:i)=[-1, 1]/dy;
end
partial_y_lat_to_slat_single_lev(1,1)=1/dy;
partial_y_lat_to_slat_single_lev(ny+1,ny)=-1/dy;
partial_y_lat_to_slat=kron(partial_y_lat_to_slat_single_lev,eye(26));

%compute w from dudz, dvdz
m_to_int_u=-1j*wn_k*eye(26*ny);
m_to_int_v=-1*partial_y_slat_to_lat;
%build vertical laplacian
for k=2:nT
  adz(k)=z(k)-z(k-1);
end
adz(1)=2*z(1);
adz(nT+1)=(z(nT+1)-z(nT));
adz=adz';
for k=1:nT
  aa(k)=adz(k+1)/(adz(k)+adz(k+1));
  bb(k)=-1.;
  cc(k)=adz(k)/(adz(k)+adz(k+1));
end
%symmetric lower BC
aa(1)=0.;
bb(1)=-(2*adz(2)+adz(1))/(adz(1)+adz(2));
%symmetric upper BC
bb(nT)=-(2*adz(nT)+adz(nT+1))/(adz(nT)+adz(nT+1));
cc(nT)=0.;
invL=inv(diag(aa(2:nT),-1)+diag(bb,0)+diag(cc(1:nT-1),1));
m_int=kron(eye(ny), invL*diag(adz(1:nT).*adz(2:nT+1)/2));
W_u = m_int*m_to_int_u;
W_v = m_int*m_to_int_v;

% 4.2 main parts for matrix
%x->x (memory)
Lx_x = kron(eye(ny),A);
damp_x = damping*kron(eye(ny),B*C*epsdt);
%d(rhou)dz,d(rhov)dz->x (input)
coef_t=t_wavebg(1:26)*0.;
coef_q=q_wavebg(1:14)*0.;
coef_t(1)=(t_wavebg(2)-t_wavebg(1))/(z(2)-z(1));
coef_t(2:26)=(t_wavebg(3:27)-t_wavebg(1:25))./(z(3:27)-z(1:25));
coef_t=-(coef_t+ggr/Cp)./rho(1:26);
coef_q(1)=(q_wavebg(2)-q_wavebg(1))/(z(2)-z(1));
coef_q(2:14)=(q_wavebg(3:15)-q_wavebg(1:13))./(z(3:15)-z(1:13));
coef_q=-coef_q./rho(1:14);
coef_matrix=diag([coef_t; coef_q]);
trans_w=eye(40, 26);
trans_w(27:41:14*40)=1;
matrix_1_5=B*coef_matrix*trans_w*dt;
Lx_u = kron(eye(ny),matrix_1_5)*W_u;
Lx_v = kron(eye(ny),matrix_1_5)*W_v;
%d(rhov)dz->d(rhou)dz
beta_v_on_u_single_lev=zeros(ny, ny+1);
for i=1:ny
    beta_v_on_u_single_lev(i, i:i+1)=[0.5,0.5]*beta*lat(i)*dt;
end
Lu_v = kron(beta_v_on_u_single_lev, eye(26));
%d(rhou)dz->d(rhov)dz
beta_u_on_v_single_lev=zeros(ny+1, ny);
for i=2:ny
    beta_u_on_v_single_lev(i, i-1:i)=[-0.5,-0.5]*beta*slat(i)*dt;
end
beta_u_on_v_single_lev(1, 1)=-0.5*beta*slat(1)*dt;
beta_u_on_v_single_lev(ny+1, ny)=-0.5*beta*slat(ny+1)*dt;
Lv_u = kron(beta_u_on_v_single_lev, eye(26));
%x->d(rhou)dz,d(rhov)dz
tmp=diag(ggr*rho(1:26)./t_wavebg(1:26))*C(1:26,:);
Lu_x = -1j*wn_k*kron(eye(ny), tmp)*dt;
Lv_x = -1.*partial_y_lat_to_slat*kron(eye(ny), tmp)*dt;
%d(rhou)dz->d(rhou)dz
damp_u = damping*eye(ny*26)*epsdt;
damp_bound_u_single_lev = zeros(ny, 1);
for i=1:ny
    if abs(lat_deg(i))>=damping_bound_lat
        pos_coef = (abs(lat_deg(i))-damping_bound_lat) / ...
                    (max_lat - damping_bound_lat);
        damp_bound_u_single_lev(i) = pos_coef/damping_bound_time*dt;
    end
end
damp_u = damp_u + ...
         damping_bound * kron(diag(damp_bound_u_single_lev), eye(26));
%d(rhov)dz->d(rhov)dz
damp_v = damping*eye(ny*26+26)*epsdt;
damp_bound_v_single_lev = zeros(ny+1, 1);
for i=1:ny+1
    if abs(slat_deg(i))>=damping_bound_lat
        pos_coef = (abs(slat_deg(i))-damping_bound_lat) / ...
                   (max_lat - damping_bound_lat);
        damp_bound_v_single_lev(i) = pos_coef/damping_bound_time*dt;
    end
end
damp_v = damp_v + ...
         damping_bound * kron(diag(damp_bound_v_single_lev), eye(26));
% add diffusion
diff_u = diffusion * nu * dt * ...
         (-wn_k^2*eye(26*ny) + partial2_y_lat_to_lat);
diff_v = diffusion * nu * dt * ...
         (-wn_k^2*eye(26*ny+26) + partial2_y_slat_to_slat);

% 4.3 build the matrix from components
% option 1: fully explicit (forward)
% MATRIX = [Lx_x - damp_x, Lx_u, Lx_v;
%     Lu_x, eye(ny*26) - damp_u + diff_u, Lu_v;
%     Lv_x, Lv_u, eye(ny*26+26) - damp_v + diff_v];
% option 2: make it half implicit time-stepping,
%           forward for state vector, backward for u and v
MATRIX_LHS = sparse([eye(ny*120, ny*172+26);
    -Lu_x, eye(ny*26)+damp_u-diff_u, -Lu_v;
    -Lv_x, -Lv_u, eye(ny*26+26)+damp_v-diff_v]);
MATRIX_RHS = sparse([Lx_x - damp_x, Lx_u, Lx_v;
    zeros(ny*52+26, ny*120), eye(ny*52+26)]);
MATRIX = MATRIX_LHS \ MATRIX_RHS;
% option 3: ???
