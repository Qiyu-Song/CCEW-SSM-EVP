%% plot the eigenstate of a certain eigenmode
wavenumber_factor = 20;
dlat = 2;
%load the eigenvalues if needed
%load(['linear_wave_ssm_results/test_implicit_bound40_14400_' ...
%    num2str(wavenumber_factor, '%02d') '_latbnd_60_dlat_' ...
%    num2str(dlat) '_rigidlid.mat'])
%rebuild the matrix and key auxiliary components if needed
%build_matrix;

%% load or compute data for plotting
ind=1; % choose the eigenmode to be plotted
D(ind,:)
state=V(:,ind);
ny=numel(lat);
dUdz_eig=reshape(state(120*ny+1:146*ny), 26, ny);
dVdz_eig=reshape(state(146*ny+1:172*ny+26), 26, ny+1);
O_eig=reshape(kron(eye(ny),C)*state(1:120*ny), 40, ny);
T_eig=O_eig(1:26,:);
Q_eig=O_eig(27:40,:);
W_eig=reshape(W_uv*state(120*ny+1:172*ny+26), 26, ny);

%% plot crosssection for T and w
figure("Name","Eigenstate_crosssection_TW", 'Position', [0 0 1200 600])
% plot meridional crosssection for T and w, real part and imaginary part
subplot('Position',[0.1,0.56,0.25,0.37])
contourf(lat_deg, z(1:26)/1000, real(T_eig));colorbar;
title("T Real");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");
subplot('Position',[0.1,0.08,0.25,0.37])
contourf(lat_deg, z(1:26)/1000, real(W_eig));colorbar;
title("W Real");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");

subplot('Position',[0.40,0.56,0.25,0.37])
contourf(lat_deg, z(1:26)/1000, imag(T_eig));colorbar;
title("T Imag");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");
subplot('Position',[0.40,0.08,0.25,0.37])
contourf(lat_deg, z(1:26)/1000, imag(W_eig));colorbar;
title("W Imag");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");

% also plot a zonal crosssection for one wavelength along the equator
xmax=2*pi/wn_k;
x=0:xmax/100:xmax;
[~,latind]=min(abs(lat_deg));
T_xz=real(T_eig(:,latind)*exp(1j*wn_k*x));
W_xz=real(W_eig(:,latind)*exp(1j*wn_k*x));
subplot('Position',[0.70,0.56,0.25,0.37])
contourf(x/1000, z(1:26)/1000, T_xz);colorbar;
title(['T at ', num2str(lat_deg(latind)), 'N']);
xlabel("Zonal distance (km)");
ylabel("Altitude (km)");
subplot('Position',[0.70,0.08,0.25,0.37])
contourf(x/1000, z(1:26)/1000, W_xz);colorbar;
title(['W at ', num2str(lat_deg(latind)), 'N']);
xlabel("Zonal distance (km)");
ylabel("Altitude (km)");

%% plot crosssection for U and V (without vertical integration)
figure("Name","Eigenstate_crosssection_UV", 'Position', [0 0 1000 500])
% plot meridional crosssection for U and V, real part and imaginary part
subplot(2,2,1)
contourf(lat_deg, z(1:26)/1000, real(dUdz_eig));colorbar;
title("dUdz Real");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");
subplot(2,2,3)
contourf(slat_deg, z(1:26)/1000, real(dVdz_eig));colorbar;
title("dVdz Real");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");

subplot(2,2,2)
contourf(lat_deg, z(1:26)/1000, imag(dUdz_eig));colorbar;
title("dUdz Imag");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");
subplot(2,2,4)
contourf(slat_deg, z(1:26)/1000, imag(dVdz_eig));colorbar;
title("dVdz Imag");
xlabel("Latitude (degN)");
ylabel("Altitude (km)");