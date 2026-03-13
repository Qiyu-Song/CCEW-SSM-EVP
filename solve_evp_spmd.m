%%
%-------------------------------------------------------------------------
%Find eigenvalues and eigenstates of state-space model coupled linear wave
%-------------------------------------------------------------------------

parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
% Start of parallel region...........................................
spmd
    nproc = numlabs;  % get total number of workers
    iproc = labindex; % get lab ID
    wavenumber_factor = iproc-1;
    dlat = 1;
    
    % 1-4 build matrix
    t = tic;
    [MATRIX, wn_k, lat, lat_deg, z, ny, dt, second] = build_matrix(wavenumber_factor, dlat);
    t_build = toc(t);
    fprintf('wn=%02d | build_matrix: %.2f s\n', wavenumber_factor, t_build);
    
    % 5 eigs
    disp(['Start solving eigenvalues for ' num2str(wavenumber_factor)])
    t = tic;
    [V, D] = eigs(MATRIX, 172*ny+26);
    D = diag(D);
    t_eigs = toc(t);
    fprintf('wn=%02d | eigs: %.2f s\n', wavenumber_factor, t_eigs);

    % 6. postprocessing
    % 6.1 plot eigenvalues
    %figure;
    %plot(real(D), imag(D), 'bo','MarkerSize',6);
    %hold on;
    %theta = linspace(0, 2*pi, 200);
    %plot(cos(theta), sin(theta), 'k--','LineWidth',1.5);  % unit circle
    %xlabel('Real Part');
    %ylabel('Imaginary Part');
    %title(['Eigenvalues for ' num2str(wavenumber_factor)]);
    %axis equal;
    %grid on;
    %hold off;

    % 6.2 save results
    t = tic;
    %growth rate (per day)
    D(:,2)=86400*second/dt*log(abs(D(:,1)));
    %frequency (cycle per day)
    %positive means westward propagation
    D(:,3)=86400*second/dt*angle(D(:,1))/(2*pi);
    disp(['Start saving for ' num2str(wavenumber_factor)])
    filename=['linear_wave_ssm_results/semi_implicit/bound40_6hr_damp_0_diff_0_wn_' ...
        num2str(wavenumber_factor, '%02d') '_latbnd_60_dlat_' ...
        num2str(dlat) '_rigidlid'];
    spmdsave(filename, wavenumber_factor,wn_k,lat,lat_deg,z,MATRIX,D,V);
    t_save = toc(t);
    fprintf('wn=%02d | save: %.2f s\n', wavenumber_factor, t_save);
    disp(['Finished for ' num2str(wavenumber_factor)])
% End of parallel region.............................................
end
delete(gcp);
exit;
