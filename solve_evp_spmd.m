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
    dlat = 2;
    
    % 1-4
    build_matrix;
    
    % 5. compute eigenvalues
    disp(['Start solving eigenvales for ' num2str(wavenumber_factor)])
    % Compute the top k eigenvalues and eigenvectors of A using eigs()
    [V, D]=eigs(MATRIX, 172*ny+26); % currently compute all eigenvalues
    %if working with dense matrix, use eig instead of eigs
    %[V, D]=eig(MATRIX, 'nobalance');
    D = diag(D);

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
    %growth rate (per day)
    D(:,2)=86400*second/dt*log(abs(D(:,1)));
    %frequency (cycle per day)
    %positive means westward propagation
    D(:,3)=86400*second/dt*angle(D(:,1))/(2*pi);
    disp(['Start saving for ' num2str(wavenumber_factor)])
    filename=['linear_wave_ssm_results/test_implicit_bound40_14400_' ...
        num2str(wavenumber_factor, '%02d') '_latbnd_60_dlat_' ...
        num2str(dlat) '_rigidlid'];
    spmdsave(filename, wavenumber_factor,wn_k,lat,lat_deg,z,MATRIX,D,V);
    disp(['Finished for ' num2str(wavenumber_factor)])
% End of parallel region.............................................
end
delete(gcp);
exit;
