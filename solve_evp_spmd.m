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

    % ---- options (single source of truth) ----
    opt = struct();

    % grid / domain
    opt.dlat    = 2;
    opt.max_lat = 60;

    % time step / stepping choices (mostly for naming/repro)
    opt.stepper = "CN";           % "CN" or "forward" or "semi-implicit"

    % upper boundary condition
    opt.bc.upper = "rigidlid";    % "rigidlid" or "radiation"
    opt.bc.merid = "Neumann";     % "Neumann"

    % damping (uniform)
    opt.damping.enable = true;
    opt.damping.tau_day = 1;     % 1 day -> eps = 1/(tau_day*86400)

    % sponge (high-lat u/v)
    opt.sponge.enable = true;
    opt.sponge.lat0_deg = 30;
    opt.sponge.tau_hr   = 2;     % 2 hours at boundary
    opt.sponge.ramp_power = 2;

    % diffusion
    opt.diffusion.enable = true;
    opt.diffusion.nu_visc = 4e4;          % m^2/s

    % hyperdiffusion
    opt.hyperdiff.scale_step = 1.0;  % <1 means on; 1 means off

    % output
    opt.out_root = "linear_wave_ssm_results";
    opt.out_case = "test_bc";        % subdir

    outdir = fullfile(opt.out_root, opt.out_case);
    if labindex == 1 && ~exist(outdir, 'dir'); mkdir(outdir); end
    labBarrier;
    
    % 1-4 build matrix
    t = tic;
    [MATRIX, wn_k, lat, lat_deg, z, ny, dt, second] = build_matrix(wavenumber_factor, opt);
    t_build = toc(t);
    fprintf('wn=%02d | build_matrix: %.2f s\n', wavenumber_factor, t_build);
    
    % 5 eigs
    disp(['Start solving eigenvalues for ' num2str(wavenumber_factor)])
    t = tic;
    [V, D] = eig(full(MATRIX), 'vector');
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

    tag = sprintf('%s_%s_%s_latbnd%d_dlat%g_damp%gday_spg%d_%ghr_p%d_diff%.1e_hydiff%g', ...
        string(opt.stepper), string(opt.bc.merid), string(opt.bc.upper), ...
        opt.max_lat, opt.dlat, ...
        opt.damping.tau_day, ...
        opt.sponge.lat0_deg, opt.sponge.tau_hr, opt.sponge.ramp_power, ...
        opt.diffusion.nu_visc, ...
        opt.hyperdiff.scale_step);

    filename = fullfile(outdir, sprintf('%s_wn_%02d.mat', tag, wavenumber_factor));
    spmdsave(filename, wavenumber_factor,wn_k,lat,lat_deg,z,MATRIX,D,V,opt);
    t_save = toc(t);
    fprintf('wn=%02d | save: %.2f s\n', wavenumber_factor, t_save);
    disp(['Finished for ' num2str(wavenumber_factor)])
% End of parallel region.............................................
end
delete(gcp);
exit;
