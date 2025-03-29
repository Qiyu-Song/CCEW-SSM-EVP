figure;
for wavenumber_factor=0:40
    wavenumber_factor
    filepath='linear_wave_ssm_results/';
    filename=['test_implicit_bound40_14400_' num2str(wavenumber_factor, '%02d') '_latbnd_60_dlat_2_rigidlid.mat'];
    load([filepath filename],"D");
    D = D(:,1);
    indices = ((abs(D) > exp(0.3/96)) & (abs(angle(D)*96/(2*pi)) > 0.03*wavenumber_factor)) | ((abs(D) <= exp(0.3/96)) & (abs(D) > 1) & (abs(angle(D)*96/(2*pi)) < 0.4));
    D = D(indices);
    Dp = D(angle(D)*96/(2*pi)>0.3);
    Dn = D(angle(D)<0);

    [~, idxSort] = sort(imag(Dn), 'descend');
    kelvin_d = Dn(idxSort(1));
    MRG_d = Dn(idxSort(2));
    [~, idxSort] = sort(imag(Dp), 'ascend');
    WIG1_d = Dp(idxSort(1));
    
    % frequency (cycle per day)
    subplot(1,3,1)
    scatter(-wavenumber_factor*ones(size(D)).*sign(angle(D)), abs(angle(D))*96/(2*pi), max(500*log(abs(D))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'k');
    scatter(wavenumber_factor, -angle(kelvin_d)*96/(2*pi), max(500*log(abs(kelvin_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'r');
    scatter(wavenumber_factor, -angle(MRG_d)*96/(2*pi), max(500*log(abs(MRG_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'g');
    scatter(-wavenumber_factor, angle(WIG1_d)*96/(2*pi), max(500*log(abs(WIG1_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'b');
    hold on;
    xlim([-41, 41]);
    ylim([0, 1])
    title('Frequency (cpd)');
    xlabel('Zonal Wavenumber')
    grid on
    box on
    
    % growth rate (per day)
    subplot(1,3,2)
    scatter(wavenumber_factor*ones(size(Dn)), log(abs(Dn))*96, max(500*log(abs(Dn))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'k');
    scatter(-wavenumber_factor*ones(size(Dp)), log(abs(Dp))*96, max(500*log(abs(Dp))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'k');
    scatter(wavenumber_factor, log(abs(kelvin_d))*96, max(500*log(abs(kelvin_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'r');
    scatter(wavenumber_factor, log(abs(MRG_d))*96, max(500*log(abs(MRG_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'g');
    scatter(-wavenumber_factor, log(abs(WIG1_d))*96, max(500*log(abs(WIG1_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'b');
    hold on;
    xlim([-41, 41]);
    ylim([0.9, 1.4]);
    title('Growth Rate (1/day)');
    xlabel('Zonal Wavenumber')
    grid on
    box on
    
    % phase speed (m/s)
    subplot(1,3,3)
    scatter(wavenumber_factor*ones(size(Dn)), -angle(Dn)/900 / (wavenumber_factor/6370e3), max(500*log(abs(Dn))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'k');
    scatter(-wavenumber_factor*ones(size(Dp)), angle(Dp)/900 / (wavenumber_factor/6370e3), max(500*log(abs(Dp))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'k');
    scatter(wavenumber_factor, abs(angle(kelvin_d))/900 / (wavenumber_factor/6370e3), max(500*log(abs(kelvin_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'r');
    scatter(wavenumber_factor, abs(angle(MRG_d))/900 / (wavenumber_factor/6370e3), max(500*log(abs(MRG_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'g');
    scatter(-wavenumber_factor, abs(angle(WIG1_d))/900 / (wavenumber_factor/6370e3), max(500*log(abs(WIG1_d))*96, 0), 'Marker', '.', 'MarkerEdgeColor', 'b');
    hold on;
    xlim([-41, 41]);
    ylim([13, 25]);
    title('Phase Speed (m/s)');
    xlabel('Zonal Wavenumber')
    grid on
    box on
end