% Shaojun Yu
% shaojun.yu@emory.edu
% given a .trc file f, makes a plot of the heal markers and
% calculate the percentage of freeze
function [R_HEEL_FOG, L_HEEL_FOG] = calculate_fog_SJY(trc_file, make_plot_yesno)

% read the data
if class(trc_file)=="table"
    d = trc_file;
else
    d = read_trc(trc_file);
end

% Truncate to first 30 seconds of file 
d = d(d.Time < 30,:);

% plot
if make_plot_yesno
    subplot(2,1,1)
    plot(d.Time, d.Top_Head_X)
    xlabel("Time, Seconds")
    ylabel("X, mm")
    legend("Top_Head_X")

    subplot(2,1,2)
    hold on
    plot(d.Time, d.Top_Head_Z)
    plot(d.Time, d.L_ASIS_Z)
    plot(d.Time, d.R_ASIS_Z)
    plot(d.Time, d.L_Heel_Z)
    plot(d.Time, d.R_Heel_Z)
    xlabel("Time, Seconds")
    ylabel("Z, mm")
    legend("Top_Head_Z", "L_ASIS_Z", "R_ASIS_Z", "L_Heel_Z", "R_Heel_Z")
end


% sample rate, Hz
Fs = 120;
z = d{:,["L_Heel_Z" "R_Heel_Z"]};
% bandpass filter the Z data between 5 and 15 Hz using the Matlab default
% filter
z = bandpass(z,[5 15],Fs);

NFFT = nrow(z);

% FFT frequencies
F = ((0:1/NFFT:1-1/NFFT)*Fs).';

% FFT values
Z = fft(z,NFFT);

magnitudeZ = abs(Z);        % Magnitude of the FFT
% phaseZ = unwrap(angle(Z));  % Phase of the FFT
% Fmax = 15;                   
if make_plot_yesno
    figure;
    plot(F, 20*log10(magnitudeZ)); %rescale the 
    xlim([5 15])
    xlabel('Frequency in Hz')
    ylabel('dB')
    legend(["L_Heel_Z" "R_Heel_Z"]);
end

zdot = sgolayderiv(d.L_Heel_Z, d.Time);
zddot = abs(sgolayderiv(zdot, d.Time));
L_HEEL_FOG = apply_freeze_thresh(zddot);

zdot = sgolayderiv(d.R_Heel_Z, d.Time);
zddot = abs(sgolayderiv(zdot, d.Time));
R_HEEL_FOG = apply_freeze_thresh(zddot);
end