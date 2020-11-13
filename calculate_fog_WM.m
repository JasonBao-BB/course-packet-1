function [R_HEEL_FOG, L_HEEL_FOG] = calculate_fog_WM(file_location, make_plot_yesno)
    % @author: Wenjing Ma (wenjing.ma@emory.edu)

    %%% Input parameters:
    % file_location: .trc file
    % make_plot_yesno: true or false whether to plot
    %%% Output_parameters:
    % R_HEEL_FOG/L_HEEL_FOG: right/left proportion of FOG

    data = read_trc(file_location);  % read data
    trunc_data = data(data.Time <= 30, :); % truncate data to first 30 s
    Fs = 125; % sample frequency
    
    if make_plot_yesno == true
        figure(1)
        h1 = subplot(2, 1, 1);
        set(h1, 'OuterPosition', [0,0.7, 1, 0.3]);

        % plot x coordinate of Top_Head
        plot(trunc_data.Time, trunc_data.Top_Head_X, 'LineWidth', 1);
        xlabel("Time (s)");
        ylabel("X, mm")
        legend("Top_Head_X", "Location", "northeast");
        
        h2 = subplot(2, 1, 2);
        set(h2, 'OuterPosition', [0,0, 1, 0.7]);

        % plot z coordinate of Top_Head, ASIS and Heel
        plot(trunc_data.Time, trunc_data.Top_Head_Z, "b", 'LineWidth', 1);
        hold on
        plot(trunc_data.Time, trunc_data.R_ASIS_Z, "r", 'LineWidth', 1);
        plot(trunc_data.Time, trunc_data.L_ASIS_Z, "y", 'LineWidth', 1);
        plot(trunc_data.Time, trunc_data.R_Heel_Z, "Color", [0.5 0 0.5], 'LineWidth', 1);
        plot(trunc_data.Time, trunc_data.L_Heel_Z, "g", 'LineWidth', 1);
        hold off
        xlabel("Time (s)");
        ylabel("Z, mm")
        legend("Top_Head_Z", "R_ASIS_Z", "L_ASIS_Z", "R_Heel_Z", "L_Heel_Z");
    end
    
    % estimate time-varying power in 5-15Hz window
    % a. bandpower function (refer to analysis-tools/bandpowerwrapper.m)
    frameLen = 250;  % 250 dt sliding window ~ 2s sliding window
    for i = 1:(length(trunc_data.L_Heel_Z)-(frameLen-1))
        inds = i:(i+frameLen-1);
        xdata = trunc_data.L_Heel_Z(inds);
        ptot_L(i) = bandpower(xdata, Fs,[0 Fs/2]);  % total power
        pband_L(i) = bandpower(xdata, Fs,[5 15]);   % [5 15] window
        prop_L(i) = pband_L(i)/ptot_L(i);           % percentage
    end
    
    for i = 1:(length(trunc_data.R_Heel_Z)-(frameLen-1))
        inds = i:(i+frameLen-1);
        xdata = trunc_data.R_Heel_Z(inds);
        ptot_R(i) = bandpower(xdata, Fs,[0 Fs/2]);
        pband_R(i) = bandpower(xdata, Fs,[5 15]);
        prop_R(i) = pband_R(i)/ptot_R(i);
    end

%     power_R_heel = bandpower(trunc_data.R_Heel_Z, Fs, [5 15]);
%     power_L_heel = bandpower(trunc_data.L_Heel_Z, Fs, [5 15]);
    
    % b. trc-tools/sgolayderiv for smooth
    % smooth with Savitsky-Golay
    if ~mod(frameLen,2)
        frameLen = frameLen+1;  % need odd framelen
    end
    
    order = 5;
    % edge effects 
    [~,g] = sgolay(5,frameLen);  % frameLen is 250+1, order 5
    dt = 1/Fs;

    % convolve, p=1, first derivative, get smoothed output
    p = 0;
    % calculate the percentage of bandpower in total
    smoothed_R = conv(prop_R, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    smoothed_L = conv(prop_L, factorial(p)/(-dt)^p * g(:,p+1), 'same');

    % plot spectrum along with smoothed power
    if make_plot_yesno == true
        figure(2)
        subplot(2, 1, 1);
        % visualize spectrum of each of heal markers
        % Right heel
        pspectrum(trunc_data.R_Heel_Z, Fs, 'spectrogram', 'OverlapPercent', 99, ...
        'MinThreshold', -10, 'FrequencyResolution', 1, 'Reassign', true);
        ylim([0 20]);
        title("R_Heel_Z");
        hold on
        plot(trunc_data.Time(1: length(smoothed_R)), smoothed_R*2000, "Color", [0.5 0 0.5], 'LineWidth', 2);
        legend("R_Heel_Z Freeze Band Power (nu)")
        hold off


        % Left heel power spectrum
        subplot(2, 1, 2)
        pspectrum(trunc_data.L_Heel_Z, Fs, 'spectrogram', 'OverlapPercent', 99, ...
        'MinThreshold', -10, 'FrequencyResolution', 1, 'Reassign', true);
        ylim([0 20]);
        title("L_Heel_Z")
        hold on
        plot(trunc_data.Time(1: length(smoothed_L)), smoothed_L*2000, "g", 'LineWidth', 2);
        legend("L_Heel_Z Freeze Band Power (nu)")
        hold off
        
        % Note: smoothed_L/smoothed_R * 2000 is to adjust to the y axis
        % which is at most 20. But I actually don't understand why do we
        % plot the percentage and the power spectrum in one plot. Is it
        % just for visualization?
    end
    
    % calculate frozen percentage
    frozen_R = all(smoothed_R >= 0.001, 1);
    frozen_L = all(smoothed_L >= 0.001, 1);
    
    R_HEEL_FOG = round(100*sum(frozen_R)/length(frozen_R),1);
    L_HEEL_FOG = round(100*sum(frozen_L)/length(frozen_L),1); 
end