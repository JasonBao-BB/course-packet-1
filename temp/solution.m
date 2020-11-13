% testing_code.m

close all
clear all
clc
% add the directory of helper functions to the path. note that addpath
% checks for duplicates before appending, so this is only done once per
% session.
addpath trc-tools

% add the directory of analysis methods to the path.
addpath analysis-tools

% make a change for plotting
set(groot, 'DefaultTextInterpreter', 'none')
set(groot, 'DefaultLegendInterpreter', 'none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'DefaultFigureColor',[1 1 1])

% load a sample trc file. don't forget the semicolon; it's a lot of data.
% f = "study-data/trc-files/FOG-056/offmed-TUG-cognitive1-TP.trc";

% set this flag to enable file i/o as binary .mat
LEGACY = true;

if ~LEGACY
    % load a sample trc file. don't forget the semicolon; it's a lot of data.
    f = 'offmed-TUG-standard1-TP.trc';
    d = read_trc(f);
else
    % load a pre-loaded .mat file containing the trc data.
    f = 'offmed-TUG-standard1-TP.mat';
    d = read_mat(f);
end



%% Part_2 Truncate 30 seconds of Data

Data = table2array(d);
k = Data(:,2);
FindIndex=find(k==30);
TruncatedData = Data(1: FindIndex, :);

%% Part_3 Plot X-Coordinate of TOP HEAD
X_TopHead = TruncatedData(: , 3);
figure (1)
plot(X_TopHead,'LineWidth',2)
title('X-Coordinate of TOP HEAD');

%% Part_4 Plot Z-Coordinates of TOP HEAD, ASIS and HEEL markers
ZTopHead = TruncatedData(: , 5);
Z_Asis    = TruncatedData(: , 125);
Z_LHeel  = TruncatedData(: , 149);
Z_RHeel  = TruncatedData(: , 176);

figure (2)
plot(ZTopHead,'LineWidth',2)
hold on
plot(Z_Asis,'LineWidth',2)
plot(Z_LHeel,'LineWidth',2)
plot(Z_RHeel,'LineWidth',2)
hold off
title('Z-Coordinates of various markers');
legend('Z-Coordinate of TOP-HEAD' , 'Z-Coordinate of ASIS' , 'Z-Coordinate of L-HEEL' , 'Z-Coordinate of R-HEEL');

%% Part_5 Visualize spectrum of each of HEEL marker

X_RHeel = TruncatedData(: , 147);
Y_RHeel = TruncatedData(: , 148);
Z_RHeel = TruncatedData(: , 149);

X_LHeel = TruncatedData(: , 174);
Y_LHeel = TruncatedData(: , 175);
Z_LHeel = TruncatedData(: , 176);

timeMarker=TruncatedData(: , 2);
dt = mean(diff(timeMarker));
fs = 1 / dt; % sample frequency (Hz)        
n = length(X_RHeel);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range


sX_RHeel = fft(X_RHeel);
sY_RHeel = fft(Y_RHeel);
sZ_RHeel = fft(Z_RHeel);

sX_LHeel = fft(X_LHeel);
sY_LHeel = fft(Y_LHeel);
sZ_LHeel = fft(Z_LHeel);


pX_RHeel = abs(sX_RHeel).^2/n;
pY_RHeel = abs(sY_RHeel).^2/n;
pZ_RHeel = abs(sZ_RHeel).^2/n;

pX_LHeel = abs(sX_LHeel).^2/n;
pY_LHeel = abs(sY_LHeel).^2/n;
pZ_LHeel = abs(sZ_LHeel).^2/n;

figure (3)
subplot(2,3,1)
plot(f,pX_RHeel,'LineWidth',2)
title('X-Coordinate of R-HEEL');
xlabel('Frequency')
ylabel('Power')

subplot(2,3,2)
plot(f,pY_RHeel,'LineWidth',2)
title('Y-Coordinate of R-HEEL');
xlabel('Frequency')
ylabel('Power')

subplot(2,3,3)
plot(f,pZ_RHeel,'LineWidth',2)
title('Z-Coordinate of R-HEEL');
xlabel('Frequency')
ylabel('Power')

subplot(2,3,4)
plot(f,pX_LHeel,'LineWidth',2)
title('X-Coordinate of L-HEEL');
xlabel('Frequency')
ylabel('Power')

subplot(2,3,5)
plot(f,pY_LHeel,'LineWidth',2)
title('Y-Coordinate of L-HEEL');
xlabel('Frequency')
ylabel('Power')

subplot(2,3,6)
plot(f,pZ_LHeel,'LineWidth',2)
title('Z-Coordinate of L-HEEL');
xlabel('Frequency')
ylabel('Power')


%% Part_6 Estimate Time varying power in 5-15 Hz window of each of HEEL marker

psdestx = psd(spectrum.periodogram,X_RHeel,'Fs',fs,'NFFT',length(X_RHeel));
pwr_X_RHeel= avgpower(psdestx,[5 15]);

psdestx = psd(spectrum.periodogram,Y_RHeel,'Fs',fs,'NFFT',length(Y_RHeel));
pwr_Y_RHeel= avgpower(psdestx,[5 15]);

psdestx = psd(spectrum.periodogram,Z_RHeel,'Fs',fs,'NFFT',length(Z_RHeel));
pwr_Z_RHeel= avgpower(psdestx,[5 15]);

psdestx = psd(spectrum.periodogram,X_LHeel,'Fs',fs,'NFFT',length(X_LHeel));
pwr_X_LHeel= avgpower(psdestx,[5 15]);

psdestx = psd(spectrum.periodogram,Y_LHeel,'Fs',fs,'NFFT',length(Y_LHeel));
pwr_Y_LHeel= avgpower(psdestx,[5 15]);

psdestx = psd(spectrum.periodogram,Z_LHeel,'Fs',fs,'NFFT',length(Z_LHeel));
pwr_Z_LHeel= avgpower(psdestx,[5 15]);



% Estimate band power

bp_X_RHeel = bandpower(X_RHeel,fs,[5 15]);
bp_Y_RHeel = bandpower(Y_RHeel,fs,[5 15]);
bp_Z_RHeel = bandpower(Z_RHeel,fs,[5 15]);

bp_X_LHeel = bandpower(X_LHeel,fs,[5 15]);
bp_Y_LHeel = bandpower(Y_LHeel,fs,[5 15]);
bp_Z_LHeel = bandpower(Z_LHeel,fs,[5 15]);

% Estimate NORM
l2norm_X_RHeel = norm(X_RHeel,2)^2/numel(X_RHeel);
l2norm_Y_RHeel = norm(Y_RHeel,2)^2/numel(Y_RHeel);
l2norm_Z_RHeel = norm(Z_RHeel,2)^2/numel(Z_RHeel);

l2norm_X_LHeel = norm(X_LHeel,2)^2/numel(X_LHeel);
l2norm_Y_LHeel = norm(Y_LHeel,2)^2/numel(Y_LHeel);
l2norm_Z_LHeel = norm(Z_LHeel,2)^2/numel(Z_LHeel);

% Perform smoothing
% Apply lowpass filter at 5 Hz frequency

Clean_X_RHeel = lowpass(X_RHeel,5,fs);
Clean_Y_RHeel = lowpass(Y_RHeel,5,fs);
Clean_Z_RHeel = lowpass(Z_RHeel,5,fs);

Clean_X_LHeel = lowpass(X_LHeel,5,fs);
Clean_Y_LHeel = lowpass(Y_LHeel,5,fs);
Clean_Z_LHeel = lowpass(Z_LHeel,5,fs);


% Apply one-dimensional median filter to the input vector
Clean_X_RHeel = medfilt1(X_RHeel);
Clean_Y_RHeel = medfilt1(Y_RHeel);
Clean_Z_RHeel = medfilt1(Z_RHeel);

Clean_X_LHeel = medfilt1(X_LHeel);
Clean_Y_LHeel = medfilt1(Y_LHeel);
Clean_Z_LHeel = medfilt1(Z_LHeel);


% Savitzky-Golay Filter
% Specify a polynomial order of 3 and a frame length of 11.
order = 3;
framelen = 11;


Clean_X_RHeel = sgolayfilt(X_RHeel,order,framelen);
Clean_Y_RHeel = sgolayfilt(Y_RHeel,order,framelen);
Clean_Z_RHeel = sgolayfilt(Z_RHeel,order,framelen);

Clean_X_LHeel = sgolayfilt(X_LHeel,order,framelen);
Clean_Y_LHeel = sgolayfilt(Y_LHeel,order,framelen);
Clean_Z_LHeel = sgolayfilt(Z_LHeel,order,framelen);

% It requires the current directory of "course-packet-main" folder
% Auto setup the location
my_path = mfilename('fullpath');
[filepath, name, ext] = fileparts(my_path);
slocation = filepath;
% Change the location to run the code
% slocation = 'C:\Users\bb\Documents\GitHub\course-packet-1';
% Enter 0 or 1 for the plotting of a graph
ifplot = 1;
[R_HEEL_FOG , L_HEEL_FOG] = calculate_fog_BBY(slocation , ifplot);
