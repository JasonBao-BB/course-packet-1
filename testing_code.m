% testing_code.m

close all
clear all

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

% look at a small section
peek(d)

% if the LEGACY flag has been set, send the dataset to the calculate_fog_JLM method.
% otherwise send the name of the .trc file.

[R_HEEL_FOG, L_HEEL_FOG] = calculate_fog_BBY(f, true);
