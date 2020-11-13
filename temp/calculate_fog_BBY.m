function [R_HEEL_FOG , L_HEEL_FOG] = calculate_fog_BBY(slocation , ifplot)

% Input Arguments:
%     slocation = Enter the current directory of "course-packet-main"
%     ifplot    = Enter 0 or 1 for the plotting of a graph

% Output Arguments:
%     R_HEEL_FOG = A scalar number 
%     L_HEEL_FOG = A scalar number 



cd(slocation)
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


marker_names = names(d);
markers = marker_names(contains(marker_names,["Heel"]));
z = d{:,markers};
t = d.Time;



RHeel_zdot = sgolayderiv(z(:,1:3),t);
LHeel_zdot = sgolayderiv(z(:,4:6),t);


RHeel_zddot = abs(sgolayderiv(RHeel_zdot,t));
LHeel_zddot = abs(sgolayderiv(LHeel_zdot,t));


R_HEEL_FOG = apply_freeze_thresh(RHeel_zddot, ifplot);
L_HEEL_FOG = apply_freeze_thresh(LHeel_zddot, ifplot);

end