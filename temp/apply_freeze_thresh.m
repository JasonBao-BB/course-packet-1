function percent_frozen = apply_freeze_thresh(zddot, ifplot)
% function percent_frozen = apply_freeze_thresh(zddot)
%
% applies a median filter and a threshold operation to heel marker
% acceleration magnitude. returns percent frozen.

filter_order = 500;
magnitude_thresh = 100;
zddot_med = medfilt1(zddot,filter_order,[],1);

frozen = all(zddot_med<magnitude_thresh,2);
percent_frozen = round(100*sum(frozen)/length(frozen),1);

if ifplot==1
    close all
    subplot(3,1,1)
    plot(zddot)
    subplot(3,1,2)
    plot(zddot_med)
    subplot(3,1,3)
    plot(frozen)
    ylim([-0.25 1.25])
end

end