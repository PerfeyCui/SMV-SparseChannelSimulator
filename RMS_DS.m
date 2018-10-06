function [rms_rec_sig,Delay]= RMS_DS(rec_sig,timeResolution)
% 计算均方根时延扩展，利用输入的时间分辨率 timeResolution和信号幅度值等间隔采样序列，rec_sig
Seg = [1:length(rec_sig)];
Delay = timeResolution*Seg*rec_sig.^2./(sum(rec_sig.^2));
rms_rec_sig = sqrt(timeResolution^2*Seg.^2*rec_sig.^2./(sum(rec_sig.^2))-Delay^2);
