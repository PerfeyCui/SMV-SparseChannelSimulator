function [rms_rec_sig,Delay]= RMS_DS(rec_sig,timeResolution)
% ���������ʱ����չ�����������ʱ��ֱ��� timeResolution���źŷ���ֵ�ȼ���������У�rec_sig
Seg = [1:length(rec_sig)];
Delay = timeResolution*Seg*rec_sig.^2./(sum(rec_sig.^2));
rms_rec_sig = sqrt(timeResolution^2*Seg.^2*rec_sig.^2./(sum(rec_sig.^2))-Delay^2);
