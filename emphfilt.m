function F = emphfilt(cf, bw, A, Ts)
% EMPHFILT Emphasis filter.
% 
% Emphasis filter to be used as additional open loop amplification within
% a give frequency band in the feedback controller.
%
% F = emphfilt(cf, bw, A, Ts)
%
% INPUTS:
%   cf:   Center frequency [Hz].
%   bw:   Bandwidth (-3 dB cutoff) [Hz].
%   A:    Target peak amplification at center frequency 'cf'
%
% OUTPUTS:
%   F:    Filter transfer function (tf object). If 'Ts' is 0 or unspecified 
%         F is a continuous-time transfer function, otherwise it is a
%         discrete-time transfer function with sampling time 'Ts'.

fl2 = 2*pi*cf/sqrt(bw);
fh1 = fl2*bw;
fh2 = fh1*A;
fl1 = fl2/A;
if nargin < 4 || isempty(Ts) || Ts == 0
    F1 = tf([1 fl1], [1 fl2]);
    F2 = tf([1 fh2], [1 fh1]);
else
    F1 = tf([1 -exp(-fl1*Ts)], [1 -exp(-fl2*Ts)], Ts);
    F2 = tf([1 -exp(-fh2*Ts)], [1 -exp(-fh1*Ts)], Ts);
end
F = F1*F2;
