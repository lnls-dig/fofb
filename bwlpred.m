function p = bwlpred(dly, Amax, Ts)
% BWLPRED Band-limited discrete-time predictor.
%
% p = BWLPRED(dly, Amax, Ts)
%
% INPUTS:
%   dly:    Delay (in samples unit) to be predicted (multiple of 0.5).
%   Amax:   (Optional) Maximum amplification (or gain at 1/Ts/2). Default
%           value is 10^dly.
%   Ts:     (Optional) Sample time. If not specified -1 is assumed.
%
% OUTPUTS:
%   p:   Transfer function object of the band-limited predictor.

if nargin < 2 || isempty(Amax)
    Amax = 10^dly;
end
if Amax <= 1
    error('''Amax'' must be greater than 1 to result in a stable predictor.');
end

N = 2*dly;
if rem(N,1) ~= 0 || N < 0
    error('Delay ''dly'' must be positive and multiple of 0.5');
end

if nargin < 3 || isempty(Ts)
    Ts = -1;
end

b = Amax^(1/N);
c = (b+1)/(b-1);
p = tf([1+c 0], [c 1], Ts)^N;
