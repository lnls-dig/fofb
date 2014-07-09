function [argout1, fcmode] = fcsweep(varargin)

% Parameters
amplitude = 0.1; % RMS value (randn) or sinusoid amplitude
expduration = 100; % each experiment's duration (packet count)
profile_mode = 'corr';
excitation_mode = 'randn';
freqs = [1:10 20:10:100 200:100:500];
Ts = 320e-6;
% ---

ncols = 50;

if ischar(varargin{1}) && strcmpi(varargin, 'size')
    argout1 = ncols;
else
    i = varargin{1};
    npts_packet = varargin{2};
    
    fprintf('iteration #%d: ', i);
    
    expnumber = floor(i/(expduration+1));

    if rem(i,(expduration+1)) == 0
        argout1 = zeros(npts_packet, ncols);
        fcmode = uint32(0);
        fprintf('pause');
    else
        repeatprofiles_after = select_profile('size', profile_mode);
        profile_number = rem(expnumber, repeatprofiles_after);
        [profile, fcmode] = select_profile(profile_number, profile_mode);
        fprintf('profile #%d', profile_number);
        if strcmpi(excitation_mode, 'randn')
            packet = amplitude*randn(npts_packet, size(profile, 1))*profile;
            fprintf(' - std = %f (Gaussian white noise)', amplitude);
        elseif strcmpi(excitation_mode, 'sin_sweep')
            f = freqs(rem(floor(expnumber/repeatprofiles_after), length(freqs))+1);
            t = (((i-1)*npts_packet+(0:npts_packet-1))*Ts)';
            packet = amplitude*repmat(cos(2*pi*f*t), 1, size(profile, 1))*profile;
            fprintf(' - amplitude = %f, frequency = %f Hz (sinusoid)', amplitude, f);
        end
        % Zero-padding
        argout1 = [packet zeros(npts_packet, ncols-size(packet,2))];
    end
    fprintf('\n');
end

function [argout1, fcmode] = select_profile(argin1, mode)

if ischar(argin1) && strcmpi(argin1, 'size')
    if strcmpi(mode, 'corr')
        argout1 = 43;
    elseif  strcmpi(mode, 'corr_profiles') || strcmpi(mode, 'orb_profiles')
        profiles = evalin('base', 'profiles');
        argout1 = size(profiles,2);
    end
else
    i = argin1;
    if strcmpi(mode, 'corr')
        if i == 0
            argout1 = diag(ones(1,42));
        else
            argout1 = [zeros(1,(i-1)) 1 zeros(1,42-i)];
        end
    elseif  strcmpi(mode, 'corr_profiles') || strcmpi(mode, 'orb_profiles')
        profiles = evalin('base', 'profiles');
        argout1 = profiles(:,i+1)';
    end

    if strcmpi(mode, 'corr') || strcmpi(mode, 'corr_profiles')
        fcmode = uint32(2);
    elseif strcmpi(mode, 'orb_profiles')
        fcmode = bitshift(uint32(2), 16);
    else
        fcmode = uint32(0);
    end
end
