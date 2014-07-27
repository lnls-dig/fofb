function [fcdata, expout] = fcsend(ip, fcn, expinfo, npts_packet, stopat)

if ischar(ip) || ~isfield(ip,'port')
    if ischar(ip)
        ip_.address = ip;
    end
    ip_.port = 3604;
else
    ip_ = ip;
end

if nargin < 4 || isempty(npts_packet)
    npts_packet = 1000;
end

if nargin < 5 || isempty(stopat)
    stopat = expinfo.duration;
end

[fcdata, expout] = fcn('init', npts_packet, expinfo);

i=0;
failure = false;

conn = tcpip(ip_.address, ip_.port, 'OutputBufferSize', 10*npts_packet*(expinfo.ncols+1)*4);
fopen(conn);
while i < stopat
    [packet, expinterval] = fcn(i, npts_packet, expinfo, fcdata);
    
    %data((1+i*npts_packet):(i+1)*npts_packet, :) = packet;
    
    % Ensure data is encoded in 4-byte floating point representation
    packet = single(packet);
    
    if expinterval
        fcmode = uint32(0);
    else
        fcmode = fcbuildmode(expinfo.mode);
    end
    
    while true
        try
            if fread(conn, 1, 'uint8')
                packet = [packet repmat(typecast(fcmode, 'single'), npts_packet, 1)];
                subdata = packet';
                subdatainfo = whos('subdata');
                fwrite(conn, subdatainfo.bytes+8, 'uint32');
                fwrite(conn, uint32(size(packet)), 'uint32');
                fwrite(conn, subdata(1:end), 'single');
                i=i+1;
                pause(0.001);
                break
            end
        catch err
            failure = true;
            break
        end
    end
    
    if failure
        break
    end
end
fclose(conn);
