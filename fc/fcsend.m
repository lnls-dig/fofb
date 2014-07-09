function data = fcsend(ip_address, fcn, npts_packet, ip_port, stopat)

if nargin < 3 || isempty(npts_packet)
    npts_packet = 1000;
end

if nargin < 4 || isempty(ip_port)
    ip_port = 3604;
end

if nargin < 5 || isempty(stopat)
    stopat = npts_packet*100;
end

nvars = fcn('size');

i=0;
failure = false;

ncols = 50;
data = zeros(npts_packet*stopat, ncols);

conn = tcpip(ip_address, ip_port, 'OutputBufferSize', 10*npts_packet*(nvars+1)*4);
fopen(conn);
while i < stopat
    [subdata, fcmode] = fcn(i, npts_packet);

    data((1+i*npts_packet):(i+1)*npts_packet, :) = subdata; 
    
    % Ensure data is encoded in 4-byte floating point representation
    subdata = single(subdata);

    while true
        try
            if fread(conn, 1, 'uint8')
                subdata = [subdata repmat(typecast(fcmode, 'single'), npts_packet, 1)];
                subdatat = subdata';
                subdatainfo = whos('subdatat');

                fwrite(conn, subdatainfo.bytes+8, 'uint32');
                fwrite(conn, uint32(size(subdata)), 'uint32');
                fwrite(conn, subdatat(1:end), 'single');

                i=i+1;
                break
            end
        catch
            failure = true;
            break
        end
        pause(0.001);
    end

    if failure
        break
    end
end
fclose(conn);
