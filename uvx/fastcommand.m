function fastcommand(ip_address, data, packetsize, pauseperiod, ip_port)

if nargin < 3 || isempty(packetsize)
    packetsize = 1000;
end

if nargin < 4 || isempty(pauseperiod)
    pauseperiod = 0.01;
end

if nargin < 5 || isempty(ip_port)
    ip_port = 4007;
end

[npts, ncols] = size(data);
npackets = ceil(npts/packetsize);

conn = tcpip(ip_address, ip_port, 'OutputBufferSize', 10*packetsize*ncols*4);
fopen(conn);
i=1;
while i<=npackets
    try
        if fread(conn, 1, 'uint8')
            if i == npackets
                subdata = data(((npackets-1)*packetsize+1):end, :);
            else
                subdata = data((1:packetsize) + (i-1)*packetsize, :);
            end
            subdatainfo = whos('subdata');
            subdatat = subdata';
            fwrite(conn, subdatainfo.bytes+8, 'uint32');
            fwrite(conn, uint32(size(subdata)), 'uint32');
            fwrite(conn, subdatat(1:end), 'single');
            i=i+1;
        end
    catch
        fclose(conn);
        break
    end
    pause(pauseperiod);
end
fclose(conn);