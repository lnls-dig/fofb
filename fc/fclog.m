function fclog(filename, expinfo, npts_packet, stopat)

fileid = fopen(filename, 'a+');
fprintf(fileid, ['[' fatimestr(now) '] ']);
fprintf(fileid, 'excitation = %s; amplitude = %g; duration = %d; mode = %s; uncorrelated = %d; npts_packet = %d; stopat = %d', expinfo.excitation, expinfo.amplitude, expinfo.duration, expinfo.mode, expinfo.uncorrelated, npts_packet, stopat);
fprintf(fileid, '\n');
fclose(fileid);