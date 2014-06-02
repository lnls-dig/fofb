function time_matlab = time_labviewrt2matlab(time_labview)

time_matlab = (double(bitsrl(time_labview, 19))*2^19/1e6)/1000/60/60/24 + datenum('1904/1/1') - java.util.Date().getTimezoneOffset()/60/24;
