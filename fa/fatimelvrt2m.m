function time_matlab = fatimelvrt2m(time_labviewrt)

time_matlab = (double(bitsrl(time_labviewrt, 19))*2^19/1e6)/1000/60/60/24 + datenum('1904/1/1') - java.util.Date().getTimezoneOffset()/60/24;
