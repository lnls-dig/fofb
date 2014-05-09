function time_labview = time_matlab2labviewrt(time_matlab)

time_matlab = time_matlab - datenum('1904/1/1') + java.util.Date().getTimezoneOffset()/60/24 ;

time_labview = uint64(time_matlab*24*60*60*1e9);