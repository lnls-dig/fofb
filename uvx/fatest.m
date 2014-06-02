path = '/media/fofb-archiver';
date_start = '2014/5/9 8:59:31';
date_stop  = '2014/5/9 9:00:00';

filenames = fafind(path, date_start, date_stop)

fa_data_9_p1 = faload(filenames)
