filenames = fafind('/media/fofb-archiver', '2014/07/15 12:59', '2014/07/15 13:12:30');
data = faload(filenames,[],[],100000);
[M, corr_steps, orb_std, corr_std] = farespm(data, 1000);