filenames = fafind('/media/DIG/Projetos/Ativos/DT_FOFB_LNLS/Dados_FOFB/dados_FOFB_2014.06.23_11h38_12h03');
data = faload(filenames,[],[],100000);

start = 145494;
npts_exp = 100000;
npts_interval = 1000;

data_array = fafcexps(data, start, npts_exp, npts_interval);