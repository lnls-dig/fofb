function fadata_array = fafcexps(fadata, start, npts_exp, npts_interval)

fadata_array{1} = facutdata(fadata, start-npts_exp + (0:npts_exp-1),[],[]);
fadata_array{2} = facutdata(fadata, start + (0:npts_exp-1),[],[]);

for i=1:42
    fadata_array{i+2} = facutdata(fadata, i*(npts_exp +npts_interval) + start + (0:npts_exp-1),[],i);
end
