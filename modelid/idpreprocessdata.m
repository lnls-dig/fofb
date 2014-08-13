function [datae, datav] = idpreprocessdata(datae, datav, ui, yi, sysarray, fixedsys, period)

nu = length(ui);
um = zeros(period, nu);

for i=1:nu
    u = datae.InputData(:,ui(i));
    sys = fixedsys;
    if iscell(sysarray)
        sys = sys*sysarray{ui(i)};
    else
        sys = sys*sysarray;
    end
    um_ = lsim(sys, [u; u]);
    um(:,i) = um_(period+1:2*period);
end

datae = datae(:,yi);
datae.InputData = um;
datav = datav(:,yi);
datav.InputData = um;