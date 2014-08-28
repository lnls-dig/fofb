function [datae, datav] = idpreprocessdata(datae, datav, ui, yi, sysarray, fixedsys, inoutsel)

if nargin < 7 || isempty(inoutsel)
    inoutsel = true;
end

nptse = size(datae,1);
nptsv = size(datav,1);

if inoutsel
    nsig = length(ui);
    sigi = ui;
else
    nsig = length(yi);
    sigi = yi;
end

sigme = zeros(nptse, nsig);
sigmv = zeros(nptsv, nsig);

for i=1:nsig
    if inoutsel
        sige = datae.InputData(:,sigi(i));
        sigv = datav.InputData(:,sigi(i));
    else
        sige = datae.OutputData(:,sigi(i));
        sigv = datav.OutputData(:,sigi(i));
    end
    sys = fixedsys;
    if iscell(sysarray)
        sys = sys*sysarray{sigi(i)};
    else
        sys = sys*sysarray;
    end
    sigme_ = lsim(sys, [sige; sige]);
    sigme(:,i) = sigme_(nptse+1:2*nptse);
    sigmv_ = lsim(sys, [sigv; sigv]);
    sigmv(:,i) = sigmv_(nptsv+1:2*nptsv);
end

if inoutsel
    datae = datae(:,yi,:);
    datav = datav(:,yi,:);
    datae.InputData = sigme;
    datav.InputData = sigmv;
else
    datae = datae(:,:,ui);
    datav = datav(:,:,ui);
    datae.OutputData = sigme;
    datav.OutputData = sigmv;
end