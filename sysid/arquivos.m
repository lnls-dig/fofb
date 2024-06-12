clc;clear;
lista_H = dir('*H.mat');
lista_V = dir('*V.mat');
ncorr_h=size(lista_H,1);
ncorr_v=size(lista_V,1);
ncorr=ncorr_h+ncorr_v;

sys=cell(ncorr, 1);

for k=1:ncorr_h
    arquivos_H(k,:) = lista_H(k).name;
end
for k=1:ncorr_v
    arquivos_V(k,:) = lista_V(k).name;
end
for i=1:ncorr_h
    arquivo=arquivos_H(i,:);
    [sys_it,fit_it] = modelo(arquivo, [10 10 2], "X");
    sys{i}=sys_it;
end
for i=1:ncorr_v
    arquivo=arquivos_V(i,:);
    [sys_it,fit_it] = modelo(arquivo, [10 10 2], "Y");
    sys{i+ncorr_h}=sys_it;
end
save('sistemas','sys');