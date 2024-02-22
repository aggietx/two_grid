delete(gcp('nocreate'))
clear parapool
if exist('parpool')~=1
c = parcluster('local');
c.NumWorkers =48/3;
parpool(c, c.NumWorkers);
% matlabpool(c, c.NumWorkers);
parpool.IdleTimeout = 1200000;
end
return
Ncore=12;parpoola = parpool;
 parpoola.IdleTimeout = 12000;
 
%  Ncore=12;parpoola = matlabpool;
%  parpoola.IdleTimeout = 12000;
a=zeros(10,1);
parfor (i=1:10,Ncore)
    a(i)=i;
end
