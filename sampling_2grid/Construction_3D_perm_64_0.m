
function perm=Construction_3D_perm_64_0(maxval)

%%%% corresponding to model 1 in 2grid paper

% 
perm=ones(16,16,16);
perm(3:4,3:4,1:end)=maxval;
perm(6,6,1:end)=maxval;
% perm(11:14,11:14,2:end)=maxval;

perm=repmat(perm,4,4,4);




