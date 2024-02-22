function [nbddof,bddof]=getbddof(ny,nz,nx)
%     
nfacezx=(ny+1)*nz*nx;nfaceyx=(nz+1)*ny*nx;
doffacezx=reshape(1:(ny+1)*nz*nx,ny+1,nz,nx);bddofzx=doffacezx([1,end],:,:);bddofzx=bddofzx(:);
doffaceyx=reshape(1:ny*(nz+1)*nx,ny,nz+1,nx)+nfacezx;bddofyx=doffaceyx(:,[1,end],:);bddofyx=bddofyx(:);
doffaceyz=reshape(1:ny*nz*(nx+1),ny,nz,nx+1)+nfacezx+nfaceyx;bddofyz=doffaceyz(:,:,[1,end]);bddofyz=bddofyz(:);
bddof=[bddofzx;bddofyx;bddofyz];
nbddof=size(bddof,1);
