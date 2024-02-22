function [doffacezx,doffaceyx,doffaceyz,bddofzx,bddofyx,bddofyz,bddof,freedof,...
    bddofzx1,bddofzx2,bddofyx1,bddofyx2,bddofyz1,bddofyz2]=getallfacedof(nx,ny,nz)
%     
nfacezx=(ny+1)*nz*nx;nfaceyx=(nz+1)*ny*nx;nfaceyz=ny*nz*(nx+1);nface=nfacezx+nfaceyx+nfaceyz;
doffacezx=reshape(1:(ny+1)*nz*nx,ny+1,nz,nx);%bddofzx=doffacezx([1,end],:,:);bddofzx=bddofzx(:);
doffaceyx=reshape(1:ny*(nz+1)*nx,ny,nz+1,nx)+nfacezx;%bddofyx=doffaceyx(:,[1,end],:);bddofyx=bddofyx(:);
doffaceyz=reshape(1:ny*nz*(nx+1),ny,nz,nx+1)+nfacezx+nfaceyx;%bddofyz=doffaceyz(:,:,[1,end]);bddofyz=bddofyz(:);


bddofzx1=doffacezx(1,:,:);bddofzx1=bddofzx1(:);
bddofzx2=doffacezx(end,:,:);bddofzx2=bddofzx2(:);bddofzx=[bddofzx1;bddofzx2];
bddofyx1=doffaceyx(:,1,:);bddofyx1=bddofyx1(:);
bddofyx2=doffaceyx(:,end,:);bddofyx2=bddofyx2(:);bddofyx=[bddofyx1;bddofyx2];
bddofyz1=doffaceyz(:,:,1);bddofyz1=bddofyz1(:);
bddofyz2=doffaceyz(:,:,end);bddofyz2=bddofyz2(:);bddofyz=[bddofyz1;bddofyz2];
bddof=[bddofzx;bddofyx;bddofyz];
freedof=setdiff(1:nface+nx*ny*nz,bddof);



