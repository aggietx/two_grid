function [nodv]=ldof2gdofmixed3d(nx,ny,nz)
ne=nx*ny*nz;
% nodp=1:ne;
nodv=zeros(6,ne);

% nfacezx=(ny+1)*nz*nx;
% nfaceyx=(nz+1)*ny*nx;
% doffacezx=reshape(1:(ny+1)*nz*nx,ny+1,nz,nx);bddofzx=doffacezx([1,end],:,:);bddofzx=bddofzx(:);
% doffaceyx=reshape(1:ny*(nz+1)*nx,ny,nz+1,nx)+nfacezx;bddofyx=doffaceyx([1,end],:,:);bddofyx=bddofyx(:);
% doffaceyz=reshape(1:ny*nz*(nx+1),ny,nz,nx+1)+nfacezx+nfaceyx;bddofyz=doffaceyz([1,end],:,:);bddofyz=bddofyz(:);

[doffacezx,doffaceyx,doffaceyz]=getallfacedof(nx,ny,nz);
y1=doffacezx(1:end-1,:,:); y2=doffacezx(2:end,:,:); clear doffacezx
z1=doffaceyx(:,1:end-1,:); z2=doffaceyx(:,2:end,:); clear doffaceyx
x1=doffaceyz(:,:,1:end-1); x2=doffaceyz(:,:,2:end); 
% nodv=[y1(:)';y2(:)';z1(:)';z2(:)';x1(:)';x2(:)'];
nodv(1,:)=y1(:); nodv(2,:)=y2(:); 
nodv(3,:)=z1(:); nodv(4,:)=z2(:); 
nodv(5,:)=x1(:); nodv(6,:)=x2(:); 

% for ix=1:nx
%     for iz=1:nz
%         for iy=1:ny
%             ne=(ix-1)*( nz*ny)+(iz-1)*ny+iy;
%             nodv(1,ne)= doffacezx(iy,iz,ix); nodv(2,ne)= doffacezx(iy+1,iz,ix);
%             nodv(3,ne)= doffaceyx(iy,iz,ix); nodv(4,ne)= doffaceyx(iy,iz+1,ix);
%             nodv(5,ne)= doffaceyz(iy,iz,ix); nodv(6,ne)= doffaceyz(iy,iz,ix+1);
%         end
%     end
% end
% nodp=(1:ne)+max(max(nodv));
% nod=[nodv;nodp];
