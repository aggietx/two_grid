function [A]=assemble_fd_fullneu2(nx,ny,nz,h,ne,K,ix,iy)
%%%% do not include augment matrix, ix iy are known
pdof=reshape(1:nx*ny*nz,ny,nz,nx);
% ix=ones(ne,7);
% iy=repmat( (1:ne)',1,7);
ivalue=zeros(ne,7);

%% y+1
freedof=pdof(1:end-1,:,:);freedof=freedof(:);
freedof1=pdof(2:end,:,:);freedof1=freedof1(:);
cent=1./K(freedof);
% ix(freedof,1)=freedof1;
ivalue(freedof,1)=2./(cent+1./K(freedof1 ));
%% y-1
freedof=pdof(2:end,:,:);freedof=freedof(:);
freedof1=pdof(1:end-1,:,:);freedof1=freedof1(:);
cent=1./K(freedof);
% ix(freedof,2)=freedof1;
ivalue(freedof,2)=2./(cent+1./K(freedof1 ));
%% z+1
freedof=pdof(:,1:end-1,:);freedof=freedof(:);
freedof1=pdof(:,2:end,:);freedof1=freedof1(:);
% ix(freedof,3)=freedof1;
cent=1./K(freedof);
ivalue(freedof,3)=2./(cent+1./K(freedof1 ));
%% z-1
freedof=pdof(:,2:end,:);freedof=freedof(:);
freedof1=pdof(:,1:end-1,:);freedof1=freedof1(:);
% ix(freedof,4)=freedof1;
cent=1./K(freedof);
ivalue(freedof,4)=2./(cent+1./K(freedof1 ));
%% x+1
freedof=pdof(:,:,1:end-1);freedof=freedof(:);
freedof1=pdof(:,:,2:end);freedof1=freedof1(:);
ix(freedof,5)=freedof1;
cent=1./K(freedof);
ivalue(freedof,5)=2./(cent+1./K(freedof1 ));
%% x-1
freedof=pdof(:,:,2:end);freedof=freedof(:);
freedof1=pdof(:,:,1:end-1);freedof1=freedof1(:);
% ix(freedof,6)=freedof1;
cent=1./K(freedof);
ivalue(freedof,6)=2./(cent+1./K(freedof1));

% ix(:,7)=1:ne;
ivalue(:,7)=-sum(ivalue(:,1:6),2);


A=sparse(iy,ix,-ivalue*h,ne,ne);
%A=fsparse(iy,ix,-ivalue*h,[ne,ne]);
clear ix iy ivalue
