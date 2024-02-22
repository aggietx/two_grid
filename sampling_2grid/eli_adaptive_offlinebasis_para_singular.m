disp('compute ms spectral basis...');
%%%% for singular source
basismap=nmaxbasis*ones(Ny,Nz,Nx);
% % kappa=K*10^6;
overelement=0;
p_dof=reshape(1:nx*ny*nz,ny,nz,nx);
Ncb=Nx*Ny*Nz;
alllocalbddodf=cell(Ncb,1);
alltruedof=cell(Ncb,1);
 alllocalsize=zeros(3,Ncb);
 iivaluepcell=cell(Ncb,1);
 alllocalk=cell(Ncb,1);
 alllocal_pdof=cell(Ncb,1);
  alllocalmasscoeff=cell(Ncb,1);
for iixp=1:Nx
    for iiz=1:Nz
        for iiyp=1:Ny
    iie=(iixp-1)*Nz*Ny+(iiz-1)*Ny+iiyp;
    localk=K( max( (iiyp-1)*n+1-overelement,1):min(iiyp*n+overelement,ny),...
        max( (iiz-1)*n+1-overelement,1):min(iiz*n+overelement,nz),...
    max( (iixp-1)*n+1-overelement,1):min(iixp*n+overelement,nx));
[localny,localnz,localnx]=size(localk);
    alllocalsize(:,iie)=[localny,localnz,localnx];
% [true_dof]=get_3dtruedof_element(iixp,iiyp,iiz,Nx,Ny,Nz,overelement,n);
             global_p_dof=p_dof( (iiyp-1)*n+1:iiyp*n,(iiz-1)*n+1:iiz*n,(iixp-1)*n+1:iixp*n);
             global_p_dof=global_p_dof(:);
% alltruedof{iie}=true_dof;    
alllocalk{iie}=localk;
alllocal_pdof{iie}=global_p_dof;
% alllocalmasscoeff{iie}=ktilda((iiyp-1)*n+1:iiyp*n,(iiz-1)*n+1:iiz*n,(iixp-1)*n+1:iixp*n);
        end
    end
end
% % localnodv=ldof2gdofmixed3d(n,n,n);
% % [~,~,~,~,~,~,localbddof]=getallfacedof(n,n,n);
% %     [nlocalface,nlocalbdface]=getnface(n,n,n);
% % localAvp=assemblemixed3dvp(hx,hy,hz,n,n,n,localnodv);
% % % localAvp(localbddof,:)=0;
% % localAvp(localbddof,:)=[];
[~,~,ixfd,iyfd]=assemble_fd_fullneu1(n,n,n,h,n^3,ones(n,n,n));
%  LASTN = maxNumCompThreads(1);
tic
parfor iie=1:Nz*Ny*Nx
% for iie=1:Nz*Ny*Nx
%     localF=zeros(nlocalface,nlocalbdface);
% locals=alllocalsize(:,iie);
localny=n;
localnz=n;
localnx=n;

nlocale=localnx*localny*localnz;


localk=alllocalk{iie};
localAeli0=assemble_fd_fullneu2(localnx,localny,localnz,h,nlocale,localk,ixfd,iyfd);


% % localmasscoeff=alllocalmasscoeff{iie};Meig=sparse(1:nlocale,1:nlocale,localmasscoeff(:))/(n^2);
Aeig=(localAeli0+localAeli0.');
Aeig=Aeig+regvalue*diag(diag(Aeig));
Meig=sparse(1:nlocale,1:nlocale,localk(:)/(n^2));
%Meig=sparse(1:nlocale,1:nlocale,diag(Aeig))/(n^2);

Meig=Meig+regvalue*diag(diag(Meig));Meig=(Meig+Meig.');
[eigfun,eigvalue]=eigs(Aeig,Meig,nmaxbasis+1,'sm');
% [eigfun,eigvalue]=eig(full(Aeig),full(Meig));
[d, order] = sort(diag(eigvalue), 'ascend'); 
        eigfun=eigfun(:,order);eigvalue=eigvalue(:,order);
        thredhold=find(d>=eigvalue_tol);%% adaptive local number of basis
       if size(thredhold,2)*size(thredhold,1)==0  %%%% all eigvalue less than tolerance
        nlocalbasis=nmaxbasis;
       else
       nlocalbasis=thredhold(1)-1;
       end 
       basis=eigfun(:,1:nlocalbasis);
       basismap(iie)=nlocalbasis;
alleigval{iie}=diag(eigvalue);
%          [eigvalue, order] = sort(diag(D), 'ascend'); %%%%%no oversampling part
%          eigfun=eigfun(:,order);basis=eigfun(:,1:noffline_basis);   
% % true_dof=alltruedof{iie};basis=basis(true_dof,:);
% % allbasis(:,1:nb,iie)=basis; 
iivaluepcell{iie}=basis;
end
toc
dim_pc=sum(basismap(:));
ixp=zeros(dim_pc+2,n^3);ixp(end-1:end,:)=ne+1;
ivaluep=zeros(dim_pc+2,n^3);
ivaluep(end,end)=1/(n^3);
% iyp=repmat((1:dim_pc)',1,n^3);
iyp=(1:dim_pc+2)'*ones(1,n^3);
ind=1;
for iixp=1:Nx
    for iiz=1:Nz
        for iiyp=1:Ny
    iie=(iixp-1)*Nz*Ny+(iiz-1)*Ny+iiyp;
    basis=iivaluepcell{iie};
    nlocal_basis=basismap(iie);
% %     allbasis(:,1:nlocal_basis,iie)=basis; 
    ivaluep(ind:ind+nlocal_basis-1,:)=basis';   
    global_p_dof=alllocal_pdof{iie};
    ixp(ind:ind+nlocal_basis-1,:)= repmat( global_p_dof,1,nlocal_basis)';
% %     ixp(ind:ind+nlocal_basis-1,:)=  (global_p_dof*ones(1,nlocal_basis))';
ind=ind+nlocal_basis;  
        end
    end
end
% if max(imag(allbasis(:)))>10^(-8)
%     disp('complex eigenfunction')
% end
% return
disp('assembling coarse matrix and add correction...')
clear alllocal_pdof allbasis p_dof
phicorrm=sparse((1+dim_pc)*ones(ne,1),(1:ne)',phicorr,dim_pc+2,ne+1);
 phimatrix_pa=sparse(iyp,ixp,ivaluep,dim_pc+2,ne+1)+phicorrm;clear ixp iyp ivaluep
 fprintf('dim of the coarse system is %d\n',dim_pc);
 fprintf('average basis number is %2.2f\n',dim_pc/(Nx*Ny*Nz));
 disp('........')