function phi=compute_corr_fun1(n,Nx,Ny,Nz,kappa,overelement,hx,hy,hz,force)
%%%% compute correction function
disp('compute correction function...')
Ncb=Nx*Ny*Nz;
nx=Nx*n;
ny=Ny*n;
nz=Nz*n;
p_dof=reshape(1:nx*ny*nz,ny,nz,nx);
alltruedof=cell(Ncb,1);
 localvp=[-hx*hz,hx*hz,-hy*hx,hy*hx,-hy*hz,hy*hz];

iivaluecell=cell(Ncb,1);
 alllocalk=cell(Ncb,1);
 alllocal_pdof=cell(Ncb,1);
 alllocalr=cell(Ncb,1);

 
for iix=1:Nx
    for iiz=1:Nz
        for iiy=1:Ny
    iie=(iix-1)*Nz*Ny+(iiz-1)*Ny+iiy;
        localr=force( max( (iiy-1)*n+1-overelement,1):min(iiy*n+overelement,ny),...
        max( (iiz-1)*n+1-overelement,1):min(iiz*n+overelement,nz),...
    max( (iix-1)*n+1-overelement,1):min(iix*n+overelement,nx));
if norm(localr(:))>0
       localk=kappa( max( (iiy-1)*n+1-overelement,1):min(iiy*n+overelement,ny),...
        max( (iiz-1)*n+1-overelement,1):min(iiz*n+overelement,nz),...
    max( (iix-1)*n+1-overelement,1):min(iix*n+overelement,nx));

[true_dof]=get_3dtruedof_element(iix,iiy,iiz,Nx,Ny,Nz,overelement,n);
             global_p_dof=p_dof( (iiy-1)*n+1:iiy*n,(iiz-1)*n+1:iiz*n,(iix-1)*n+1:iix*n);
             global_p_dof=global_p_dof(:);
             alllocal_pdof{iie}=global_p_dof;
alltruedof{iie}=true_dof;    
alllocalk{iie}=localk;

alllocalr{iie}=localr(:);

end
        end
    end
end

parfor iie=1:Nz*Ny*Nx
% for iie=1:Nz*Ny*Nx   
% 
if norm(alllocalr{iie})>0
localr=alllocalr{iie};localr=localr(:);
localk=alllocalk{iie};
 [localny,localnz,localnx]=size(localk);
 

ne=localnx*localny*localnz;
nface=(localnx+1)*localny*localnz+(localny+1)*localnx*localnz+(localnz+1)*localny*localnx;


nodv=ldof2gdofmixed3d(localnx,localny,localnz);
 [~,bddof]=getbddof(localny,localnz,localnx);
ivalue=repmat(localvp',1,ne);
ix=repmat(1:ne,6,1);
Avp=sparse(nodv,ix,ivalue,nface,1+ne);
Avp(bddof,:)=0;

diagAvv=sparse(nodv,ones(size(nodv)),ones(6,1)*(1./localk(:))')*hx*hy*hz*1/2;
diagAvv=1./diagAvv';diagAvv(bddof)=1;
Aeli=(Avp'.*diagAvv)*Avp;
ix1=[1:ne,(1+ne)*ones(1,ne),ne+1];iy1=[(1+ne)*ones(1,ne),1:ne,ne+1];ivalue1=ones(1,2*ne+1);
B=sparse(ix1,iy1,ivalue1);
Aeli=Aeli-B;
Feli=[ localr*hx*hy*hz;0];
localp=Aeli\Feli;localp(end)=[];

true_dof=alltruedof{iie};
iivaluecell{iie}=localp(true_dof);
end

end
phi=zeros(n^3*Nx*Ny*Nz,1);

for iix=1:Nx
    for iiz=1:Nz
        for iiy=1:Ny
             iie=(iix-1)*Nz*Ny+(iiz-1)*Ny+iiy;
             phi(alllocal_pdof{iie})=iivaluecell{iie};
        end
    end
end
