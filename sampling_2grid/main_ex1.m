%%%% two grid, full Neumann, ilu smoother, coarse problem initial
%%%% average p=0 is imposed
clear;close all;
% parastart; %%start parapool
%% mesh
q1=0;nmaxbasis=8;maxvalue=10^(4);overelement=0;overcorr=2;fprintf('overcorr is %d\n',overcorr);
maxit=2000;pretol=10^(-6);
high=1;eigvalue_tol=.5;
regularc=0;regvalue=10^(-14);%%%%%%%%%%%%%%%%
n=16;Ny=8;h=.1; Nz=Ny;Nx=Ny;nx=Nx*n;ny=Ny*n;nz=Nz*n;Lx=h*nx;Ly=h*ny;Lz=h*nz;hx=h;hy=h;hz=h;  
nfacezx=(ny+1)*nz*nx;nfaceyx=(nz+1)*ny*nx;nfaceyz=ny*nz*(nx+1);nnx=nx+1;nny=ny+1;nnz=nz+1;
nface=nfacezx+nfaceyx+nfaceyz;ne=nx*ny*nz;ndof=nface+ne;fprintf('nx and n are %d %d\n',nx,n);
tic;[~,~,~,~,~,~,bddof]=getallfacedof(nx,ny,nz);toc;closepara=0;
nlocalface=3*n^2*(n+1);

%% model
K=Construction_3D_perm256(maxvalue,nx/16);  

contrast=max(K(:))/min(K(:));maxk=max(K(:));
fprintf('contrast is %2.2e \n',contrast);

%% fine problem
disp('assemble Avp...');
tfinestart= tic;
nodv=ldof2gdofmixed3d(nx,ny,nz);localnodv=ldof2gdofmixed3d(n,n,n);

localvp=[-hx*hz,hx*hz,-hy*hx,hy*hx,-hy*hz,hy*hz];
ivalue=repmat(localvp',1,ne);
ix=repmat(1:ne,6,1);
tic
Avp=sparse(nodv,ix,ivalue,nface,1+ne);clear ivalue ix 
Avp(bddof,:)=[];
disp('assemble Avv p1...');

diagAvv=sparse(nodv,ones(size(nodv)),ones(6,1)*(1./K(:))')*hx*hy*hz*1/2;
diagAvv=1./diagAvv';
diagAvv(bddof)=[];
clear nodv bddof

disp('assemble Avv p2...');
Aeli=(Avp'.*diagAvv)*Avp;clear Avp diagAvv

% return
disp('assemble Avv p3...');

ix1=(1+ne)*ones(2*ne+1,1);ix1(1:ne)=1:ne;
iy1=(1+ne)*ones(2*ne+1,1);iy1(1+ne:end)=1:ne+1;
ivalue1=ones(2*ne+1,1);
B=sparse(ix1,iy1,ivalue1);clear ix1 iy1 ivalue1
toc

Aeli=Aeli-B;clear Avp diagAvv B
% return
%% force
disp('get force')
force=zeros(ny,nz,nx);sourcevalue=1;
        nyw=floor(ny/2)+1;nzw=1:nz;nxw=floor(nx/2)+1;
    force(1,:,1)=sourcevalue;force(end,:,end)=sourcevalue;
    force(1,:,end)=sourcevalue;force(end,:,1)=sourcevalue;
    force(nyw,nzw,nxw)=-4*sourcevalue;
Feli=zeros(ne+1,1);Feli(1:ne)=-force(:)*hx*hy*hz;clear force_vec forcev 

tfineass=toc(tfinestart);
fprintf('fine assemble cost %2.1f seconds\n',tfineass);

%% ms basis
  tbasisstart=tic; 
tic;phicorr=compute_corr_fun1(n,Nx,Ny,Nz,K,overcorr,hx,hy,hz,force);toc
eli_adaptive_offlinebasis_para_singular;

% return
tbasis=toc(tbasisstart);%%%%%%%%%%%%
fprintf('tmsbasis are %2.1f  \n',tbasis);
%% coarse matrix

clear Avv A Avp F1 doffaceyx doffaceyz doffacezx force freedof 
clear localAvp dofp  bddof alpha K0 p_dof 
disp('form coarse matrix....')

tcms=tic;

Amsa=phimatrix_pa*Aeli*phimatrix_pa';Fmsa=phimatrix_pa*Feli;
if regularc~=0
Amsa=Amsa+regularc*diag(diag(Amsa));
end

    disp('lu for coarse matrix');
tic;[L,U,P,Q] = lu(Amsa);toc
pmsac=Q*(U\(L\(P*(Fmsa))));
pmsa=phimatrix_pa'*pmsac;%%%% ms initial condition
pmsa=zeros(ne+1,1);

disp('fine ilu...')
   tic; [Lilu,Uilu]=ilu(Aeli);toc

clear phicorr phicorrm iivaluepcell  iivalue K  localnodv basis global_p_dof
%% preconditioner
 precond=@(inc)pre2grid2_ilu(Aeli, inc, zeros(ne+1,1),1,1, 1,L,U,P,Q,...
    phimatrix_pa,1,Lilu,Uilu);

tcms=toc(tcms);

disp('preconditioner solve....')
if 1
 [ppre, error, iter, flag,condn,tpre] = pcg0_cond(Aeli, pmsa, Feli, maxit,pretol,precond); 
else
    disp('solve with bicgstab...');
 t2grid1=tic;[ppre,flag,relres2,iter2]=bicgstab(Aeli,Feli,pretol, 5000,precond,[],pmsa);tpre=toc(t2grid1);iter=iter2*2;condn=0;
     disp('solve with gmres...');
 t2grid1=tic;[ppre,flag,relres2,iter2]=gmres(Aeli,Feli,10,pretol, 5000,precond,[],pmsa);tpre=toc(t2grid1);iter=iter2*2;condn=0;
    disp('solve with cg...');
 t2grid1=tic;[ppre,flag,relres2,iter2]=pcg(Aeli,Feli,pretol, 5000,precond,[],pmsa);tpre=toc(t2grid1);iter=iter2*2;condn=0;

end
 ndofc=size(Amsa,1)-2;%norm(Aeli*ppre-Feli)/norm(Feli)
fprintf('tfinemat, tmsbasis, tcoarsemat are %2.1f  & %2.1f & %2.1f\n',tfineass,tbasis,tcms);
fprintf('contrast, dof (Nb), Iter, condn, toff, ton are %2.1e & %d (%2.2f)& %d & %2.1f  & %2.1f & %2.1f\n',contrast,ndofc,ndofc/Nx/Ny/Nz,iter,condn,tfineass+tbasis+tcms,tpre);


%  fprintf('pre solve iter number is %d \n',iter);
    fprintf('contrast  is %2.2e \n',contrast);
if q1==1
    disp('polynomial')
else
    disp('ms basis')
end


 fprintf('dim are %d %d\n',dim_pc,ne);
 fprintf('average basis number is %2.2f\n',dim_pc/(Nx*Ny*Nz));

 disp('...............................END........................................');
