function [true_dof]=get_3dtruedof_element(ix,iy,iz,Nx,Ny,Nz,overelement,n0)
%%% try to classify all the oversampling case
%%% assume square coarse elements
overelement=min(overelement,n0);
% nxtotal=n0*Nx;nytotal=n0*Ny;nztotal=n0*Nz;
left=1:n0;
middle=overelement+1:overelement+n0;
right=overelement+1:overelement+n0;
n1=n0+overelement;
n2=n0+overelement*2;
if ix==1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if  iy==1 && iz==1             %%%
      all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);

true_dof=all_dof(left,left,left);

  elseif ismember(iy,2:Ny-1) &&  iz==1 
        all_local_dof=1:n2*n1*n1;
all_dof=reshape(all_local_dof,[n2,n1,n1]);
true_dof=all_dof(middle,left,left);
  elseif iy==Ny && iz==1  
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(right,left,left);

  elseif iy==1 && ismember(iz,2:Nz-1) %%%
        all_local_dof=1:n1*n2*n1;
all_dof=reshape(all_local_dof,[n1,n2,n1]);
true_dof=all_dof(left,middle,left);

  elseif ismember(iy,2:Ny-1) && ismember(iz,2:Nz-1) 
        all_local_dof=1:n2*n2*n1;
all_dof=reshape(all_local_dof,[n2,n2,n1]);
true_dof=all_dof(middle,middle,left);

  elseif iy==Ny && ismember(iz,2:Nz-1)  
        all_local_dof=1:n1*n2*n1;
all_dof=reshape(all_local_dof,[n1,n2,n1]);
true_dof=all_dof(right,middle,left);

  elseif iy==1 && iz==Nz               %%%
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(left,right,left);

  elseif ismember(iy,2:Ny-1)  &&  iz==Nz
        all_local_dof=1:n2*n1*n1;
all_dof=reshape(all_local_dof,[n2,n1,n1]);
true_dof=all_dof(middle,right,left);

  elseif iy==Ny && iz==Nz
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(right,right,left);
  end

elseif ix==Nx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if  iy==1 && iz==1                %%%
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(left,left,right);

  elseif ismember(iy,2:Ny-1) &&  iz==1
        all_local_dof=1:n2*n1*n1;
all_dof=reshape(all_local_dof,[n2,n1,n1]);
true_dof=all_dof(middle,left,right);
  elseif iy==Ny && iz==1  
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(right,left,right);

  elseif iy==1 && ismember(iz,2:Nz-1) %%%
        all_local_dof=1:n1*n2*n1;
all_dof=reshape(all_local_dof,[n1,n2,n1]);
true_dof=all_dof(left,middle,right);

  elseif ismember(iy,2:Ny-1) && ismember(iz,2:Nz-1)
        all_local_dof=1:n2*n2*n1;
all_dof=reshape(all_local_dof,[n2,n2,n1]);
true_dof=all_dof(middle,middle,right);

  elseif iy==Ny && ismember(iz,2:Nz-1)  
        all_local_dof=1:n1*n2*n1;
all_dof=reshape(all_local_dof,[n1,n2,n1]);
true_dof=all_dof(right,middle,right);

  elseif iy==1 && iz==Nz            %%%
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(left,right,right);

  elseif ismember(iy,2:Ny-1)  &&  iz==Nz
        all_local_dof=1:n2*n1*n1;
all_dof=reshape(all_local_dof,[n2,n1,n1]);
true_dof=all_dof(middle,right,right);

  elseif iy==Ny && iz==Nz
        all_local_dof=1:n1*n1*n1;
all_dof=reshape(all_local_dof,[n1,n1,n1]);
true_dof=all_dof(right,right,right);
  end

else %%% middle x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  if  iy==1 && iz==1            %%%
      all_local_dof=1:n1*n1*n2;
all_dof=reshape(all_local_dof,[n1,n1,n2]);
true_dof=all_dof(left,left,middle);

  elseif ismember(iy,2:Ny-1) &&  iz==1 
      all_local_dof=1:n2*n1*n2;
all_dof=reshape(all_local_dof,[n2,n1,n2]);
true_dof=all_dof(middle,left,middle);
  elseif iy==Ny && iz==1  
      all_local_dof=1:n1*n1*n2;
all_dof=reshape(all_local_dof,[n1,n1,n2]);
true_dof=all_dof(right,left,middle);

  elseif iy==1 && ismember(iz,2:Nz-1) %%%
      all_local_dof=1:n1*n2*n2;
all_dof=reshape(all_local_dof,[n1,n2,n2]);
true_dof=all_dof(left,middle,middle);

  elseif ismember(iy,2:Ny-1) && ismember(iz,2:Nz-1) 
      all_local_dof=1:n2*n2*n2;
all_dof=reshape(all_local_dof,[n2,n2,n2]);
true_dof=all_dof(middle,middle,middle);

  elseif iy==Ny && ismember(iz,2:Nz-1)  
      all_local_dof=1:n1*n2*n2;
all_dof=reshape(all_local_dof,[n1,n2,n2]);
true_dof=all_dof(right,middle,middle);

  elseif iy==1 && iz==Nz            %%%
      all_local_dof=1:n1*n1*n2;
all_dof=reshape(all_local_dof,[n1,n1,n2]);
true_dof=all_dof(left,right,middle);

  elseif ismember(iy,2:Ny-1)  &&  iz==Nz
      all_local_dof=1:n2*n1*n2;
all_dof=reshape(all_local_dof,[n2,n1,n2]);
true_dof=all_dof(middle,right,middle);

  elseif iy==Ny && iz==Nz
      all_local_dof=1:n1*n1*n2;
all_dof=reshape(all_local_dof,[n1,n1,n2]);
true_dof=all_dof(right,right,middle);
  end


end
true_dof=true_dof(:);
% [nylocal,nzlocal,nxlocal]=size(all_dof);