function  x=pre2grid2_ilu(Apre, Fpre, x,step1,step2, omegadinvA,L,U,P,Q,...
   phimatrix_enriched,ncycle,Lilu,Uilu)
for ii=1:ncycle
%    for k = 1:step1
%     x = x + omegadinvA.*(Fpre - Apre*x);
%    end

%  for k = 1:step1
%    r=Fpre - Apre*x;
%    x=x+Uilu\(Lilu\r);


   x=Uilu\(Lilu\Fpre);%%%% input x is zero
%  end
%  for k=1:2
%    r=Fpre - Apre*x;
%    x=x+Uilu\(Lilu\r);
%  end    
 
   r=Fpre - Apre*x;
   r0=phimatrix_enriched*r;     
% %    uc=Ams1\(Ams2\r0);  
   uc=Q*(U\(L\(P*(r0))));
   umsfine=phimatrix_enriched.'*uc; 
   x=x + umsfine;
   
% r0=phimatrix_enriched*(Fpre - Apre*x);  
%   x=x + phimatrix_enriched.'*(Q*(U\(L\(P*(r0)))));         
%  for k = 1:step2         
    r=Fpre - Apre*x;
   x=x+Uilu\(Lilu\r);
%  end
%   for k=1:2
%    r=Fpre - Apre*x;
%    x=x+Uilu\(Lilu\r);
%   end  
% %    for k = 1:step2
% %     x = x + omegadinvA.*(Fpre - Apre*x);
% %    end

end