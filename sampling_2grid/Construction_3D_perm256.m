function K=Construction_3D_perm256(maxvalue,N)
K0=ones(16,16,16);

% K0(2:6,2:6,2:6)=maxvalue;
% K0(2:6,10:14,2:6)=maxvalue;
% K0(2:6,2:6,10:14)=maxvalue;
% K0(2:6,10:14,10:14)=maxvalue;
% K0(10:14,2:6,2:6)=maxvalue;
% K0(10:14,10:14,2:6)=maxvalue;
% K0(10:14,2:6,10:14)=maxvalue;
% K0(10:14,10:14,10:14)=maxvalue;

K0(:,2:4,2:4)=maxvalue;
% K0(:,2:14,2:14)=maxvalue;

% K0(:,2:4,5:6)=maxvalue;
% K0(:,2:4,8:11)=maxvalue;
% K0(:,2:4,13:15)=maxvalue;
% 
% K0(:,5:6,2:4)=maxvalue;
K0(:,5:6,5:6)=maxvalue;
% K0(:,5:6,8:11)=maxvalue;
% K0(:,5:6,13:15)=maxvalue;
% K0(:,7:11,7:11)=maxvalue;
% K0(:,8:11,2:4)=maxvalue;
% K0(:,8:11,5:6)=maxvalue;
% K0(:,8:11,8:11)=maxvalue;
% K0(:,8:11,13:15)=maxvalue;
% 
% K0(:,13:15,2:4)=maxvalue;
% K0(:,13:15,5:6)=maxvalue;
% K0(:,13:15,8:11)=maxvalue;
% K0(:,13:15,13:15)=maxvalue;

K=repmat(K0,[N,N,N]);
% K=repmat(K0,[32,32,32]);