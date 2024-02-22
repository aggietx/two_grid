
function perm=Construction_3D_perm_64_1(maxval)

perm=ones(64,64,64);
indexs=[3,6];

for iy=1:2:7
for iz=1:2:7
for ix=1:2:7
aa=[(iy-1)*8+indexs; (ix-1)*8+indexs;(iz-1)*8+indexs];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

end
end
end

for iy=[2,6]
for iz=[2,6]
aa=[(iy-1)*8+indexs; 6,59;(iz-1)*8+indexs];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

end
end

for iy=[2,6]
for ix=[2,6]
aa=[(iy-1)*8+indexs; (ix-1)*8+indexs;6,59];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

end

end

for iz=[2,6]
for ix=[2,6]
aa=[6,59; (ix-1)*8+indexs;(iz-1)*8+indexs];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

end
end

aa=[6,59;30,34;30,34];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

aa=[30,34;30,34;6,59];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

aa=[30,34;6,59;30,34];%1
perm(aa(1,1):aa(1,2),aa(3,1):aa(3,2),aa(2,1):aa(2,2))=maxval;

