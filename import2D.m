function import(f)
global RAW2D 
pwrite= 2;
fname= f;
if exist(fname), disp('loading file..'), eval(['load ' f])
else disp('FILE DOES NOT EXIST'), break, end 
eval(['RAW=' f ';'])
[m,n]= size(RAW);
k=m/n
RAW2D=zeros(n,n,k);
for i=1:k
  RAW2D(:,:,i)=RAW(1+(i-1)*n:n+(i-1)*n,:);
end
