%open stl file
%function [abt]=slicesuccess(triangles,slice_height)

% clear all
% clc
% close all
% [newbodi,colnormal]=buildunbuild ()      %change
% bdi=newbodi;              %change
 %[vow,holefilledges] = holefillalgo(bdi)  %change
%fid=fopen('cylinderp.stl','r');
%fid=fopen('cstlwr3.stl','r');
% if fid==-1
%     error('error opening file');
% end
% i=1;
% s1='solidascii';
 s2='facetnormal';
% s3='outerloop';
 s4='vertex';
% s5='endloop';
 n1=[];
  v=[];
  v1=[];
while feof(fid)==0
    l=fgetl(fid);
    f1=sscanf(l,'%s%*s%*s%*s');
    c=strcmp(s4,f1);                             %comparing the string wih string s4
    if c~=0         
        f=90;                                  %initiating the variable f to verify whether loop is going into it
        v1=[v1,(sscanf(l,'%*s%f%f%f'))] ;       %extracting the values of vertex
    else
        n1=[n1,(sscanf(l,'%*s%*s%f%f%f'))];     %extracting the values of normal
    end 
  
    disp(l);
    
end
n2=n1';
%v1=v1';
v2 = [];
v3=[];
v4=[];
v5=[];
% v1=holefilledges; %change
% for i=1:1:size(holefilledges,1)
%     v3 = cat(2,v1(i,1),v1(i,2),v1(i,3));
%     v4 = cat(2,v1(i,4),v1(i,5),v1(i,6));
%     v5 = cat(2,v1(i,7),v1(i,8),v1(i,9));
%     v2 = cat(1,v2,v3(1,:),v4(1,:),v5(1,:));
% end
% v1=vow; %change
% for i=1:2:size(vow,1)-1
%     v3 = cat(2,v1(i,1),v1(i,2),v1(i,3));
%     v4 = cat(2,v1(i,4),v1(i,5),v1(i,6));
%     v5 = cat(2,v1(i+1,4),v1(i+1,5),v1(i+1,6));
%     v2 = cat(1,v2,v3(1,:),v4(1,:),v5(1,:));
% end
%v1=triangles; %change 
%v1 = trianglesone;
%for i=1:1:size(triangles,1)
 for i=1:1:size(v1,1)
    v3 = cat(2,v1(i,1),v1(i,2),v1(i,3));
    v4 = cat(2,v1(i,4),v1(i,5),v1(i,6));
    v5 = cat(2,v1(i,7),v1(i,8),v1(i,9));
    v2 = cat(1,v2,v3(1,:),v4(1,:),v5(1,:));
end
%v2=v1';
%v2 = v1; %change
%[r1 k1]=size(n2); 
[r1]=(size(v2)/3); %change
[r2 k2]=size(v2);
%close file
%fclose(fid);
m=max(abs(v2),[],1) %getting 1X3 matrix of maximum values in x,y and z co-ordinates
%sd=input('choose in which direction slicing has to do 1: X-direction, 2: Y-direction, 3: Z-direction:::'); %choosing the slice direction
sd = 3; %change
disp('Maximum height in choosen direction is:'); 
mx=m(:,sd)                               % extracting the maximum value in that direction
%st=input('Enter slice thickness:'); % slice thickness value
%st = slice_height;
st=0.3;
nsl=mx/st;                           % No of slices
z1(:,1)=0:st:mx;                         % Creating layer height matrix
[rl cl]=size(z1);
ec(:,:)=v2(:,sd);%extracting all co-ordinate values in choosen direction
k=1;
ex=[];
figure
for zv=1:rl
    
    j=1;
    for i=1:r1
        v=v2(j:j+2,sd);
        z=z1(zv,:); % Taking the first layer thickness
        miv=min(abs(v));
        mav=max(abs(v));
        
        if (miv<=z)&&(mav>=z)
            cvvx=[];
            cvvy=[];
            cvvz=[];
            ex=[(v2(j:j+2,:))]; %taking the values of x,y and z co-ordinates which satisfy the rule
            if ((((ex(1,3)>z)&&(ex(2,3)<z)))||((ex(2,3)>z)&&(ex(1,3)<z)))
                xi1=(ex(2,1))+(((z-ex(2,3))/(ex(1,3)-ex(2,3)))*(ex(1,1)-ex(2,1)));
                yi1=(ex(2,2))+(((z-ex(2,3))/(ex(1,3)-ex(2,3)))*(ex(1,2)-ex(2,2)));
                cvvx=[xi1];
                cvvy=[yi1];
                cvvz=[z];
            end
            if ((((ex(1,3)>z)&&(ex(3,3)<z)))||((ex(3,3)>z)&&(ex(1,3)<z)));
                xi2=(ex(3,1))+(((z-ex(3,3))/(ex(1,3)-ex(3,3)))*(ex(1,1)-ex(3,1)));
                yi2=(ex(3,2))+(((z-ex(3,3))/(ex(1,3)-ex(3,3)))*(ex(1,2)-ex(3,2)));
                cvvx=[cvvx,xi2];
                cvvy=[cvvy,yi2];
                cvvz=[cvvz,z];
            end
            if ((((ex(2,3)>z)&&(ex(3,3)<z)))||((ex(3,3)>z)&&(ex(2,3)<z)))
                xi3=(ex(3,1))+(((z-ex(3,3))/(ex(2,3)-ex(3,3)))*(ex(2,1)-ex(3,1)));
                yi3=(ex(3,2))+(((z-ex(3,3))/(ex(2,3)-ex(3,3)))*(ex(2,2)-ex(3,2)));
                cvvx=[cvvx,xi3];
                cvvy=[cvvy,yi3];
                cvvz=[cvvz,z];
            end 
            plot3(cvvx,cvvy,cvvz,'r');
            hold on      
            view(3);  axis vis3d
        end
        
        j=j+3;
        i=i+1;
    end
    rl=rl+1;
end
xlabel('X');
ylabel('Y');
zlabel('Z');
abt = cvvx+cvvy+cvvz;
%end