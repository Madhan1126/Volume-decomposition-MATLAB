function [vow,holefilledges] = holefillalgo(bdi)         %[vow,cd] 
% clc 
% clear all
% close all
% impmodel = 'holsq.stl';
% [nn,vv,n,v] = readasciifilestl(impmodel)
% [newbodi,colnormal]=buildunbuild (nn,vv,n,v) 

%% hole detection
%bdi=newbodi;
count=0;
vow=[];
abcnew=[];
chd=0;
for i=1:3:size(bdi,1)
abc= cat(1,bdi((i:i+1),:), bdi((i+1:i+2),:), bdi([i+2 i],:));
 for k=1:2:5
    abcnew = cat(2,abc(k,:),abc(k+1,:));
    vow=cat(1,vow,abcnew(1,:));%collected edges
 end
end
detedge=[]; 
ds=0;
for i=1:1:size(vow,1)
    for j=1:3:size(bdi,1)
        abc= cat(1,bdi((j:j+1),:), bdi((j+1:j+2),:),bdi([j+2 j],:));
for k=1:2:5
         abcnew = cat(2,abc(k,:),abc(k+1,:));
   if ((vow(i,1)-abcnew(1,4))==0)&((vow(i,2)-abcnew(1,5))==0)&((vow(i,3)-abcnew(1,6))==0)&((vow(i,4)-abcnew(1,1))==0)&((vow(i,5)-abcnew(1,2))==0)&((vow(i,6)-abcnew(1,3))==0)
        detedge=cat(1,detedge,abcnew(1,:));
   end
end
end
end
pd=0;
ad=0;
bd=0;
gapedge=[];
dt=[];
for i=1:1:size(vow,1)
    dt=vow(i,:);
    pd=pd+1;
    for j=1:1:size(detedge,1)
       if ((dt(1,1)-detedge(j,1))==0) & ((dt(1,2)-detedge(j,2))==0) & ((dt(1,3)-detedge(j,3))==0)&((dt(1,4)-detedge(j,4))==0) & ((dt(1,5)-detedge(j,5))==0) & ((dt(1,6)-detedge(j,6))==0)
            ad=pd;
       else
            bd=pd;
       end
       end
       if bd>ad
           gapedge = cat(1,gapedge,dt(1,:)); %% correctly detecting the gaps
       else
           gapedge=gapedge;
       end
end

%% sorting
j=1;
empty=[0 0 0 0 0 0];
sortedge=[];
nextedge=[];
newedge = [];
firstedge=gapedge(1,:);
first=gapedge(1,:);
sortvertices=[];
sortvertice=[];
filledge=[];
filledges=[];
trsort=[];
trsorted=[];
% convert this sorting into function

for j=1:1:size(gapedge,1)
    %firstedge=firstedge;
for i=1:1:size(gapedge,1)
    if ((firstedge(1,4)-gapedge(i,1))==0)&((firstedge(1,5)-gapedge(i,2))==0)&((firstedge(1,6)-gapedge(i,3))==0)
         nextedge=gapedge(i,:);  
    end
end
    sortedge=cat(1,sortedge,nextedge(1,:));
        sortvertice = cat(2,nextedge(:,1),nextedge(:,2),nextedge(:,3));
        sortvertices = cat(1,sortvertices,sortvertice(1,:));
    n= size(sortedge,1);
if ((first(1,1)-nextedge(1,1))==0)&((first(1,2)-nextedge(1,2))==0)&((first(1,3)-nextedge(1,3))==0)
    if n < size(gapedge,1)
    first = gapedge(n+1,:);
    else
        break
    end
    if size(sortedge,1) == size(gapedge,1)
    else
    m= size(sortvertices,1)
    js=1;
    for i=1:m
        filledge = cat(2,sortvertices(js+1,:),sortvertices(m,:));
        filledges = cat(1,filledges,filledge(1,:));
        trsort= cat(2,sortvertices(js,:),sortvertices(js+1,:),sortvertices(m,:));
        trsorted = cat(1,trsorted,trsort(1,:));
        if mod(i,2) == 0
            js=js+1;
        else
            m = m-1; 
        end
        if m == js+2
        trsort= cat(2,sortvertices(js+1,:),sortvertices(js+3,:),sortvertices(m,:));
        trsorted = cat(1,trsorted,trsort(1,:));
            break
        end
    end
        ot = size(sortedge,1)+1;
        newedge = gapedge(ot,:);
        sortvertices=[];
        sortvertice=[];
       
  for i=1:1:size(sortedge,1)
    if (newedge(1,1)-sortedge(i,1)==0)&(newedge(1,2)-sortedge(i,2)==0)&(newedge(1,3)-sortedge(i,3)==0)&(newedge(1,4)-sortedge(i,4)==0)&(newedge(1,5)-sortedge(i,5)==0)&(newedge(1,6)-sortedge(i,6)==0)     
    else
        firstedge = newedge(1,:);
    end
   end
    end
else
    firstedge = nextedge(1,:);
end
end
m= size(sortvertices,1)
js=1;
    for i=1:m
        filledge = cat(2,sortvertices(js+1,:),sortvertices(m,:));
        filledges = cat(1,filledges,filledge(1,:));
        trsort= cat(2,sortvertices(js,:),sortvertices(js+1,:),sortvertices(m,:));
        trsorted = cat(1,trsorted,trsort(1,:));
        if mod(i,2) == 0
           js=js+1;
        else
            m = m-1; 
        end
        if m == js+2
        trsort= cat(2,sortvertices(js+1,:),sortvertices(js+3,:),sortvertices(m,:));
        trsorted = cat(1,trsorted,trsort(1,:));
            break
        end
    end

holefilledge=[];
% for i=1:1:size(detedge,1)
%     holefilledge = cat(1,holefilledge,detedge(i,:));
% end
% for i=1:1:size(sortedge,1)
%     holefilledge = cat(1,holefilledge,sortedge(i,:));
% end
% for i=1:1:size(filledges,1)
%     holefilledge = cat(1,holefilledge,filledges(i,:));
% end

% for i=1:1:size(vow,1)
%      holefilledge = cat(1,holefilledge,vow(i,:));
%  end
% for i=1:1:size(filledges,1)
%      holefilledge = cat(1,holefilledge,filledges(i,:));
%  end
for i=1:1:size(trsorted,1)
     holefilledge = cat(1,holefilledge,trsorted(i,:));
 end
 holefilledges = unique(holefilledge,'rows','stable');
% %to get vertices multiples of 3
% n=size(holefilledges,1)
% if mod(n,3)==0
%     holefilledges = holefilledges;
% else
%     holefilledges= cat(1,holefilledges,holefilledges(n,:));
% end
% n=size(holefilledges,1)
% if mod(n,3)==0
%      holefilledges = holefilledges;
%     else
%      holefilledges= cat(1,holefilledges,holefilledges(n,:));
% end

figure

for i=1:1:size(vow,1)
  plot3 (vow(i,[1 4]),vow(i,[2 5]),vow(i,[3 6]))
  hold on
end
title('unique  edges')

figure

for i=1:1:size(detedge,1)
  plot3 (detedge(i,[1 4]),detedge(i,[2 5]),detedge(i,[3 6]))
  hold on
end
title('correct edges')

figure
for i=1:1:size(gapedge,1)
  plot3 (gapedge(i,[1 4]),gapedge(i,[2 5]),gapedge(i,[3 6]))
  hold on
end
title('extracted edges')

figure
for i=1:1:size(sortedge,1)
  plot3 (sortedge(i,[1 4]),sortedge(i,[2 5]),sortedge(i,[3 6]))
  hold on
end
for i=1:1:size(filledges,1)
  plot3 (filledges(i,[1 4]),filledges(i,[2 5]),filledges(i,[3 6]))
  hold on
end
title('Filled edges')

% figure
% for i=1:1:size(detedge,1)
%   plot3 (detedge(i,[1 4]),detedge(i,[2 5]),detedge(i,[3 6]))
%   hold on
% end
% title('Filled edges')
% for i=1:1:size(sortedge,1)
%   plot3 (sortedge(i,[1 4]),sortedge(i,[2 5]),sortedge(i,[3 6]))
%   hold on
% end
% for i=1:1:size(filledges,1)
%   plot3 (filledges(i,[1 4]),filledges(i,[2 5]),filledges(i,[3 6]))
%   hold on
% end

% figure
% for i=1:1:size(holefilledges,1)
%   plot3 (holefilledges(i,[1 4]),holefilledges(i,[2 5]),holefilledges(i,[3 6]))
%   hold on
% end
% title('HoleFilled geometry')


%% centroid
%               P=filledges;
%         %polyin = polyshape(P);
%         [cx,cy] = centroid(P);
end