%function [newbodi,colnormal,newubi]=buildunbuild (nn,vv,n,v)
function [bdi,newbdi,ubdtrik]=buildunbuildnew1(nn,vv,n,v)
%[gal]=buildunbuild ()
%[nn,vv,n,v] = readasciifilestl()

nnw=nn';
vvw=vv';

botz=min(vvw(:,3));
topz=max(vvw(:,3));

mg=[];
%% layer thickness
%thickness=input('enter the slice thickness')
% thickness=0.2;
% %sd=input('1.X-direction, 2.Y-direction, 3.Z-direction'); %directionofslice
% sd=3;
% minimumvalue=min(vv(sd,:))
% maximumvalue=max(vv(sd,:))
% zlayer=minimumvalue
% count =0;
% [nrow,ncol]=size(vv(3,:));
% vu=0;
% for i=zlayer:thickness:maximumvalue
%     count=count+1;
% end
% % for i=0:thickness:maximumvalue
% %     x_centroid(i)=xy(1,i).centroid(1);
% %     y_centroid(i)=xy(2,i).centroid(2);
% %     z_centroid(i)=xy(3,i).centroid(3);
% % end
% disp('the number of layers is')
% disp(count)


 %% possible edge collection
%take first triangle set, compare with other sets, take 2 common points set
%compare its normals for separating the edges.

 kl=n+v;
 u=1;
 x=1;
      for j=1:kl
        mg= cat(1,mg,nnw(u,:),vvw((x:x+2),:));
        %has to join all the normals and its corresponding vertices in one loop
         u=u+1;
        x=x+3;
        if (u>n)&&(x>v)
            break
        end
      end
 exv=[];
for j=2:4:kl
    for l=1:1:4
      exv= cat(1,exv,mg(l,:));    %possible combinations for edge 1-2,2-3,3-1
        if (l>4)
            break
        end
    end
end

 %% test
cnt=0;
edf=[];
exx=[];
bd =[0 0 1];
ubld=[];
angdot=0;
angcolc=[];
chkubik=[];
ubdtrik=[];
bdi=[];
newbodi=[];
colnormal=[];
chkunbi=[];
%take the edges and normal first 
%then remove the duplicated vertices
 for i=2:4:kl
         exx=mg((i:i+2),:);
         edf=mg((i-1),:);  
         ubd = dot(bd,edf);
         angdot= acosd(ubd);                               %this is for dot product
         ubld=cat(1,ubld,ubd(1,:));
         angcolc=cat(1,angcolc,angdot(1,:)); 
         %idea: downward facing normals are negative, extract its triangles
         % taking that as base, extract the unbuildable features.
    %check the if condition once and its logic
%             if angdot > 90
%                  ubdtrik=cat(1,ubdtrik,mg((i:i+2),:));        %unbuildable triangles               
%                 %chkubik=cat(1,chkubik,angcolc(k,:));
%             end
%              if angdot <= 90
%                  bdi=cat(1,bdi,mg((i:i+2),:));                %buildable triangles  
%                 %chkunbi=cat(1,chkunbi,angcolc(k,:));
%              end
if (i>kl)
            break
        end
 end
 %% separation
 % this logic wrong because the first data set collects in first loop 
 %there is no effect of for loop. this is waste
 
 %% new section
 
 j=2;
chkubik=[];
ubdtrik=[];
bdi=[];
chkunbi=[];
nr=[];
 for k=1:1:size(angcolc,1)
            if angcolc(k,:) >  108 %108  %96 also working fine %108    %standard 108
                 ubdtrik=cat(1,ubdtrik,mg((j:j+2),:)); %unbuildable triangles               
                chkubik=cat(1,chkubik,angcolc(k,:));
                j=j+4;
            end
             if angcolc(k,:) <= 108
                 nr=cat(1,nr,mg(j-1,:));
                 bdi=cat(1,bdi,mg((j:j+2),:)); %buildable triangles  
                chkunbi=cat(1,chkunbi,angcolc(k,:));
                j=j+4;
            end
            if k>size(angcolc,1)
            break
    end  
 end
 %% plotting build and unbuild faces
 
 figure
 view(3); camlight; axis off
 axis equal
    pq=1;
 for i=1:size(ubdtrik,1)
        va=ubdtrik(pq:pq+2,:);
        f=[1 2 3];
        p=patch('Faces',f,'Vertices',va);
        set(p,'facecolor','m','facealpha',0.5);
        set(p,'Edgecolor','k','linewidth',2);
        hold on
        pq=pq+3;
        i=i+1;
        if pq>(size(ubdtrik,1))
            break 
        end
 end
%  
figure
view(3); camlight; axis off
axis equal
    qp=1;
 for i=1:size(bdi,1)
        vb=bdi(qp:qp+2,:);
        f=[1 2 3];
        p=patch('Faces',f,'Vertices',vb);
        set(p,'facecolor','m','facealpha',0.5);
        set(p,'Edgecolor','k','linewidth',2);
        hold on
        qp=qp+3;
        i=i+1;
        if qp> size(bdi,1)
            break 
        end
 end
 
 %% attach the base part to the buildable volume
newbdi=[];
for i=1:1:size(ubdtrik,1)
   if ubdtrik(i,3)==botz
       newbdi = cat(1,newbdi,ubdtrik(i,:));
   end
end

newbdi = cat(1,newbdi,bdi);
newbodi=[];
for i=1:1:size(newbdi,1)
   if newbdi(i,3)==topz
       newbdi = newbdi;
   else
       newbodi = cat(1,newbodi,newbdi(i,:));
   end
end


figure
view(3); camlight; axis off
axis equal
    qp=1;
 for i=1:size(newbodi,1)
        vb=newbodi(qp:qp+2,:);
        f=[1 2 3];
        p=patch('Faces',f,'Vertices',vb);
        set(p,'facecolor','m','facealpha',0.5);
        set(p,'Edgecolor','k','linewidth',2);
        hold on
        qp=qp+3;
        i=i+1;
        if qp> size(newbodi,1)
            break 
        end
 end
 
 
 %% remove the top from buildable and move to unbuildable
 
 newunbdi=[];
for i=1:1:size(bdi,1)
   if bdi(i,3)==topz
       newunbdi = cat(1,newunbdi,bdi(i,:));
   end
end

newunbdi = cat(1,newunbdi,ubdtrik);
newunbodi=[];
for i=1:1:size(newunbdi,1)
   if newunbdi(i,3)==botz
       newunbdi = newunbdi;
   else
       newunbodi = cat(1,newunbodi,newunbdi(i,:));
   end
end

figure
view(3); camlight; axis off
axis equal
    qp=1;
 for i=1:size(newunbodi,1)
        vb=newunbodi(qp:qp+2,:);
        f=[1 2 3];
        p=patch('Faces',f,'Vertices',vb);
        set(p,'facecolor','m','facealpha',0.5);
        set(p,'Edgecolor','k','linewidth',2);
        hold on
        qp=qp+3;
        i=i+1;
        if qp> size(newunbodi,1)
            break 
        end
 end
end

 %% removing the base part
% newubi=[];
% for i=1:1:size(ubdtrik,1)
%    if ubdtrik(i,3)==botz
%       % newubi=newubi;
%    else
%        newubi=cat(1,newubi,ubdtrik(i,:));
%    end
% end
% chkarray = isempty(newubi)
% if chkarray == 1
%     
% elseif chkarray == 0
% 
% %% old
% zmax=0;
% newunbi=newubi; %%% finding error
% newxbi=[];
% ac=newubi(:,3);
% av=unique(ac,'rows','stable');
% as=sort(av,'ascend');
% %find the maximum of the bounding box
% % by=[0 0 1];  %% have to correct something
% % 
% for j=1:1:size(newubi,1)
%     for i=1:4:size(mg,1)
%       % if (mg(i,:)==by)%|
%        %if (mg(i,3)<=0)%|(mg(i,3)>0)
%          if (newubi(j,1) == mg(i+1,1))&(newubi(j,2) == mg(i+1,2))
%              newxbi=cat(1,newxbi,mg((i+1:i+3),:)); 
%          else if (newubi(j,1) == mg(i+2,1))&(newubi(j,2) == mg(i+2,2))
%              newxbi=cat(1,newxbi,mg((i+1:i+3),:)); 
%              else if (newubi(j,1) == mg(i+3,1))&(newubi(j,2) == mg(i+3,2))
%                 newxbi=cat(1,newxbi,mg((i+1:i+3),:));     
%                  end
%             %end
%          end
%          end
%        if i>size(mg,1)
%             break
%     end  
%     end
% end
% 
% %for i=1:1:size(newxbi,1)
%     zmax= max(newxbi(:,3));
% %end
% 
% % rt=min(as);
% % xyt=max(as);                    %1008;  %21.5
% % rad = xyt-rt;
% % yt = rt+(2*rad);
% %yt=zmax;
% angn=[];
% yt=zmax;
%  for k=1:1:size(as,1)
%       rt= as(k,1);
% for j=1:1:size(newubi,1)
%     for i=1:4:size(mg,1)
%         if (mg(i+1,3) >= rt)&(mg(i+2,3) >= rt)&(mg(i+3,3) >= rt)
%             if yt == topz
%             if (mg(i+1,3) < yt)&(mg(i+2,3) < yt)&(mg(i+3,3) < yt)
%        if (newubi(j,1) == mg(i+1,1))&(newubi(j,2) == mg(i+1,2))  
%            newunbi=cat(1,newunbi,mg((i+1:i+3),:));
%            angn=cat(1,angn,mg(i,:));
%        elseif (newubi(j,1) == mg(i+2,1))&(newubi(j,2) == mg(i+2,2))
%           newunbi=cat(1,newunbi,mg((i+1:i+3),:));
%            angn=cat(1,angn,mg(i,:));
%        elseif (newubi(j,1) == mg(i+3,1))&(newubi(j,2) == mg(i+3,2))
%             newunbi=cat(1,newunbi,mg((i+1:i+3),:));
%                   angn=cat(1,angn,mg(i,:));
%        end
% %         end
% %        end
%        end
%             else
%           if (mg(i+1,3) <= yt)&(mg(i+2,3) <= yt)&(mg(i+3,3) <= yt)
%        if (newubi(j,1) == mg(i+1,1))&(newubi(j,2) == mg(i+1,2))  
%            newunbi=cat(1,newunbi,mg((i+1:i+3),:));
%            angn=cat(1,angn,mg(i,:));
%        elseif (newubi(j,1) == mg(i+2,1))&(newubi(j,2) == mg(i+2,2))
%            newunbi=cat(1,newunbi,mg((i+1:i+3),:));
%            angn=cat(1,angn,mg(i,:));
%        elseif (newubi(j,1) == mg(i+3,1))&(newubi(j,2) == mg(i+3,2))
%             newunbi=cat(1,newunbi,mg((i+1:i+3),:));
%                  angn=cat(1,angn,mg(i,:));
%        end
%        end
%       end
%         end
%    if i>size(mg,1)
%             break
%     end 
%     end
%     if j>size(newubi,1)
%             break
%     end  
%   end 
%     
%  end
% kal=[];
% newkal=[];
% for i=1:3:size(newunbi,1)
%     kal=cat(2,newunbi(i,:),newunbi(i+1,:),newunbi(i+2,:));
%     newkal=cat(1,newkal,kal(1,:));
%     if i>size(newunbi,1)
%             break
%     end  
% end
% %end
% newuqi=unique(newkal,'rows','stable');
% gal=[];
% for i=1:1:size(newuqi,1)
%     gal=cat(1,gal,newuqi(i,(1:3)),newuqi(i,(4:6)),newuqi(i,(7:9)));
%         if i>size(newuqi,1)
%             break
%     end 
% end
% % 
% nal=[];
% 
% % Acommon = intersect(mg,gal)
% % Arr3 = setxor(mg,Acommon)
% % %% new
% % for i=1:1:size(ubdtrik,1)
% % newubi=cat(1,newubi,ubdtrik(i,:));
% % end
% 
% %% Buildable geometries
% p=0;
% q=0;
% s=0;
% colonenormal=[];
% 
% for i=1:4:size(mg,1)
%     nal= mg((i+1:i+3),:);
%     colonenormal = mg(i,:);
%     p=p+1;
%     for j=1:3:size(gal,1)
%         if ((nal(1,1)-gal(j,1))==0) & ((nal(2,1)-gal(j+1,1))==0) & ((nal(3,1)-gal(j+2,1))==0)&((nal(1,2)-gal(j,2))==0) & ((nal(2,2)-gal(j+1,2))==0) & ((nal(3,2)-gal(j+2,2))==0)&((nal(1,3)-gal(j,3))==0) & ((nal(2,3)-gal(j+1,3))==0) & ((nal(3,3)-gal(j+2,3))==0)
%                     %newbdi=newbdi;
%                    % j=j+3;
%                   q=p;
%                 else
%                   s=p;
%                     %j=j+3;
%                    
%         end
%         if j>size(gal,1)
%             break
%     end 
% %                if s>p
% %                        break
% %                    end
%     end
%      if s>q
%         newbodi=cat(1,newbodi,nal((1:3),:));
%         colnormal= cat(1,colnormal,colonenormal(1,:));
%        else 
%            newbodi=newbodi;
%        end
% %        q=p;
%        if i>size(mg,1)
%             break
%     end 
% end
% 
% %% next       
% % j=1;
% %  %sx=150; %%%problem
% %  for k =1:1:15
% %      ubdtrik=cat(1,ubdtrik,bdi((j:j+2),:));
% %      j=j+3;
% %      if k> 15 %4
% %             break
% %      end  
% %  end
% 
% %% plotting figures
% figure
%     view(3); camlight; %axis vis3d
%     axis off
%     axis equal
%     qp=1;
%  for i=1:size(newubi,1)
%         vc=newubi(qp:qp+2,:);
%         f=[1 2 3];
%      tic; p=patch('Faces',f,'Vertices',vc); toc;
%         set(p,'facecolor','m','facealpha',0.5);
%         set(p,'Edgecolor','k','linewidth',2);
%         hold on
%         qp=qp+3;
%         i=i+1;
%         if qp> size(newubi,1)
%             break 
%         end
%  end
%  
% %  figure
% %     view(3); camlight; %axis vis3d
% %     axis off
% %     axis equal
% %     qp=1;
% %  for i=1:size(newunbi,1)
% %         vc=newunbi(qp:qp+2,:);
% %         f=[1 2 3];
% %      tic; p=patch('Faces',f,'Vertices',vc); toc;
% %         set(p,'facecolor','m','facealpha',0.5);
% %         set(p,'Edgecolor','k','linewidth',2);
% %         hold on
% %         qp=qp+3;
% %         i=i+1;
% %         if qp> size(newunbi,1)
% %             break 
% %         end
% %  end
%  
% figure
%     view(3);  %axis vis3d
%     axis off
%     axis equal
%     qp=1;
%  for i=1:size(gal,1)
%         vd=gal(qp:qp+2,:);
%         f=[1 2 3];
%         p=patch('Faces',f,'Vertices',vd);
%         set(p,'facecolor','m','facealpha',0.5);
%         set(p,'Edgecolor','k','linewidth',2);
%         hold on
%         qp=qp+3;
%         i=i+1;
%         if qp> size(gal,1)
%             break 
%         end
%  end
% 
%  figure
%     view(3);  %axis vis3d
%     axis off
%     axis equal
%     qp=1;
%  for i=1:size(newbodi,1)
%         ve=newbodi(qp:qp+2,:);
%         f=[1 2 3];
%         p=patch('Faces',f,'Vertices',ve);
%         set(p,'facecolor','m','facealpha',0.5);
%         set(p,'Edgecolor','k','linewidth',2);
%         hold on
%         qp=qp+3;
%         i=i+1;
%         if qp> size(newbodi,1)
%             break 
%         end
%  end
% end
%  
% 
% %% tips
% 
% % If no overhang detected, combine bdi and ubdtrik for figure base plot
% %  figure
% %     view(3); camlight; axis vis3d
% %     qp=1;
% %  for i=1:size(newxbi,1)
% %         vf=newxbi(qp:qp+2,:);
% %         f=[1 2 3];
% %         p=patch('Faces',f,'Vertices',vf);
% %         set(p,'facecolor','b','facealpha',0.5);
% %         set(p,'Edgecolor','k','linewidth',2);
% %         hold on
% %         qp=qp+3;
% %         i=i+1;
% %         if qp> size(newxbi,1)
% %             break 
% %         end
% %  end
%  
% %end
% %  
% % x1 = [0 1 2];
% % z1 = [0 0 3];
% % x2 = [2 3 4];
% % z2=[0 0 3];
% % polyin = polyshape({x1,x2},{z1,z2});
% % [xlim,ylim] = boundingbox(polyin);
% % plot(polyin)
% % hold on
% % plot(xlim,ylim,'r*',xlim,fliplr(ylim),'r*')
% 
% % figure
% %     axis equal
% %     pq=1;
% %  for i=1:size(ubdrik,1)
% %         vb=ubdrik(pq:pq+2,:);
% %         f=[1 2 3];
% %         p=patch('Faces',f,'Vertices',vb);
% %         set(p,'facecolor','b','facealpha',0.5);
% %         set(p,'Edgecolor','k','linewidth',2);
% %         hold on
% %         pq=pq+3;
% %         i=i+1;
% %         if pq>(size(ubdrik,1))
% %             break 
% %         end
% %  end
% 
%  %  if ((newubi(j,3)) == rt)
% %        if (newubi(j,1) == vvw(i,1))&(newubi(j,2) == vvw(i,2))%&(newubi(j,3) == vvw(i,3)))
% %            if (newubi(j,1) == vvw(i+1,1))&(newubi(j,2) == vvw(i+1,2))%&(newubi(j,3) == vvw(i+1,3))
% %                     newunbi=cat(1,newunbi,vvw((i:i+2),:));
% %            end
% %        end
% %        
% %        if (newubi(j,1) == vvw(i+1,1))&(newubi(j,2) == vvw(i+1,2))%&(newubi(j,3) == vvw(i,3)))
% %            if (newubi(j,1) == vvw(i+2,1))&(newubi(j,2) == vvw(i+2,2))%&(newubi(j,3) == vvw(i+1,3))
% %                     newunbi=cat(1,newunbi,vvw((i:i+2),:));
% %            end
% %        end
% %        if (newubi(j,1) == vvw(i+2,1))&(newubi(j,2) == vvw(i+2,2))%&(newubi(j,3) == vvw(i,3)))
% %            if (newubi(j,1) == vvw(i,1))&(newubi(j,2) == vvw(i,2))%&(newubi(j,3) == vvw(i+1,3))
% %                     newunbi=cat(1,newunbi,vvw((i:i+2),:));
% %            end
% %        end