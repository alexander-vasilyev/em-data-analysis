addpath('amd_bin');
addpath('amd_mon');
addpath('normal_easy');
addpath('normal_hard');
% % h=load('6.txt');
% h0=load('axreM.txt');
% h1=load('chstM.txt');
% h2=load('dileM.txt');
% h3=load('ebmeM.txt');
% h4=load('frklM.txt');
% h5=load('hemaM.txt');
% h6=load('hemeM.txt');
% h7=load('hescM.txt');
% h8=load('hewaM.txt');
% h9=load('inkrM.txt');
% h10=load('maneM.txt');
% h11=load('mastM.txt');
% h12=load('mikuM.txt');
% h13=load('reriM.txt');
% h14=load('rohaM.txt');
% h15=load('ulgeM.txt');
% hcell={h0,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15};
%hcell={h0,h2,h4,h8,h10,h11,h13,h14,h15};
% h1=load('chst.txt');
% h2=load('dile.txt');
% h3=load('ebme.txt');
% h4=load('frkl.txt');
% h5=load('hema.txt');
% h6=load('inkr.txt');
% h7=load('mane.txt');
% h8=load('mast.txt');
% h9=load('miku.txt');
% h10=load('reri.txt');
% h11=load('roha.txt');
% h12=load('ulge.txt');
% hcell={h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12};
% 
h1=load('alsaCC.txt');
h2=load('bescCC.txt');
h3=load('edbeCC.txt');
h4=load('erbuCC.txt');
h5=load('gegrCC.txt');
h6=load('heweCC.txt');
h7=load('ilpuCC.txt');
h8=load('jolaCC.txt');
h9=load('mahrCC.txt');
h10=load('maraCC.txt');
h11=load('paszCC.txt');
h12=load('roraCC.txt');
h13=load('sidoCC.txt');
h14=load('vokuCC.txt');
h15=load('wotsCC.txt');
%hcell={h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15};
hcell={h1,h2,h3,h4,h6,h8,h9,h11,h12};

bord=-200;


trajcl={};
timetrajcl={};
avlengh=repmat(1,1,length(hcell));
dim=64;
for pen=1:length(hcell)
h=hcell{pen};
sz=size(h);
sxtraj=[];
sytraj=[];
sttraj=[];
sltraj=[];
i=1;
ixtraj=[];
iytraj=[];
tttraj=[];
lttraj=[];
while i<sz(1)
    if h(i,4)>0 
        ixtraj(end+1)=h(i,4);
        iytraj(end+1)=h(i,5);
        tttraj(end+1)=h(i,3);        
    end
%     if h(i,4)>0 && h(i+1,4)>0
%         sttraj(end+1)=h(i+1,1)-h(i,2);    
%     end
%     if h(i,4)==0 
%         sttraj(end+1)=0;
%     end
    if h(i,4)==0 && h(i+1,4)>0    
        ixtraj=ixtraj-250+bord;
        iytraj=iytraj+bord;
        msx=max(ixtraj);
        msy=max(iytraj);
        mtx=min(ixtraj);
        mty=min(iytraj);
        if ((mty>0)&&(msy<1000+2*bord))&&((mtx>0)&&(msx<1150+2*bord))
%             if length(ixtraj)>5 && length(ixtraj)<15
%                 sxtraj=[sxtraj,0,ixtraj(1:5)];
%                 sytraj=[sytraj,0,iytraj(1:5)];
%             elseif length(ixtraj)>15
%                 sxtraj=[sxtraj,0,ixtraj(1:15)];
%                 sytraj=[sytraj,0,iytraj(1:15)];
%             end;
%                 sxtraj=[sxtraj,0,ixtraj];
%                 sytraj=[sytraj,0,iytraj];
%                 sttraj=[sttraj,0,tttraj];
             if length(ixtraj)>8
                sxtraj=[sxtraj,0,ixtraj(1:7)];
                sytraj=[sytraj,0,iytraj(1:7)];
                sttraj=[sttraj,0,tttraj(1:7)];
              end;
        end;
        
        ixtraj=[];
        iytraj=[];
        tttraj=[];
    end;
    i=i+1;
end;
% sytraj(1)=[];
% min(sytraj)
% sxtraj(1)=[];
% max(sytraj)
difxtraj=(sxtraj(1:length(sxtraj)-1)-sxtraj(2:length(sxtraj)));
difytraj=(sytraj(1:length(sytraj)-1)-sytraj(2:length(sytraj)));
sltraj=sqrt((difxtraj.*difxtraj)+(difytraj.*difytraj));

% x=sxtraj;
% y=sytraj;
% 
% nbins=[70 70];
% his = histogram2(x,y,nbins);
% min_x = min(x); max_x = max(x); step_x = (max_x - min_x)/nbins(1);
% min_y = min(y); max_y = max(y); step_y = (max_y - min_y)/nbins(2);
% 
% % make grid
% surf_z = his.Values;
% surf_x = [min_x + step_x/2 : step_x : max_x - step_x/2];
% surf_y = [min_y + step_y/2 : step_y : max_y - step_y/2];
% [xx, yy] = meshgrid(surf_x, surf_y);
% 
% % plot 3D surface
% figure
% surf(xx',yy',surf_z);
% hold on;
% xlabel('X axis');
% ylabel('Y axis');


for i=1:length(sxtraj)
    if (sxtraj(i)~=0)
    sxtraj(i)=ceil(sxtraj(i)/((1150+2*bord)/dim));
    sytraj(i)=ceil(sytraj(i)/((1000+2*bord)/dim));
        if sxtraj(i)==0
            sxtraj(i)=1;
        end;
        if sytraj(i)==0
            sytraj(i)=1;
        end;
    end;
end;
traj=[sxtraj;sytraj];
trajcl{end+1}=traj;
timetrajcl{end+1}=sttraj;
end;

ltrj=length(trajcl);
for i=1:ltrj
    k=ltrj+1-i;
    if length(trajcl{k})<500
        trajcl(k)=[];
    end;
end;