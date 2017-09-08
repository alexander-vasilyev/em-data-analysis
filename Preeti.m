addpath('C:/Users/av306/Desktop/coding/mf/multifractal paper/mf-dfa-infomax/code/data')
addpath('C:/Users/av306/Desktop/coding/mf/multifractal paper/mf-dfa-infomax/code/MFDFA')
addpath('C:/Users/av306/Desktop/coding/mf/multifractal paper/mf-dfa-infomax/code/MDP')
addpath('C:/Users/av306/Desktop/coding/mf/multifractal paper/mf-dfa-infomax/code/preeti')
% h=load('6.txt');
h1=load('2.txt');
h2=load('3.txt');
h3=load('4.txt');
h4=load('5.txt');
h5=load('6.txt');
h6=load('7.txt');
h7=load('8.txt');
h8=load('9.txt');
h9=load('10.txt');
h10=load('11.txt');
h11=load('12.txt');
hcell={h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11};
traj=repmat(1,2,3000);
trajcl={traj,traj,traj,traj,traj,traj,traj,traj,traj,traj,traj};
avlengh=repmat(1,1,11);
avpen=repmat(1,1,11);
for pen=1:11
h=hcell{pen};    
amplitude=0;
distangle=0;
dirangles=0;
shiftsx=zeros(1,1);
shiftsy=zeros(1,1);
shiftst=zeros(1,1);
shiftsa=zeros(1,1);
angles=zeros(1,1);
time=0;
lhist= repmat( zeros, [1 181]);
Xtraj = 0;
Ytraj = 0;
anghist= repmat( zeros, [1 361]);
lhist=repmat(zeros,[1 361]);
sz=size(h);

xarray=h(:,7);
yarray=h(:,8);
minx=min(xarray);
miny=min(yarray);
sxtraj=0;
sytraj=0;
for i=1:(sz(1)-1)
   
    if (h(i+1,2)==h(i,2)+1)
        %if h(i,7)>0 && h(i,8)>0
%         Xtraj(end+1)=h(i,7)-minx;
%         Ytraj(end+1)=h(i,8)-miny;
       % else 
%         Xtraj(end+1)=1;
%         Ytraj(end+1)=1;
%         end;
        saclength=sqrt((h(i+1,7)-h(i,7))*(h(i+1,7)-h(i,7))+(h(i+1,8)-h(i,8))*(h(i+1,8)-h(i,8)));
        %if saclength<1000 && saclength>20
        shiftsx(end+1)=h(i+1,7)-h(i,7);
        shiftsy(end+1)=-h(i+1,8)+h(i,8);
      
        %end;
        if  (saclength>20)
%         Xtraj(end+1)=h(i,7)-minx+1;
%         Ytraj(end+1)=h(i,8)-miny+1;
        sxtraj(end+1)=h(i,7);
        sytraj(end+1)=h(i,8);
        
        end;
        
        if (saclength>20) && (saclength<150)
        shiftst(end+1)=h(i+1,3)-h(i,3);
        shiftsa(end+1)=saclength;
        end;
        sacangle=atan2((h(i+1,8)-h(i,8)),(h(i+1,7)-h(i,7)));
        if sacangle<0            
            sacangle=sacangle+2*pi;
        end;
        angles(end+1)=sacangle;
        anghist(ceil(sacangle*180/pi+0.01))=anghist(ceil(sacangle*180/pi+0.01))+1;
        lhist(ceil(sacangle*180/pi+0.01))=lhist(ceil(sacangle*180/pi+0.01))+h(i,5);
    else
        
        minsx=min(sxtraj(1:length(sxtraj)));
        minsy=min(sytraj(1:length(sxtraj)));
        Xtraj=[Xtraj,0,sxtraj(2:length(sxtraj))-minsx+10];
        Ytraj=[Ytraj,0,sytraj(2:length(sytraj))-minsy+10];
        sxtraj=0;
        sytraj=0;
    end;
    amplitude(end+1)=h(i,5);
    time(end+1)=h(i,4)-h(i,3);
end;
Xtraj(1)=[];
Ytraj(1)=[];
p=polyfit(shiftsa,shiftst,1);
figure
hist(angles,200);
figure
plot(lhist./anghist);
lnx=length(Xtraj);
for i=1:lnx
    if Xtraj(i)>550
        Xtraj(i)=Xtraj(i)-550;
    end;
    if Ytraj(i)>550
        Ytraj(i)=Ytraj(i)-550;
    end;
end;

x=Xtraj;
y=Ytraj;

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





mxx=max(Xtraj);
mxy=max(Ytraj);
if mxx>mxy
    crx=mxx;
else 
    crx=mxy;
end;
p1=p(1)/p(2);
avpen(pen)=p1/(128/crx);
avpen=abs(avpen);
avlengh(pen)=ceil(sqrt(mean(shiftsa.*shiftsa))/(crx/128));

%healthy settings
lnx=length(Xtraj);

for i=1:lnx
    
    if (Xtraj(i)~=0)
    
    Xtraj(i)=ceil(Xtraj(i)/(crx/128));
    Ytraj(i)=ceil(Ytraj(i)/(crx/128));
        if Xtraj(i)==0
        Xtraj(i)=1;
        end;
        if Ytraj(i)==0
        Ytraj(i)=1;
        end;
        if Xtraj(i)>128
            Xtraj(i)=[];
            Ytraj(i)=[];
        end;
        if Ytraj(i)>128
            Ytraj(i)=[];
            Xtraj(i)=[];
        end; 
    end;
end;
lnx=length(Xtraj);

% for i=1:lnx-1
%    if (i<lnx)&&(Xtraj(i)~=0) 
%        if (i<lnx)&&((Xtraj(i+1)-Xtraj(i))*(Xtraj(i+1)-Xtraj(i))+(Ytraj(i+1)-Ytraj(i))*(Ytraj(i+1)-Ytraj(i)))<14
%            Xtraj(i)=[];
%            Ytraj(i)=[];
%            lnx=lnx-1;
%        end;
%    end;
% end;
% lnx=length(Xtraj);
% for i=1:lnx-1
%    if (i<lnx)&&(Xtraj(i)~=0)   
%        if (i<lnx)&&((Xtraj(i+1)-Xtraj(i))*(Xtraj(i+1)-Xtraj(i))+(Ytraj(i+1)-Ytraj(i))*(Ytraj(i+1)-Ytraj(i)))>100
%            Xtraj(i)=[];
%            Ytraj(i)=[];
%            lnx=lnx-1;
%        end;  
%    end;
% end;

for id=1:length(Xtraj)-2
    if abs(Xtraj(id))>128
        Xtraj(id)=[];
        Ytraj(id)=[];
    end;
    if abs(Ytraj(id))>128
        Ytraj(id)=[];
        Xtraj(id)=[];
    end; 
end;

traj=[Xtraj;Ytraj];
trajcl{pen}=traj;
end;