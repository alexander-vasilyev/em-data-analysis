addpath('dutch');
addpath('p1');
addpath('p2');
addpath('preeti');
addpath('MFDFA');
%h=load('6.txt');
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
% hcell={h1,h2,h5,h7,h8,h9,h10,h11};
% hcell={h3,h4,h6};
hcell={h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11};
%hcell={h1,h2,h7,h8,h11};
traj=repmat(1,2,3000);
trajcl={};
avlength=repmat(1,1,11);
avpen=repmat(1,1,11);
dim=36;
ltrajectory=[];
ttrajectory=[];
for pen=1:length(hcell)
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
limits=true;
for i=1:(sz(1)-1)    
    if (h(i+1,2)==h(i,2)+1)
        saclength=sqrt((h(i+1,7)-h(i,7))*(h(i+1,7)-h(i,7))+(h(i+1,8)-h(i,8))*(h(i+1,8)-h(i,8)));
 
        shiftsx(end+1)=h(i+1,7)-h(i,7);
        shiftsy(end+1)=-h(i+1,8)+h(i,8);
      
        sxtraj(end+1)=h(i,7);
        sytraj(end+1)=h(i,8);
        
        if h(i,7)>1000 || h(i,7)<0
            limits=false;            
        end;
        if h(i,8)>800 || h(i,8)<0
            limits=false;
        end;
            
        shiftst(end+1)=h(i+1,3)-h(i,4);
        shiftsa(end+1)=saclength;
        
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
        if limits==true           
                Xtraj=[Xtraj,0,sxtraj(2:length(sxtraj))];
                Ytraj=[Ytraj,0,sytraj(2:length(sytraj))];            
        end;
        sxtraj=0;
        sytraj=0;
        limits=true;
    end;
    amplitude(end+1)=h(i,5);
    time(end+1)=h(i,4)-h(i,3);
end;
Xtraj(1)=[];
Ytraj(1)=[];
p=polyfit(shiftsa,shiftst,1);
% figure
% hist(angles,200);
%  figure
%  plot(lhist./anghist);
lnx=length(Xtraj);
% for i=1:lnx
%     if Xtraj(i)>500
%         Xtraj(i)=Xtraj(i)-500;
%     end;
%     if Ytraj(i)>400
%         Ytraj(i)=Ytraj(i)-400;
%     end;
% end;
[R, P]= corrcoef(shiftsa,shiftst);



x=Xtraj;
y=Ytraj;

nbins=[70 70];
his = histogram2(x,y,nbins);
min_x = min(x); max_x = max(x); step_x = (max_x - min_x)/nbins(1);
min_y = min(y); max_y = max(y); step_y = (max_y - min_y)/nbins(2);

% make grid
surf_z = his.Values;
surf_x = [min_x + step_x/2 : step_x : max_x - step_x/2];
surf_y = [min_y + step_y/2 : step_y : max_y - step_y/2];
[xx, yy] = meshgrid(surf_x, surf_y);

% plot 3D surface
figure
surf(xx',yy',surf_z);
hold on;
xlabel('X axis');
ylabel('Y axis');

mxx=1000;
mxy=800;

% p1=p(1)/p(2);
% avpen(pen)=p1/(dim/crx);
% avpen=abs(avpen);
avlength(pen)=ceil(mean(shiftsa)/(mxx/dim));

%healthy settings
lnx=length(Xtraj);
i=1;
while i<(length(Xtraj)+1)    
    if (Xtraj(i)~=0)   
    xt=ceil(Xtraj(i)/(mxx/dim));
    yt=ceil(Ytraj(i)/(mxy/dim));
        %if (i>1)&&(abs(xt-Xtraj(i-1))+abs(yt-Ytraj(i-1))>0)
        Xtraj(i)=xt;
        Ytraj(i)=yt;
            if Xtraj(i)==0
                Xtraj(i)=1;
            end;
            if Ytraj(i)==0
                Ytraj(i)=1;
            end;
            if Xtraj(i)>dim
                Xtraj(i)=[];
                Ytraj(i)=[];           
            end;
            if Ytraj(i)>dim
                Ytraj(i)=[];
                Xtraj(i)=[];            
            end; 
        %else
        %    Xtraj(i)=[];
        %    Ytraj(i)=[];
        %end;
    end;
    i=i+1;
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
%        if ((Xtraj(i+1)-Xtraj(i))*(Xtraj(i+1)-Xtraj(i))+(Ytraj(i+1)-Ytraj(i))*(Ytraj(i+1)-Ytraj(i)))<1
%            Xtraj(i)=[];
%            Ytraj(i)=[];
%            lnx=lnx-1;           
%        end;  
%    end;
% end;

id=1;
while id<length(Xtraj)
    if abs(Xtraj(id))>dim
        Xtraj(id)=[];
        Ytraj(id)=[];
    end;
    if abs(Ytraj(id))>dim
        Ytraj(id)=[];
        Xtraj(id)=[];
    end; 
    id=id+1;
end;

traj=[Xtraj;Ytraj];
trajcl{pen}=traj;
ltrajectory=[ltrajectory,shiftsa];
ttrajectory=[ttrajectory,shiftst];
end;
[R, P]= corrcoef(ltrajectory,ttrajectory);