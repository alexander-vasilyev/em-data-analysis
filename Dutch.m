addpath('C:/Users/av306/Desktop/coding/mf/multifractal paper/mf-dfa-infomax/code/data')
addpath('dutch');
addpath('p1');
addpath('p2');
addpath('preeti');
% h=load('cont1.txt');
% k=load('cont2.txt');
% l=load('cont3.txt');
% m=load('cont4.txt');
% n=load('cont8.txt');
% o=load('cont9.txt');
% p=load('cont10.txt');
% r=load('cont5.txt');
% 
% 
% hcell=[h;k;l;m;n;o;p;r];
% hcell={hcell};
%  h=load('md1.txt');
%  k=load('md2.txt');
%  l=load('md3.txt');
%  n=load('md4.txt');
%   hcell={h,k,l,n};
% h=[n];
o=load('sim1.txt');
k=load('sim2.txt');
l=load('sim3.txt');
n=load('sim4.txt');
m=load('sim5.txt');
hcell={o,k,l,n,m};
lhist= repmat( zeros, [1 181]);
anghist= repmat( zeros, [1 181]);
lenghist= repmat(zeros, [1 181]);
lenghist1=repmat(zeros, [1 181]);
dirhist= repmat(zeros, [1 181]);
dirhist1= repmat(zeros, [1 181]);
dirhist2=repmat(zeros, [1 181]);

amplitude=0;
distangle=0;
dirangles=0;
shifts=zeros(1,1);
shifts1=zeros(1,1);
shifts2=zeros(1,1);
shifts3=zeros(1,1);
locs=zeros(1,1);
locs1=zeros(1,1);
sz=size(h);
avlength=repmat(1,1,length(hcell));

traj=repmat(1,2,3000);
trajcl={};
dim=64;
ltrajectory=[];
ttrajectory=[];
for pen= 1:length(hcell)
h=hcell{pen};    
sz=size(h);
shifts=zeros(1,1);
shifts1=zeros(1,1);
shifts2=zeros(1,1);
shifts3=zeros(1,1);
Xtraj = 0;
Ytraj = 0; 
    
for i= 1:sz(1)
    
    j=0;
    if (h(i,1)>4)
        while ((j<h(i,1)-1)&&(1+5*j<sz(2)))
            j=j+1;
            if (h(i,5*j-1)>0)&& (h(i,5*j)>0)
            Xtraj(end+1)=h(i,5*j-1);
            Ytraj(end+1)=h(i,5*j);
            end;
        end;
        Xtraj(end+1)=0;
        Ytraj(end+1)=0;
    end;
    
    while ((j<h(i,1)-1)&&(1+5*j<sz(2)))
        j=j+1;
        if(h(i,1+5*j)<sz(2));
            amplitude(end+1)=  h(i,1+5*j);
        end;
    end;
    j=0;
    jk=0;
    while((j<h(i,1)-1)&&(5*(j+2)<sz(2)))
        j=j+1;
        
        %if (abs(h(i,5*j-1)-h(i,5*(j+1)-1))<1000)&&(abs(h(i,5*j)-h(i,5*(j+1)))<1000)
        shifts(end+1)=h(i,5*j)-h(i,5*(j+1));
        shifts1(end+1)=-h(i,5*j-1)+h(i,5*(j+1)-1);
        slength=sqrt((h(i,5*j)-h(i,5*(j+1)))*(h(i,5*j)-h(i,5*(j+1)))+(h(i,5*j-1)-h(i,5*(j+1)-1))*(h(i,5*j-1)-h(i,5*(j+1)-1)));
        
        if (abs((h(i,5*j-2)-h(i,5*(j+1)-3)))<2000)&& abs(slength)<1000
            shifts2(end+1)=-(h(i,5*j-2)-h(i,5*(j+1)-3)); 
            shifts3(end+1)=slength;
        end;
        jk=jk+1;
        finish=h(i,5*j);
        finish1=h(i,5*j-1);
        %end;
       
        angle= acosd((h(i,5*j-1)-h(i,5*(j+1)-1))/sqrt((h(i,5*j-1)-h(i,5*(j+1)-1))*(h(i,5*j-1)-h(i,5*(j+1)-1))+(h(i,5*j)-h(i,5*(j+1)))*(h(i,5*j)-h(i,5*(j+1)))));
%         angle= atan2d((h(i,5*j)-h(i,5*(j+1))),(h(i,5*j-1)-h(i,5*(j+1)-1)));
        if (~isnan(angle))
            distangle(end+1)= angle;
        
            if(sqrt((h(i,5*j-1)-h(i,5*(j+1)-1))*(h(i,5*j-1)-h(i,5*(j+1)-1))+(h(i,5*j)-h(i,5*(j+1)))*(h(i,5*j)-h(i,5*(j+1))))<1000)
            dirhist(floor(real(angle)+1))=dirhist(floor(real(angle)+1))+1;
            dirhist1(floor(real(angle)+1))=dirhist1(floor(real(angle)+1))+sqrt((h(i,5*j-1)-h(i,5*(j+1)-1))*(h(i,5*j-1)-h(i,5*(j+1)-1))+(h(i,5*j)-h(i,5*(j+1)))*(h(i,5*j)-h(i,5*(j+1))));
            dirhist2(floor(real(angle)+1))=dirhist2(floor(real(angle)+1))+sqrt((h(i,5*j-1)-h(i,5*(j+1)-1))*(h(i,5*j-1)-h(i,5*(j+1)-1))+(h(i,5*j)-h(i,5*(j+1)))*(h(i,5*j)-h(i,5*(j+1))))/(h(i,5*j-2)-h(i,5*j-3));
           
            end;
        end;
    end;
    for il=1:jk
        locs(end+1)=1;
        locs1(end+1)=1;
    end;
    for ko=1:8
        locs(end+1)=0;
        locs1(end+1)=0;
%         shifts(end+1)=0;
%         shifts1(end+1)=0;
%         shifts2(end+1)=0;
        
    end;
%     while((j<h(i,1)-2)&&(5*(j+3)<sz(2)))
%         j=j+1;
%         difx1= h(i,5*j-1)-h(i,5*(j+1)-1);
%         difx2= h(i,5*(j+1)-1)-h(i,5*(j+2)-1);
%         dify1= h(i,5*j)-h(i,5*(j+1));
%         dify2= h(i,5*(j+1))-h(i,5*(j+2));
%         
%         sacangle= acosd((difx1*difx2+dify1*dify2)/(sqrt(difx1*difx1+dify1*dify1)*sqrt(difx2*difx2+dify2*dify2)));
%         if (~isnan(sacangle))
%             dirangles(end+1)=sacangle;  
%         
%             lenghist(floor(real(sacangle)+1))=lenghist(floor(real(sacangle)+1))+1;
%             lenghist1(floor(real(sacangle)+1))=lenghist1(floor(real(sacangle)+1))+sqrt((h(i,5*j-1)-h(i,5*(j+1)-1))*(h(i,5*j-1)-h(i,5*(j+1)-1))+(h(i,5*j)-h(i,5*(j+1)))*(h(i,5*j)-h(i,5*(j+1))));  
%         end;
%     end;
end;
ltrajectory=[ltrajectory shifts3];
ttrajectory=[ttrajectory shifts2];
mix=min(Xtraj);
miy=min(Ytraj);
% Xtraj=Xtraj-mix;
% Ytraj=Ytraj-miy;
mxx=max(Xtraj);
mxy=max(Ytraj);

%healthy settings
for i=1:length(Xtraj)
    if (Xtraj(i)~=0)
    Xtraj(i)=ceil(Xtraj(i)/(1200/dim));
    Ytraj(i)=ceil(Ytraj(i)/(1000/dim));
        if Xtraj(i)==0
            Xtraj(i)=1;
        end;
        if Ytraj(i)==0
            Ytraj(i)=1;
        end;
    end;
end;
traj=[Xtraj;Ytraj];
trajcl{end+1}=traj;
avlength(pen)=ceil(mean(shifts3)/(1100/dim));


scmin=5;
scmax=16;
%shiftx=[shifts;shifts1;shifts2];
%lox=[locs;locs1];
%amd settings
% scmin=8;
% scmax=36;
% scres=64;
% exponents=linspace(log2(scmin),log2(scmax),scres);
% scale=round(2.^exponents);
% q=linspace(-10,10,21);
% 
% [Hq,tq,hq,Dq,Fq] = MFDFA1(shifts,scale,q,2,1);
% shiftperm=shifts(randperm(length(shifts)));
% [Hq,tq,hq,Dq,Fq] = MFDFA1(shiftperm,scale,q,2,1);
             
x=Xtraj;
y=Ytraj;
nbins=[70 70];
hi = histogram2(x,y,nbins);
min_x = min(x); max_x = max(x); step_x = (max_x - min_x)/nbins(1);
min_y = min(y); max_y = max(y); step_y = (max_y - min_y)/nbins(2);

% make grid
surf_z = hi.Values;
surf_x = [min_x + step_x/2 : step_x : max_x - step_x/2];
surf_y = [min_y + step_y/2 : step_y : max_y - step_y/2];
[xx, yy] = meshgrid(surf_x, surf_y);
% 
% % plot 3D surface
figure
surf(xx',yy',surf_z);
% figure;
% hist(amplitude,400);
% title('amplitude distribution');
% figure;
% hist(distangle,50);
% title('horizontal angle distribution');
% 
% figure;
% hist(dirangles,30);
% title('directional angle distribution');
% figure;
% plot(smooth(lenghist1./lenghist,'moving'));
% title('length to directional angle');
% 
% figure;
% hist(dirangles,30);
% title('average speed');
% figure;
% plot(smooth(dirhist2./dirhist,'moving'));
% title('length to directional angle');
end;
disp('preceeding saccade');
[R P]=corrcoef(ltrajectory, ttrajectory)
disp('preceeding fixation');
[R P]=corrcoef(ltrajectory(2:length(ltrajectory)), ttrajectory(1:length(ltrajectory)-1))
