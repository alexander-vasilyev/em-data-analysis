addpath('p1');
addpath('p2');
addpath(genpath('MDP'));
addpath(genpath('MFDFA'));
h=load('me.025.txt');
k=load('me.02.txt');
l=load('me.015.txt');
m=load('me.01.txt');
n=load('miflah.01.txt');
o=load('miflah.02.txt');
p=load('miflah.015.txt');
sz=size(n);
addz=repmat(zeros,[2 sz(1)]);
n=[addz',n];
h=[h];
 %    h=[k;o];
  %    h=[l;p];
 %     h=[m;n];

dispstats=false;


lhist= repmat( zeros, [1 181]);
anghist= repmat( zeros, [1 181]);
lenghist= repmat(zeros, [1 181]);
lenghist1=repmat(zeros, [1 181]);
dirhist= repmat(zeros, [1 181]);
dirhist1= repmat(zeros, [1 181]);
dirhist2=repmat(zeros, [1 181]);

hcell={h};

amplitude=0;
ttrajectory=0;
distangle=0;
dirangles=0;
timetrial=0;
shifts=zeros(1,1);
shifts1=zeros(1,1);
shifts2=zeros(1,1);
shifts3=zeros(1,1);
locs=zeros(1,1);
locs1=zeros(1,1);
sz=size(h);
avlengh2=repmat(1,1,4);

traj=repmat(1,2,3000);
trajcl2={traj,traj,traj,traj};


for pen= 1:1
h=hcell{pen};    
sz=size(h);
shifts=zeros(1,1);
shifts1=zeros(1,1);
shifts2=zeros(1,1);
shifts3=zeros(1,1);
Xtraj = 0;
Ytraj = 0; 
    
for i= 1:sz(1)-2
    
    
    if (h(i,12)>1)
        Xtraj(end+1)=h(i,8);
        Ytraj(end+1)=h(i,9);
    end;
    
    
    
    difx1=h(i+1,8)-h(i,8);
    dify1=h(i+1,9)-h(i,9);
    difx2=h(i+2,8)-h(i+1,8);
    dify2=h(i+2,9)-h(i+1,9);
    dift=h(i+2,5)-h(i+1,6);
    sametrial=(h(i,3)==h(i+1,3))&&(h((i+1),3)==h((i+2),3));
    sametrial1=(h(i,3)==h(i+1,3));

    difm=sqrt(difx1*difx1+dify1*dify1);
    %if (sametrial1)
      % amplitude(end+1)=  h(i,12);
      ttrajectory(end+1)=dift;
      amplitude(end+1)=21*(difm/1440);
   % end;
%     else 
%        if(h(i,3)<10000)
%        timetrial(end+1)=h(i,3);
%        end;
%     end; 
    shifts(end+1)=difx1;
    shifts1(end+1)=-dify1;
    angle= acosd(difx1/sqrt(difx1*difx1+dify1*dify1));
%         angle= atan2d((h(i,5*j)-h(i,5*(j+1))),(h(i,5*j-1)-h(i,5*(j+1)-1)));
    if (~isnan(angle))
        distangle(end+1)= angle;        
        if (difm<1000)
        dirhist(floor(real(angle)+1))=dirhist(floor(real(angle)+1))+1;
        dirhist1(floor(real(angle)+1))=dirhist1(floor(real(angle)+1))+difm;         
        end;
    end;
    sacangle= acosd((difx1*difx2+dify1*dify2)/(sqrt(difx1*difx1+dify1*dify1)*sqrt(difx2*difx2+dify2*dify2)));
    if (~isnan(sacangle)&&sametrial)&&(21*(difm/1440)>0.5)
        dirangles(end+1)=sacangle;  
         
        lenghist(floor(real(sacangle)+1))=lenghist(floor(real(sacangle)+1))+1;
        lenghist1(floor(real(sacangle)+1))=lenghist1(floor(real(sacangle)+1))+difm;  
    end;

end;
disp('preceeding saccade');
[R P]=corrcoef(ttrajectory, amplitude)
disp('preceeding fixation');
[R P]=corrcoef(ttrajectory(1:length(ttrajectory)-1), amplitude(2:length(amplitude)))
%healthy settings

avlengh2(pen)=ceil(mean(shifts3)/(1080/128));
disp(mean(amplitude));
disp(sum(lenghist(1:90))/sum(lenghist));
disp(mean(timetrial)/1000);

% shiftx=[shifts;shifts1;shifts2];
% lox=[locs;locs1];
%amd settings
 scmin=4;
 scmax=500;

% scman=scale(1)
% scmax=scale(20);

% scres=128;
% exponents=linspace(log2(scmin),log2(scmax),scres);
% scale=round(2.^exponents);
% q=linspace(-10,10,21);
% shifts=abs(shifts);
% [Hq1,tq,hq,Dq,Fq] = MFDFA1(shifts,scale,q,1,1);
% correction=log(Fq(1,1:scres));
% figure
% plot((correction));
             
% x=Xtraj;
% y=Ytraj;
% nbins=[70 70];
% h = histogram2(x,y,nbins);
% min_x = min(x); max_x = max(x); step_x = (max_x - min_x)/nbins(1);
% min_y = min(y); max_y = max(y); step_y = (max_y - min_y)/nbins(2);
% 
% % make grid
% surf_z = h.Values;
% surf_x = [min_x + step_x/2 : step_x : max_x - step_x/2];
% surf_y = [min_y + step_y/2 : step_y : max_y - step_y/2];
% [xx, yy] = meshgrid(surf_x, surf_y);
% % 
% % % plot 3D surface
% figure
% surf(xx',yy',surf_z);
if dispstats
figure;
hist(amplitude,400);
title('amplitude distribution');
figure;
hist(distangle,50);
title('horizontal angle distribution');

figure;
hist(dirangles,30);
title('directional angle distribution');
% figure;
% plot(smooth(lenghist1./lenghist,'moving'));
% title('length to directional angle');

figure;
hist(dirangles,30);
title('average speed');
% figure;
% plot(smooth(dirhist2./dirhist,'moving'));
% title('length to directional angle');
end;
end;