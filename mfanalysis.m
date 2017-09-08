function [ O ] = mfanalysis( traj, scale, q, mfdfagr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
shifts=traj;

    
[Hq1,tq,hq,Dq,Fq] = MFDFA1(shifts,scale,q,1,0);
correction=Fq(1,1:10)/sum(Fq(1,1:10));
shiftperm=shifts(randperm(length(shifts)));
[Hq2,tq,hq,Dq,Fq2] = MFDFA1(shiftperm,scale,q,1,0);
shiftperm1=shifts(randperm(length(shifts)));
[Hq2,tq,hq,Dq,Fq3] = MFDFA1(shiftperm1,scale,q,1,0);
shiftperm2=shifts(randperm(length(shifts)));
[Hq2,tq,hq,Dq,Fq4] = MFDFA1(shiftperm2,scale,q,1,0);


Faver=(Fq2+Fq3+Fq4)/3;
correction2=Faver(1,1:10)/sum(Faver(1,1:10));




j=1;
for i=1:9
    if (correction(i)/correction(i+1)<0.01)||(correction2(i)/correction2(i+1)<0.01)
    j=i;
    end;
end;

scmin=scale(j+1);
scmax=scale(10);
                
                %scmax=ceil(mean(aviter));
scres=10;
exponents=linspace(log2(scmin),log2(scmax),scres);
scale=round(2.^exponents);    

[Hq1,tq,hq,Dq,Fq] = MFDFA1(shifts,scale,q,1,0 );

[Hq2,tq,hq,Dq,Fq2] = MFDFA1(shiftperm,scale,q,1,0);

[Hq3,tq,hq,Dq,Fq3] = MFDFA1(shiftperm1,scale,q,1,0);

[Hq4,tq,hq,Dq,Fq4] = MFDFA1(shiftperm2,scale,q,1,0);



figure;
plot(Hq1-(Hq2+Hq3+Hq4)/3);
drawnow;

O=Hq1-(Hq2+Hq3+Hq4)/3;
end

