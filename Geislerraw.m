addpath(genpath('p1'));
addpath(genpath('p2'));
addpath('MFDFA');
%enable graphics for md-dfa
mfdfagr=1;


%load data files
h=load('me.raw.025.txt');
k=load('me.raw.02.txt');
l=load('me.raw.015.txt');
m=load('me.raw.01.txt');
n=load('raw.mif.01.txt');
o=load('raw.mif.02.txt');
p=load('raw.mif.015.txt');
q=load('me01.new.txt');
r=load('me0.15new.txt');
s=load('me0.2new.txt');
h=[h];
h=[o;k;s];
h=[p;l;r];
h=[q;n;m];

h1=h(1:length(h)-1,1);
h2=h(2:length(h),1);
dh=h2-h1;

scmin=64;
scmax=1500;

 
%       scmin=scale(5);
%        scmax=scale(11);
scres=24;
exponents=linspace(log2(scmin),log2(scmax),scres);
scale=round(2.^exponents);
q=linspace(-10,10,21);
xtrn=randn(length(dh),1);
shifts=abs(dh+0.0*xtrn);
[Hq1,tq,hq,Dq,Fq] = MFDFA1(shifts,scale,q,1,1);
correction=log(Fq(1,1:scres));
figure
plot((correction));
shiftperm=shifts(randperm(length(shifts)));
[Hq2,tq,hq,Dq,Fq2] = MFDFA1(shiftperm,scale,q,1,mfdfagr);



%[Hqcor,Hshuf]=mfanalysis(shifts,scale,scres,q,mfdfagr);