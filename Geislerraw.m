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
h=[m;n];
%h=[p;l];
%h=[k;o];
%h=[h];

h1=h(1:length(h)-1,1);
h2=h(2:length(h),1);
dh=h2-h1;

scmin=36;
scmax=400;
scres=32;
exponents=linspace(log2(scmin),log2(scmax),scres);
scale=round(2.^exponents);
q=linspace(-10,10,21);
xtrn=randn(length(dh),1);
shifts1=abs(dh+0.0*xtrn);
[Hq1,tq,hq,Dq,Fq] = MFDFA1(shifts1,scale,q,1,1);
%HqcorX=mfanalysis(shifts1,scale,q,mfdfagr);