% This code draws the evolutionary trajectories of juvenile and adult
% resistance (and also shows how disease prevalence and total population
% density change with time). 

%% Define Parameters:

version=5;
% version = 1 : Juv res with maturation, adult res with reproduction
% version = 2 : Juv res with juv mortality, adult res with adult mortality
% version = 3 : Juv res with reproduction, adult res with adult mortality
% version = 4 : Juv res with maturation, adult res with adult mortality
% version = 5 : Juv res with juv mortality, adult res with  reproduction
% version = 6 : Juv res with reproduction, adult res with reproduction
if version~=1 && version~=2 && version~=3 && version~=4 && version~=5 && version~=6
    disp('Please choose one version of the model')
    return
end

% Parameter values:
t_max=100;
a0=5;
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
beta0=8;
resJmin=0;
resJmax=1;
resJstart=0.5;
resAmin=0;
resAmax=1;
resAstart=0.5;
h=1;
f=0.5;
alpha=0; 
res0=101;
nevol=5000;
bJ=1;
bA=1;

% Set up vectors:
RESJ=[];
RESA=[];
DISPREV=[];
N=[];

%% Set up Initial Conditions for a Monomorphic Population

% There is one strain of host initially:
strain_totalJ = 1;
strain_totalA = 1;

% Initial values of trade-off parameters:
if version==1
    astart=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)));
    gstart=g0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
    bJstart=1;
    bAstart=1;
elseif version==2
    bJstart=1+c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g));
    bAstart=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
    astart=a0;
    gstart=g0;
elseif version==3
    bAstart=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
    astart=a0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
    gstart=g0;
    bJstart=1;
elseif version==4
    gstart=g0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
    bAstart=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
    astart=a0;
    bJstart=1;
elseif version==5
    bJstart=1+c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g));
    astart=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)));
    bAstart=1;
    gstart=g0;
elseif version==6
    astart=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
    gstart=g0;
    bJstart=1;
    bAstart=1;
end

% Initial population distribution (SJ, SA, IJ, IA):
init_pop = [0.1,0.1,0.1,0.1];

% Initial conditions
ResJ = linspace(resJmin,resJmax,res0);
ResA = linspace(resAmin,resAmax,res0);
initialJ = find(ResJ>=resJstart,1);
initialA = find(ResA>=resAstart,1);
resJ_start = ResJ(initialJ);
resA_start = ResA(initialA);
indexJ_start = initialJ;
indexA_start = initialA;

%% Allow both traits to evolve
% This section can be run multiple times to extend the number of
% evolutionary timesteps without having to start the simulation again. 

% The juvenile and adult resistance traits have equal mutation rates:
resJmutationprob=0.5;
[resJ_start,resA_start,init_pop,strain_totalJ,strain_totalA,indexJ_start,indexA_start,RESJnew,RESAnew,DISPREVnew,Nnew] = simulation_function(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJmin,resJmax,resJ_start,resAmin,resAmax,resA_start,h,f,resJmutationprob,init_pop,strain_totalJ,strain_totalA,indexJ_start,indexA_start,res0,nevol,version);
RESJ=[RESJ;RESJnew];
RESA=[RESA;RESAnew];
DISPREV=[DISPREV;DISPREVnew];
N=[N;Nnew];

%% Make the Plot

% Plot the juvenile resistance trajectory:
clf
aa0=1;
RESJ0=log10(RESJ);
RESJ0(RESJ0<-aa0)=-aa0;
RESJ0=(RESJ0+aa0)/aa0;
subplot(1,3,1)
imagesc(RESJ0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evo time')
xlabel('resJ')
set(gca,'xtick',1:1:res0,'xticklabel',round(ResJ(1:1:res0)*10000)/10000);
set(gca,'xtick',1:(res0-1)/4:res0,'xticklabel',round(ResJ(1:(res0-1)/4:res0)*10000)/10000);
title('(a)')

% Plot the adult resistance trajectory:
ab0=1;
RESA0=log10(RESA);
RESA0(RESA0<-ab0)=-ab0;
RESA0=(RESA0+ab0)/ab0;
subplot(1,3,2)
imagesc(RESA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evo time')
xlabel('resA')
set(gca,'xtick',1:1:res0,'xticklabel',round(ResA(1:1:res0)*10000)/10000);
set(gca,'xtick',1:(res0-1)/4:res0,'xticklabel',round(ResJ(1:(res0-1)/4:res0)*10000)/10000);
title('(b)')

% Plot the disease prevalence and (scaled) total population density:
subplot(1,3,3)
h1=plot(DISPREV);
hold on
h2=plot(N);
xlabel('Evo Time')
legend([h1,h2],{'Disease prevalence','Total population density'})
ylim([0,1])
title('(c)')
