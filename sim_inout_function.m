function [Jout,Aout]=sim_inout_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0)

% This function determines only the endpoints of an evolutionary
% simulation.

if version~=1 && version~=2 && version~=3 && version~=4 && version~=5 && version~=6
    disp('Please choose one version of the model')
    return
end

% Set up parameters and vectors to use later:
t_max=100;
resJmin=0;
resJmax=1;
resAmin=0;
resAmax=1;
res0=51;
nevol=1000;
RESJ=[];
RESA=[];
resJmutationprob=0.5;
doneyet=0;
evotimesteps=0;

% Initial conditions:
strain_totalJ = 1;
strain_totalA = 1;
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
init_pop = [0.5*bAstart*(astart*gstart-bAstart*(bJstart+gstart))/(astart*gstart*(bJstart+gstart)),0.5*(astart*gstart-bAstart*(bJstart+gstart))/(astart*bJstart+astart*gstart),0.5*bAstart*(astart*gstart-bAstart*(bJstart+gstart))/(astart*gstart*(bJstart+gstart)),0.5*(astart*gstart-bAstart*(bJstart+gstart))/(astart*bJstart+astart*gstart)];

% Create vectors of inital conditions:
ResJ = linspace(resJmin,resJmax,res0);
ResA = linspace(resAmin,resAmax,res0);
initialJ = find(ResJ>=resJstart,1);
initialA = find(ResA>=resAstart,1);
resJ_start = ResJ(initialJ);
resA_start = ResA(initialA);
indexJ_start = initialJ;
indexA_start = initialA;


% Allow traits to evolve until they reach an evolutionary attractor:
while doneyet==0
    
    % Run an evolutionary simulation:
    [resJ_start,resA_start,init_pop,strain_totalJ,strain_totalA,indexJ_start,indexA_start,RESJnew,RESAnew,~,~] = simulation_function(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJmin,resJmax,resJ_start,resAmin,resAmax,resA_start,h,f,resJmutationprob,init_pop,strain_totalJ,strain_totalA,indexJ_start,indexA_start,res0,nevol,version);
    RESJlengthold=length(RESJ);
    RESAlengthold=length(RESA);
    RESJagain=NaN(RESJlengthold+length(RESJnew),res0);
    RESAagain=NaN(RESAlengthold+length(RESAnew),res0);
    for i=1:RESJlengthold
        RESJagain(i,:)=RESJ(i,:);
    end
    for i=1:length(RESJnew)
        RESJagain(RESJlengthold+i,:)=RESJnew(i,:);
    end
    for i=1:RESAlengthold
        RESAagain(i,:)=RESA(i,:);
    end
    for i=1:length(RESAnew)
        RESAagain(RESAlengthold+i,:)=RESAnew(i,:);
    end
    RESJ=RESJagain;
    RESA=RESAagain;
    
    % Find the endpoints of the trajectories:
    RESJend1=RESJ(end,:);
    RESJend=[0 RESJend1 0];
    [~,locsJ]=findpeaks(RESJend);
    peaknumJ=length(locsJ);
    RESAend1=RESA(end,:);
    RESAend=[0 RESAend1 0];
    [~,locsA]=findpeaks(RESAend);
    peaknumA=length(locsA);
    if peaknumJ==1 && peaknumA==1
        Jout=ResJ(locsJ-1);
        Aout=ResA(locsA-1);
    else
        % The output take the impossible value of 5 if branching has
        % occurred:
        Jout=5;
        Aout=5;
    end
    
    % Determine whether the trajectories have reached an evolutionary
    % attractor. We know that, for the phase plane which this code is being
    % applied to, the region below and to the left of (0.7,0.7) is in one
    % basin of attraction and the region above and to the right of
    % (0.9,0.9) is in the other. 
    if Jout<0.7 && Aout<0.7
        doneyet=1;
    elseif 0.9<Jout && 0.9<Aout
        doneyet=1;
    elseif evotimesteps>50000
        doneyet=1;
    end
    
    evotimesteps=evotimesteps+nevol;

end

% Record evolutionary endpoints:
RESJend1=RESJ(end,:);
RESJend=[0 RESJend1 0];
[~,locsJ]=findpeaks(RESJend);
peaknumJ=length(locsJ);
RESAend1=RESA(end,:);
RESAend=[0 RESAend1 0];
[~,locsA]=findpeaks(RESAend);
peaknumA=length(locsA);
if peaknumJ==1 && peaknumA==1
    Jout=ResJ(locsJ-1);
    Aout=ResA(locsA-1);
else
    disp('Not one endpoint')
    Jout=5;
    Aout=5;
end


end

