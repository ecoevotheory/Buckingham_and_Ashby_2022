% This code can be used to (1) determine that there is a unique, endemic
% equilibrium of our ecological system whenever R0>1 and (2) determine that
% all such endemic equilibria are linearly stable. 

%% Check for stability by running ODEs for a long time

% For many combinations of parameter values, we run the ODE system for
% 't_max' timesteps and determine the range of the values of SJ, SA, IJ and
% IA which are generated near the end of that time period. If the range is
% smaller than 'maxrange' then SJ, SA, IJ and IA have reached a stable
% equilibrium. If not, there may be numerical tolerance issues and so we
% need to use linear stability analysis to check the stability of the
% equilibrium point.

% Note that we can simply check all of the parameter sets using linear
% stability analysis, but this takes a long time. 

% Set up useful quantites and parameters which are not varied:
initvec=[0.1,0.1,0.1,0.1];
t_max=1000;
b0=1;
h=1;
maxrange=0.001;

% We consider a variety of parameter values of the orders of magnitude
% considered in the paper:
version_vec=[1,2,3,4,5,6];
a0_vec=[1,3,5];
g0_vec=[0.5,1,2];
c1a_vec=0.5;
c2a_vec=-3;
c1g_vec=0.5;
c2g_vec=-3;
f_vec=[0.5,1];
beta0_vec=[5,8];
alpha_vec=[0,1,2];
resJ_vec=[0.1,0.3,0.5,0.7,0.9];
resA_vec=[0.1,0.3,0.5,0.7,0.9];

% Create vectors of all possible combinations of these parameter values:
all_combinations=combvec(version_vec,a0_vec,g0_vec,c1a_vec,c2a_vec,c1g_vec,c2g_vec,f_vec,beta0_vec,alpha_vec,resJ_vec,resA_vec);
versionvec=all_combinations(1,:);
a0vec=all_combinations(2,:);
g0vec=all_combinations(3,:);
c1avec=all_combinations(4,:);
c2avec=all_combinations(5,:);
c1gvec=all_combinations(6,:);
c2gvec=all_combinations(7,:);
fvec=all_combinations(8,:);
beta0vec=all_combinations(9,:);
alphavec=all_combinations(10,:);
resJvec=all_combinations(11,:);
resAvec=all_combinations(12,:);

% For each set of parameter values, run the ecological dynamics for a long
% time and see if they tend towards a stable equilibrium:
problemvalues=zeros(length(a0vec),1);
for k=1:length(a0vec)
    
    % Define parameters and trade-off for the given parameter set:
    version=versionvec(k);
    a0=a0vec(k);
    g0=g0vec(k);
    c1a=c1avec(k);
    c2a=c2avec(k);
    c1g=c1gvec(k);
    c2g=c2gvec(k);
    f=fvec(k);
    beta0=beta0vec(k);
    alpha=alphavec(k);
    resJ=resJvec(k);
    resA=resAvec(k);
    betaJ=beta0*(1-resJ);
    betaA=beta0*(1-resA);
    if version==1
        a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)));
        g=g0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        bJ=1;
        bA=1;
    elseif version==2
        bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
        bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
        a=a0;
        g=g0;
    elseif version==3
        bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
        a=a0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        g=g0;
        bJ=1;
    elseif version==4
        g=g0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
        a=a0;
        bJ=1;
    elseif version==5
        bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
        a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)));
        bA=1;
        g=g0;
    elseif version==6
        a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        g=g0;
        bJ=1;
        bA=1;
    end
    
    % Run the ecological dynamics
    [t,x] = fast_ode_function(t_max,a,h,f,g,bJ,bA,betaJ,betaA,alpha,initvec);
    SJ=x(:,1);
    SA=x(:,2);
    IJ=x(:,3);
    IA=x(:,4);
    
    % Determine how much the ecological values are varying by. If they are
    % changing a lot then the equilibrium is unstable:
    nearend=round(9*length(SJ)/10);
    SJ_end=SJ(nearend:end);
    SA_end=SA(nearend:end);
    IJ_end=IJ(nearend:end);
    IA_end=IA(nearend:end);
    if range(SJ_end)>maxrange || range(SA_end)>maxrange || range(IJ_end)>maxrange || range(IA_end)>maxrange
        % Cases where the endemic equilibrium appears to be unstable are
        % recorded in this vector:
        problemvalues(k)=k;
    end

end

problemvalues=nonzeros(problemvalues);

%% Use linear stability analysis on the problematic cases.

% NB: To check that the endemic equilibrium is unique for all sets of
% parameter values, let 
% problemvalues=linspace(1,length(a0vec),length(a0vec)) and then run this
% section only. 

% For each set of parameters, we determine the endemic equilibrium
% analytically, find the Jacobian matrix of the linearised system and
% calculate its eigenvalues (which tell us linear stability).
isproblem=0;
for j=length(problemvalues)
    k=problemvalues(j);
    
    % Define parameter values and trade-offs:
    syms SJ SA IJ IA
    version=versionvec(k);
    a0=a0vec(k);
    g0=g0vec(k);
    c1a=c1avec(k);
    c2a=c2avec(k);
    c1g=c1gvec(k);
    c2g=c2gvec(k);
    f=fvec(k);
    beta0=beta0vec(k);
    alpha=alphavec(k);
    resJ=resJvec(k);
    resA=resAvec(k);

    if version==1
        a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)));
        g=g0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        bJ=1;
        bA=1;
    elseif version==2
        bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
        bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
        a=a0;
        g=g0;
    elseif version==3
        bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
        a=a0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        g=g0;
        bJ=1;
    elseif version==4
        g=g0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
        a=a0;
        bJ=1;
    elseif version==5
        bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
        a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)));
        bA=1;
        g=g0;
    elseif version==6
        a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
        g=g0;
        bJ=1;
        bA=1;
    end
    
    % The original system of equations governing the ecological dynamics of
    % our model is:
    I=IJ+IA;
    N=SJ+SA+IJ+IA;
    kJ=beta0*(1-resJ);
    kA=beta0*(1-resA);
    eq1=a*(1-N)*(SA+f*IA)-(bJ+g+kJ*I)*SJ;
    eq2=g*SJ-(bA+kA*I)*SA;
    eq3=kJ*I*SJ-(bJ+alpha+g)*IJ;
    eq4=g*IJ+kA*I*SA-(bA+alpha)*IA;
    
    % We find the roots of this sytem to get the ecological equilibrium:
    sols=solve([eq1,eq2,eq3,eq4],[SJ,SA,IJ,IA]);
    SJval=double(sols.SJ);
    SAval=double(sols.SA);
    IJval=double(sols.IJ);
    IAval=double(sols.IA);
    
    % Remove any complex roots (SJ, SA, IJ and IA must be real):
    indexvec=linspace(1,length(SJval),length(SJval));
    for i=1:length(SJval)
        if imag(SJval(i))~=0 || imag(SAval(i))~=0 || imag(IJval(i))~=0 || imag(IAval(i))~=0
            indexvec(i)=0;
        elseif SJval(i)<=0 || SAval(i)<=0 || IJval(i)<=0 || IAval(i)<=0
            indexvec(i)=0;
        end
    end
    index=nonzeros(indexvec);
    
    % Calculate the basic reproductive ratio to see if the disease can
    % survive (and hence if there is an endemic equilibrium or only a
    % disease-free one):
    R0=beta0*(a*g-bA*(bJ+g))*((1-resJ)*(bA+alpha+g)*bA+(1-resA)*g*(bJ+alpha+g))/(a*g*(bA+alpha)*(bJ+g)*(bJ+g+alpha));
    
    % Determine whether or not there is a single, endemic equilibrium.
    % Display warning messages if there are multiple endemic equilibria or
    % if there is no equilibrium when the disease should be able to spread.
    if length(index)>1
        disp("More than one endemic equilibrium at k=" +k)
    elseif isempty(index) && R0>1
        disp("Unexpected - no endemic equilibrium at k=" +k)
    elseif length(index)==1
    
        % If there is a single, endemic equilibrium then determine the
        % Jacobian matrix of the linearised system:
        SJstar=SJval(index);
        SAstar=SAval(index);
        IJstar=IJval(index);
        IAstar=IAval(index);
        Istar=IJstar+IAstar;
        Nstar=SJstar+SAstar+IJstar+IAstar;
        
        M11=-a*(SAstar+f*IAstar)-bJ-g-kJ*Istar;
        M12=a*(1-Nstar)-a*(SAstar+f*IAstar);
        M13=-a*(SAstar+f*IAstar);
        M14=a*(1-Nstar)*f-a*(SAstar+f*IAstar);
        M21=g;
        M22=-bA-kA*Istar;
        M23=-kA*SAstar;
        M24=-kA*SAstar;
        M31=kJ*Istar;
        M32=0;
        M33=kJ*SJstar-bJ-alpha-g;
        M34=kJ*SJstar;
        M41=0;
        M42=kA*Istar;
        M43=g+kA*SAstar;
        M44=kA*SAstar-bA-alpha;
        
        M=[M11,M12,M13,M14;M21,M22,M23,M24;M31,M32,M33,M34;M41,M42,M43,M44];
        
        % Find its eigenvalues:
        eigenvalues=double(eig(M));
        if sum(real(eigenvalues)<0)~=length(eigenvalues)
            % If the eigenvalues have positive real part then the
            % equilibrium is unstable:
            disp("Eigenvalues do NOT all have negative real part at k=" +k)
            isproblem=1;
        end
    
    end
    
end

% Display this message if any of the equilibria are unstable:
if isproblem==0
	disp('All equilibria are linearly stable')
end

