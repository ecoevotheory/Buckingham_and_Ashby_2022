function [vec,Jplotvec,Aplotvec,resJss,resAss,disprevvec,juvpropvec,totalpopvec,disprevvec2,juvpropvec2,totalpopvec2]=beta0_fig_data(version)

% This function determines the singular values of adult and juvenile
% resistance, and the proportion of hosts which are infected, for a variety
% of values of the baseline transmissibility.

% Description of model versions:
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

% Define parameters:
a0=5;
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
f=0.5;
h=1;
alpha=0;
initvec=[0.1,0.1,0.1,0.1];
orig_tmax=100;
maxsingstrats=1000;

% We will vary beta0 between 1 and 1000:
vec=logspace(0,2,21);
vec=[vec(1:20),logspace(2,3,5)];
numberof=length(vec);

%% Section 1: Find singular strategies

% Set up vectors to be used later:
resJss=-10+zeros(numberof,maxsingstrats);
resAss=-10+zeros(numberof,maxsingstrats);
Jfitgradcorners=NaN(numberof,4);
Afitgradcorners=NaN(numberof,4);
Jfitgradcornernearone=NaN(numberof,10);
Afitgradcornernearone=NaN(numberof,10);
JfitgradcornerBR=NaN(numberof,3);
AfitgradcornerBR=NaN(numberof,3);
JfitgradcornerTL=NaN(numberof,3);
AfitgradcornerTL=NaN(numberof,3);
JfitgradcornerBL=NaN(numberof,3);
AfitgradcornerBL=NaN(numberof,3);
SSres=100;

% Find singular strategies for each value of the baseline transmissibility:
parfor i=1:length(vec)
    
    disp(i)
    beta0=vec(i);
    
    % To find the singular strategies, we consider a mesh of values of resJ
    % and resA. For each pair, we find the endemic equilibrium and then the
    % two fitness gradients. We then look to see where the sign of these
    % fitness gradients changes.
    
    startJ=0;
    startA=0;
    finJ=1;
    finA=1;
    % Find the fitness gradients:
    [fitgradJval,fitgradAval,resJvalvec,resAvalvec,R0counter]=fitness_gradients(SSres,startJ,startA,finJ,finA,a0,g0,c1a,c2a,c1g,c2g,beta0,h,alpha,f,initvec,version,orig_tmax);
    
    % Record the fitness gradients at the corners of the mesh:
    Jfitgradcorners(i,:)=[fitgradJval(1,1),fitgradJval(end,1),fitgradJval(1,end),fitgradJval(end,end)];
    Afitgradcorners(i,:)=[fitgradAval(1,1),fitgradAval(end,1),fitgradAval(1,end),fitgradAval(end,end)];
    Jfitgradcornernearone(i,:)=[fitgradJval(end,end-1),fitgradJval(end-1,end),fitgradJval(end-1,end-1),fitgradJval(end,end-2),fitgradJval(end-2,end),fitgradJval(end-2,end-2),fitgradJval(end-1,end-2),fitgradJval(end-2,end-1),fitgradJval(end-3,end-3),fitgradJval(end-4,end-4)];
    Afitgradcornernearone(i,:)=[fitgradAval(end,end-1),fitgradAval(end-1,end),fitgradAval(end-1,end-1),fitgradAval(end,end-2),fitgradAval(end-2,end),fitgradAval(end-2,end-2),fitgradAval(end-1,end-2),fitgradAval(end-2,end-1),fitgradAval(end-3,end-3),fitgradAval(end-4,end-4)];
    JfitgradcornerBR(i,:)=[fitgradJval(end,2),fitgradJval(end-1,1),fitgradJval(end-1,2)];
    AfitgradcornerBR(i,:)=[fitgradAval(end,2),fitgradAval(end-1,1),fitgradAval(end-1,2)];
    JfitgradcornerTL(i,:)=[fitgradJval(2,end),fitgradJval(1,end-1),fitgradJval(2,end-1)];
    AfitgradcornerTL(i,:)=[fitgradAval(2,end),fitgradAval(1,end-1),fitgradAval(2,end-1)];
    JfitgradcornerBL(i,:)=[fitgradJval(2,1),fitgradJval(1,2),fitgradJval(2,2)];
    AfitgradcornerBL(i,:)=[fitgradAval(2,1),fitgradAval(1,2),fitgradAval(2,2)];
    
    % This function determines where the fitness gradients change sign:
    [markerJ,markerA,singstratpointerJ,~]=fitgrad_signchange_function(fitgradJval,fitgradAval);
    numberofsingstrats=length(singstratpointerJ);
    
    % If no sign changes are detected:
    if numberofsingstrats==0
        
        % We need to consider the cases where there are no singular
        % strategies. There could still be an evolutionary attractor when
        % one trait is at zero or one:
        [resJss(i,:),resAss(i,:)]=singstrats_at_0or1(fitgradJval,fitgradAval,resJvalvec,resAvalvec,SSres,maxsingstrats,R0counter);

    % Back to the cases where sign changes have been found:
    else
        data=(markerJ'&markerA');
        data=smoothdata(smoothdata(double(data'),1,'gaussian',5),2,'gaussian',5);
        data2 = data';
        [~, locs1, ~, ~] = findpeaks(double(data(:))); % peaks along x
        [~, locs2, ~, ~] = findpeaks(double(data2(:))); % peaks along y
        data_size = size(data); % Gets matrix dimensions
        [col2, row2] = ind2sub(data_size, locs2); % Converts back to 2D indices
        locs2 = sub2ind(data_size, row2, col2);
        ind = intersect(locs1, locs2); % Finds common peak position
        [row, column] = ind2sub(data_size, ind);
        
        % Having found the sign changes, we now 'zoom in' around them to
        % determine the singular strategies more accurately:
        zoomsize=1/SSres;
        Jsingstratvec=-10+zeros(maxsingstrats,1);
        Asingstratvec=-10+zeros(maxsingstrats,1);
        for j=1:length(row)

            resJval=resJvalvec(row(j));
            resAval=resAvalvec(column(j));
            
            startJ = max(0,resJval - zoomsize);
            finJ = min(1,resJval + zoomsize);
            startA = max(0,resAval - zoomsize);
            finA = min(1,resAval + zoomsize);
            
            % Find the fitness gradients:
            [fitgradJvalzoomed,fitgradAvalzoomed,resJvalveczoomed,resAvalveczoomed,~]=fitness_gradients(SSres,startJ,startA,finJ,finA,a0,g0,c1a,c2a,c1g,c2g,beta0,h,alpha,f,initvec,version,orig_tmax);
            
            % Find where they change sign:
            [markerJzoomed,markerAzoomed,~,~]=fitgrad_signchange_function(fitgradJvalzoomed,fitgradAvalzoomed);
            
            data = (markerJzoomed'&markerAzoomed');
            data = smoothdata(smoothdata(double(data'),1,'gaussian',5),2,'gaussian',5);
            data2 = data';
            [~, locs1, ~, ~] = findpeaks(double(data(:))); % peaks along x
            [~, locs2, ~, ~] = findpeaks(double(data2(:))); % peaks along y
            data_size = size(data); % Gets matrix dimensions
            [col2, row2] = ind2sub(data_size, locs2); % Converts back to 2D indices
            locs2 = sub2ind(data_size, row2, col2);
            ind = intersect(locs1, locs2); % Finds common peak position
            [rowzoomed, columnzoomed] = ind2sub(data_size, ind);
            
            % If no sign changes have been detected, we look at the edges
            % of the region again:
            if isempty(rowzoomed)
                [Jsingstrats,Asingstrats]=singstrats_at_0or1(fitgradJvalzoomed,fitgradAvalzoomed,resJvalveczoomed,resAvalveczoomed,SSres,maxsingstrats,R0counter);
                
                % If sign changes have been found, we record their
                % locations and this tells us the values of the singular
                % strategies:
                if ~any(Jsingstrats+10)
                    Jsingstratvec(j)=resJvalvec(round(median(row)));
                    Asingstratvec(j)=resAvalvec(round(median(column)));
                else
                    Jsingstratvec(j)=Jsingstrats(1);
                    Asingstratvec(j)=Asingstrats(1);
                end
                
            else
                Jsingstratvec(j)=resJvalveczoomed(round(mean(rowzoomed)));
                Asingstratvec(j)=resAvalveczoomed(round(mean(columnzoomed)));
            end
            
        end
        
        % Check to see if there are any additional singular strategies at
        % 0 or 1 (where there are no necessarily sign changes in both
        % directions):
        [resJssextra,resAssextra]=singstrats_at_0or1(fitgradJval,fitgradAval,resJvalvec,resAvalvec,SSres,maxsingstrats,R0counter);
        resJssextra(resJssextra==-10)=[];
        resAssextra(resAssextra==-10)=[];
        Jsingstratvec(length(row)+1:length(row)+length(resJssextra))=resJssextra;
        Asingstratvec(length(row)+1:length(row)+length(resAssextra))=resAssextra;
        
        % Create final vectors of singular strategies:
        resJss(i,:)=Jsingstratvec;
        resAss(i,:)=Asingstratvec;
    
    end

end

%% Section 2: When singular strategies are too close together, use analytical methods

% Sometimes, section 1 produces many singular strategies which are very
% close together. This occurs when the resJ and resA nullclines are almost
% parallel, and so numerical tolerances causes several singular startegies
% to be identified when there is, in fact, only one.

% In these cases, we find the singular strategy analytically.

parfor i=1:numberof
    disp(i)
    
    % Identify singular strategies which are very close together:
    Jsingstratvec=-10+zeros(maxsingstrats,1);
    Asingstratvec=-10+zeros(maxsingstrats,1);
    beta0=vec(i);
    resJtemp=resJss(i,:);
    resAtemp=resAss(i,:);
    nonemptyvec=find(resJtemp~=-10);
    ssnumbering=NaN(maxsingstrats,1);
    ssnumbering(resJtemp==-10)=0;
    markernumber=1;
    for j=1:length(nonemptyvec)
        changemarker=0;
        for k=1:length(nonemptyvec) 
            if j~=k && abs(resJtemp(nonemptyvec(j))-resJtemp(nonemptyvec(k)))<0.05 && abs(resAtemp(nonemptyvec(j))-resAtemp(nonemptyvec(k)))<0.05
                if isnan(ssnumbering(nonemptyvec(j))) && isnan(ssnumbering(nonemptyvec(k)))
                    changemarker=1;
                    ssnumbering(nonemptyvec(j))=markernumber;
                    ssnumbering(nonemptyvec(k))=markernumber;
                elseif isnan(ssnumbering(nonemptyvec(j))) && ~isnan(ssnumbering(nonemptyvec(k)))
                    ssnumbering(nonemptyvec(j))=ssnumbering(nonemptyvec(k));
                elseif ~isnan(ssnumbering(nonemptyvec(j))) && isnan(ssnumbering(nonemptyvec(k)))
                    ssnumbering(nonemptyvec(k))=ssnumbering(nonemptyvec(j));
                elseif ~isnan(ssnumbering(nonemptyvec(j))) && ~isnan(ssnumbering(nonemptyvec(k))) && ssnumbering(nonemptyvec(j))~=ssnumbering(nonemptyvec(k))
                    for m=1:length(nonemptyvec)
                        if ssnumbering(nonemptyvec(m))==ssnumbering(nonemptyvec(k)) && m~=k && m~=j
                            ssnumbering(nonemptyvec(m))=ssnumbering(nonemptyvec(j));
                        end
                    end
                    ssnumbering(nonemptyvec(k))=ssnumbering(nonemptyvec(j));
                end
            elseif isnan(ssnumbering(nonemptyvec(j)))
                changemarker=1;
                ssnumbering(nonemptyvec(j))=markernumber;
            end
        end
        if changemarker==1
            markernumber=markernumber+1;
        end
    end
    relabel_ssnumbering=nonzeros(sort(unique(ssnumbering)));
    for j=1:length(relabel_ssnumbering)
        for k=1:length(nonemptyvec)
            if ssnumbering(k)==relabel_ssnumbering(j)
                ssnumbering(k)=j;
            end
        end
    end
    
    % For each cluster of singular strategies, try to find singular
    % strategies in this area using an analytical solver:
    for j=1:max(ssnumbering)
        fallback_marker=0;
        startJ=min(resJtemp(ssnumbering==j));
        finJ=max(resJtemp(ssnumbering==j));
        startA=min(resAtemp(ssnumbering==j));
        finA=max(resAtemp(ssnumbering==j));
        if abs(startJ-finJ)<1e-10 && abs(startA-finA)<1e-10
            newresJ=startJ;
            newresA=startA;
        else
            [newresJ,newresA]=singstrat_finder_analytical(startJ-0.05,finJ+0.05,startA-0.05,finA+0.05,c1g,c2g,c1a,c2a,a0,g0,beta0,h,f,alpha,version);
            
            if isempty(newresJ)
                % Sometimes, the analytical solver doesn't work, so we try
                % it again over a slightly different range:
                [newresJ,newresA]=singstrat_finder_analytical(startJ-0.08,finJ+0.08,startA-0.08,finA+0.08,c1g,c2g,c1a,c2a,a0,g0,beta0,h,f,alpha,version);
            end
            
            if isempty(newresJ)
                disp("No singular strategy found by vpasolver at i= " +i)
                % Where the analytical solver still doesn't work, we use an
                % approximation (see below).
                fallback_marker=1;
            else
                [newresJ1,~]=singstrat_finder_analytical(newresJ+(1e-5),finJ+0.05+(1e-5),startA,finA,c1g,c2g,c1a,c2a,a0,g0,beta0,h,f,alpha,version);
                
                [newresJ2,~]=singstrat_finder_analytical(startJ-0.05-(1e-5),newresJ-(1e-5),startA,finA,c1g,c2g,c1a,c2a,a0,g0,beta0,h,f,alpha,version);
                
                if ~isempty(newresJ1) || ~isempty(newresJ2)
                    disp("Two singular strategies very near each other - maybe more to find - at i= " +i)
                    % We can see from phase planes that this eventuality
                    % should not occur. If it does, we use an approximation
                    % but should also investigate such a case further.
                    fallback_marker=1;
                end
            end
        end
        
        if ~isempty(newresJ)
            Jsingstratvec(j)=newresJ;
            Asingstratvec(j)=newresA;
        elseif fallback_marker==1
            % Where the analytical solver did not work, we use the median 
            % of the 'very close together singular strategies' determined 
            % in section 1 as an approximation.
            Jsingstratvec(j)=median(resJtemp(ssnumbering==j));
            Asingstratvec(j)=median(resAtemp(ssnumbering==j));
        end
    
    end
    
    % The vectors of singular strategies are updated accordingly:
    resJss(i,:)=Jsingstratvec;
    resAss(i,:)=Asingstratvec;
end
        

%% Section 3: Create plot of singular strategies (with classifications)

% Set up vectors to be used later:
aval=[];
gval=[];
bJval=[];
bAval=[];
classificationJ=zeros(numberof,maxsingstrats);
classificationA=zeros(numberof,maxsingstrats);
resJstab=NaN(numberof,maxsingstrats);
resJrep=NaN(numberof,maxsingstrats);
resAstab=NaN(numberof,maxsingstrats);
resArep=NaN(numberof,maxsingstrats);
Jplotvec=NaN(numberof,maxsingstrats*2);
Aplotvec=NaN(numberof,maxsingstrats*2);
deriv2J=NaN(numberof,maxsingstrats);
deriv2A=NaN(numberof,maxsingstrats);

% For each singular strategies, we determine analytically whether or not it
% is evolutionarily stable:
for i=1:numberof
    for j=1:maxsingstrats
        if resJss(i,j)~=-10 && ~isnan(resJss(i,j))
            resJval=resJss(i,j);
            resAval=resAss(i,j);
            beta0=vec(i);
            
            % Define trade-offs:
            if version==1
                aval=a0*(1-c1a*(1-exp(-c2a*resAval))/(1-exp(-c2a)));
                gval=g0*(1-c1g*(1-exp(-c2g*resJval))/(1-exp(-c2g)));
                bJval=1;
                bAval=1;
            elseif version==2
                bJval=1+c1g*(1-exp(-c2g*resJval))/(1-exp(-c2g));
                bAval=1+c1a*(1-exp(-c2a*resAval))/(1-exp(-c2a));
                aval=a0;
                gval=g0;
            elseif version==3
                bAval=1+c1a*(1-exp(-c2a*resAval))/(1-exp(-c2a));
                aval=a0*(1-c1g*(1-exp(-c2g*resJval))/(1-exp(-c2g)));
                gval=g0;
                bJval=1;
            elseif version==4
                gval=g0*(1-c1g*(1-exp(-c2g*resJval))/(1-exp(-c2g)));
                bAval=1+c1a*(1-exp(-c2a*resAval))/(1-exp(-c2a));
                aval=a0;
                bJval=1;
            elseif version==5
                bJval=1+c1g*(1-exp(-c2g*resJval))/(1-exp(-c2g));
                aval=a0*(1-c1a*(1-exp(-c2a*resAval))/(1-exp(-c2a)));
                bAval=1;
                gval=g0;
            elseif version==6
                aval=a0*(1-c1a*(1-exp(-c2a*resAval))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJval))/(1-exp(-c2g)));
                gval=g0;
                bJval=1;
                bAval=1;
            end
            
            % Find endemic equilibrium:
            [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(resJval,resAval,aval,gval,bJval,bAval,beta0,h,f,alpha,initvec,orig_tmax);
            
            % Determine conditions for evolutionary stability:
            if version==1
                deriv2J(i,j)=(a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((c1g*c2g^2*g0*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g^2*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))^2*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^3*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) + (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) + (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*c1g^2*c2g^2*g0^2*exp(-2*c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)^2*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2) - (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*c1g*c2g^2*g0^2*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2) + (a0*c1g*c2g^2*g0^2*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1g*c2g*g0^2*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) + (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2) - (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*c1g^2*c2g^2*g0^3*exp(-2*c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)^2*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^3) - (2*a0*c1g*c2g*g0^2*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2);
                deriv2A(i,j)=-(2*a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0^2*f*(IAval + IJval*h)^2*(resJval - 1))*(IAval + IJval*h)*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*beta0^2*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)^2*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^3*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0^2*f*(IAval + IJval*h)^2*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*beta0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1));
            elseif version==2
                deriv2J(i,j)=(2*a0*g0*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) + (c1g*c2g*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*g0*((c1g*c2g^2*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g^2*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*g0*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))^2*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^3*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) + (c1g*c2g*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2) - (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2) - (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*c1g^2*c2g^2*g0*exp(-2*c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)^2*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^3) + (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2);
                deriv2A(i,j)=(2*a0*g0*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*g0*((c1a*c2a^2*exp(-c2a*resAval)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))/(exp(-c2a) - 1) - (beta0*c1a*c2a^2*f*exp(-c2a*resAval)*(IAval + IJval*h)*(resJval - 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*g0*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))^2*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^3*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*(beta0*f*(IAval + IJval*h)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*c1a^2*c2a^2*g0*exp(-2*c2a*resAval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)^2*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^3*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1));
            elseif version==3
                deriv2J(i,j)=(2*a0*beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (2*a0*beta0^2*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)^2*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^3*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (2*a0*beta0^2*f*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)^2*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^2*(alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) + (2*a0*beta0*c1g*c2g*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^2*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1));
                deriv2A(i,j)=(a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((c1a*c2a^2*exp(-c2a*resAval)*(alpha + g0 + 1))/(exp(-c2a) - 1) - (beta0*c1a*c2a^2*f*exp(-c2a*resAval)*(IAval + IJval*h)*(resJval - 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (2*a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(alpha + g0 + 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha + g0 + 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (2*a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))^2*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^3*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (2*a0*c1a^2*c2a^2*g0*exp(-2*c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)^2*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^3) - (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha + g0 + 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha + g0 + 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2) + (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2) + (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2);
            elseif version==4
                deriv2J(i,j)=(a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((c1g*c2g^2*g0*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g^2*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))^2*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^3*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1g^2*c2g^2*g0^2*exp(-2*c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)^2*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2) + (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1g^2*c2g^2*g0^3*exp(-2*c2g*resJval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)^2*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^3) + (a0*c1g*c2g^2*g0^2*exp(-c2g*resJval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2) - (a0*c1g*c2g^2*g0^2*exp(-c2g*resJval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (2*a0*c1g*c2g*g0^2*exp(-c2g*resJval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2) + (2*a0*c1g*c2g*g0^2*exp(-c2g*resJval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2);
                deriv2A(i,j)=(2*a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))/(exp(-c2a) - 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((c1a*c2a^2*exp(-c2a*resAval)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))/(exp(-c2a) - 1) - (beta0*c1a*c2a^2*f*exp(-c2a*resAval)*(IAval + IJval*h)*(resJval - 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))^2*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^3*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))/(exp(-c2a) - 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1a^2*c2a^2*g0*exp(-2*c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)^2*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^3*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1));
            elseif version==5
                deriv2J(i,j)=(a0*c1g*c2g^2*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2) - (2*a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))^2*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^3*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (2*a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) - (c1g*c2g*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (2*a0*c1g^2*c2g^2*g0*exp(-2*c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)^2*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^3) - (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) - (c1g*c2g*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2) - (a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*c2g^2*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g^2*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (2*a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2);
                deriv2A(i,j)=(2*a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0^2*f*(IAval + IJval*h)^2*(resJval - 1) - beta0*f*(IAval + IJval*h)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval*h)*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (2*a0*beta0^2*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(IAval + IJval*h)^2*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^3*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*(beta0^2*f*(IAval + IJval*h)^2*(resJval - 1) - beta0*f*(IAval + IJval*h)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (2*a0*beta0*c1a*c2a*g0*exp(-c2a*resAval)*(IAval + IJval*h)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1));
            elseif version==6
                deriv2J(i,j)=(2*a0*beta0^2*f*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)^2*(IAval + IJval + SAval + SJval - 1))/((alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^2*(alpha + g0 + 1)) + (2*a0*beta0^2*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)^2*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^3*(alpha + g0 + 1)) + (a0*c1g*c2g^2*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) - (2*a0*beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(IAval + IJval*h)*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) - (2*a0*beta0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^2*(alpha + g0 + 1));
                deriv2A(i,j)=(2*a0*beta0^2*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)^2*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^3*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) + (2*a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*(beta0*f*(IAval + IJval*h)*(alpha + g0 + 1) - beta0^2*f*(IAval + IJval*h)^2*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) + (2*a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha + g0 + 1) - beta0^2*f*(IAval + IJval*h)^2*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) + (a0*c1a*c2a^2*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) + (2*a0*beta0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1));
            end
            
            % Classification 1 represents an attractor and classification 2
            % represents a repeller:
            if deriv2J(i,j)<0 && deriv2A(i,j)<0
                classificationJ(i,j)=1;
                classificationA(i,j)=1;
            elseif (resJval==0 || resJval==1 || resAval==0 ||resAval==1) && deriv2J(i,j)<40 && deriv2A(i,j)<40
                classificationJ(i,j)=1;
                classificationA(i,j)=1;
            else
                classificationJ(i,j)=2;
                classificationA(i,j)=2;
            end 
        end
        
        % It is possible, but very rare, for the juvenile and adult
        % singular strategies to have different classifications:
        if classificationJ(i,j)~=classificationA(i,j)
            disp('Different classifications for juveniles and adults')
        end
    end
    
    % Very occassionally, a repeller near the edge is at the intersection
    % of nullclines which are almost parallel to the edge, and it is
    % mis-identified as an attractor. If the derivatives at this value are
    % very positive, then this is likely:
    for j=1:maxsingstrats-1
        if isequal(classificationJ(i,j:j+1),[1,1]) && (deriv2J(i,j)>1 || deriv2A(i,j)>1)
            classificationJ(i,j)=2;
            classificationA(i,j)=2;
        elseif isequal(classificationJ(i,j:j+1),[1,1]) && (deriv2J(i,j+1)>1 || deriv2A(i,j+1)>1)
            classificationJ(i,j+1)=2;
            classificationA(i,j+1)=2;
        end
    end
    
    % Sort the vectors ready for plotting:
    stabs=find(classificationJ(i,:)==1);
    if ~isempty(stabs)
        resJstab(i,1:length(stabs))=sort(resJss(i,stabs));
        resAstab(i,1:length(stabs))=sort(resAss(i,stabs));
    end

    % There may also be evolutionary attractors or repellers at the corners
    % of our domain (the plane where juvenile resistance is between 0 and 1
    % and where adult resistance is between 0 and 1):
    if Jfitgradcorners(i,1)<0 && Afitgradcorners(i,1)<0
        resJrep(i,1)=0;
        resArep(i,1)=0;
    end
    if Jfitgradcorners(i,2)>0 && Afitgradcorners(i,2)<0
        resJrep(i,2)=1;
        resArep(i,1)=0;
    end  
    if Jfitgradcorners(i,3)<0 && Afitgradcorners(i,3)>0
        resJrep(i,1)=0;
        resArep(i,2)=1;
    end      
    if Jfitgradcorners(i,4)>0 && Afitgradcorners(i,4)>0
        resJrep(i,2)=1;
        resArep(i,2)=1;
    end 
    
    % Sort vectors ready for plotting:
    Jplotvec(i,1:length(sort(unique([resJstab(i,:),resJrep(i,:)]))))=sort(unique([resJstab(i,:),resJrep(i,:)]));
    Aplotvec(i,1:length(sort(unique([resAstab(i,:),resArep(i,:)]))))=sort(unique([resAstab(i,:),resArep(i,:)]));
    
    % Sometimes, the nullclines which are 'on top of each other' near
    % (1,1) are so close together and so close to (1,1) that their
    % attractor does not register:
    Jvecsizes=Jplotvec(i,:);
    Jvecsizes(isnan(Jvecsizes))=[];
    Avecsizes=Aplotvec(i,:);
    Avecsizes(isnan(Avecsizes))=[];
    Jnearone=Jfitgradcornernearone(i,:);
    Jnearone(isnan(Jnearone))=[];
    Anearone=Afitgradcornernearone(i,:);
    Anearone(isnan(Anearone))=[];
    JBR=JfitgradcornerBR(i,:);
    JBR(isnan(JBR))=[];
    ABR=AfitgradcornerBR(i,:);
    ABR(isnan(ABR))=[];
    JTL=JfitgradcornerTL(i,:);
    JTL(isnan(JTL))=[];
    ATL=AfitgradcornerTL(i,:);
    ATL(isnan(ATL))=[];
    JBL=JfitgradcornerBL(i,:);
    JBL(isnan(JBL))=[];
    ABL=AfitgradcornerBL(i,:);
    ABL(isnan(ABL))=[];
    % Where this occurs, we add in the attractor:
    if ~isempty(Jvecsizes) && ~isempty(Avecsizes) && ~isempty(Jnearone) && ~isempty(Anearone)
        if (all(Jvecsizes<0.9) || all(Avecsizes<0.9)) && ~all(Jnearone<0) && ~all(Anearone<0)
            resJrep(i,2)=0.99;
            resArep(i,2)=0.99;
        end
    end
    if ~isempty(Jvecsizes) && ~isempty(Avecsizes) && ~isempty(JBR) && ~isempty(ABR)
        if (all(Jvecsizes<0.9) || all(Avecsizes>0.1)) && ~all(JBR<0) && ~all(ABR>0)
            resJrep(i,2)=0.99;
            resArep(i,2)=0.01;
        end
    end
    if ~isempty(Jvecsizes) && ~isempty(Avecsizes) && ~isempty(JTL) && ~isempty(ATL)
        if (all(Jvecsizes>0.1) || all(Avecsizes<0.9)) && ~all(JTL>0) && ~all(ATL<0)
            resJrep(i,2)=0.01;
            resArep(i,2)=0.99;
        end
    end
    if ~isempty(Jvecsizes) && ~isempty(Avecsizes) && ~isempty(JBL) && ~isempty(ABL)
        if (all(Jvecsizes>0.1) || all(Avecsizes>0.1)) && ~all(JBL>0) && ~all(ABL>0)
            resJrep(i,2)=0.01;
            resArep(i,2)=0.01;
        end
    end
    
    % Sort vectors ready for plotting:
    Jplotvec(i,1:length(sort(unique([resJstab(i,:),resJrep(i,:)]))))=sort(unique([resJstab(i,:),resJrep(i,:)]));
    Aplotvec(i,1:length(sort(unique([resAstab(i,:),resArep(i,:)]))))=sort(unique([resAstab(i,:),resArep(i,:)]));

    
end
    

%% Section 4: Find disease prevalence and juvenile proportion

% Set up vectors to be used later:
disprevvec=NaN(numberof,1);
juvpropvec=NaN(numberof,1);
totalpopvec=NaN(numberof,1);
disprevvec2=NaN(numberof,1);
juvpropvec2=NaN(numberof,1);
totalpopvec2=NaN(numberof,1);
aval=[];
gval=[];
bJval=[];
bAval=[];
syms SJ SA IJ IA
        
% For each attractor, we determine the ecological equilibrium at that
% attractor and hence calculate the disease prevalence and juvenile
% proportion at the ecological equilibrium:
for j=1:numberof

    resJstart=Jplotvec(j,1);
    resAstart=Aplotvec(j,1);
    beta0=vec(j);
    
    if (~isnan(resJstart) && ~isnan(resAstart))
        
        % Set up trade-offs for different model versions:
        if version==1
            aval=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)));
            gval=g0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            bJval=1;
            bAval=1;
        elseif version==2
            bJval=1+c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g));
            bAval=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
            aval=a0;
            gval=g0;
        elseif version==3
            bAval=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
            aval=a0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            gval=g0;
            bJval=1;
        elseif version==4
            gval=g0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            bAval=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
            aval=a0;
            bJval=1;
        elseif version==5
            bJval=1+c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g));
            aval=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)));
            bAval=1;
            gval=g0;
        elseif version==6
            aval=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            gval=g0;
            bJval=1;
            bAval=1;
        end
        
        % Determine endemic equilibrium:
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(resJstart,resAstart,aval,gval,bJval,bAval,beta0,h,f,alpha,initvec,orig_tmax);
        
        % Define disease prevalence, juvenile proportion and the
        % total population density:
        disprevvec(j)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        juvpropvec(j)=(IJval+SJval)/(SJval+SAval+IJval+IAval);
        totalpopvec(j)=SJval+SAval+IJval+IAval;
    end
    
    % If there is bistability, we now determine the disease prevalence,
    % juvenile proportion and total population density for the second
    % set of co-singular strategies:
    resJstart=NaN;
    resAstart=NaN;
    if ~isnan(Jplotvec(j,2)) && ~isnan(Aplotvec(j,2))
        resJstart=Jplotvec(j,2);
        resAstart=Aplotvec(j,2);
    elseif ~isnan(Jplotvec(j,2))
        resJstart=Jplotvec(j,2);
        resAstart=Aplotvec(j,1);
    elseif ~isnan(Aplotvec(j,2))
        resJstart=Jplotvec(j,1);
        resAstart=Aplotvec(j,2);
    end
    
    if (~isnan(resJstart) && ~isnan(resAstart))
        
        % Set up trade-offs for different model versions:
        if version==1
            aval=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)));
            gval=g0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            bJval=1;
            bAval=1;
        elseif version==2
            bJval=1+c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g));
            bAval=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
            aval=a0;
            gval=g0;
        elseif version==3
            bAval=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
            aval=a0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            gval=g0;
            bJval=1;
        elseif version==4
            gval=g0*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            bAval=1+c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a));
            aval=a0;
            bJval=1;
        elseif version==5
            bJval=1+c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g));
            aval=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)));
            bAval=1;
            gval=g0;
        elseif version==6
            aval=a0*(1-c1a*(1-exp(-c2a*resAstart))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJstart))/(1-exp(-c2g)));
            gval=g0;
            bJval=1;
            bAval=1;
        end
        
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(resJstart,resAstart,aval,gval,bJval,bAval,beta0,h,f,alpha,initvec,orig_tmax);
        
        % Define disease prevalence, juvenile proportion and the
        % proportions of juveniles and adults which are infected:
        disprevvec2(j)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        juvpropvec2(j)=(IJval+SJval)/(SJval+SAval+IJval+IAval);
        totalpopvec2(j)=SJval+SAval+IJval+IAval;
    end
        

end


end
