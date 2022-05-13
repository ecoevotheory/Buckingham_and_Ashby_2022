% This code draws a phase plane for the evolution of juvenile and adult 
% resistance traits.

%% Section 1: Draw the nullclines and fitness gradient directions

% Define parameters:
version=2;
a0=5;
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
f=0.5;
h=1;
alpha=0;
beta0=1000;
initvec=[0.1,0.1,0.1,0.1];
orig_tmax=100;
matsize=1;
SSres=1000;

% Calculate the fitness gradients at different values of the resistance
% traits:
startJ=0;
startA=0;
finJ=1;
finA=1;
[fitgradJval,fitgradAval,resJvalvec,resAvalvec,R0counter]=fitness_gradients(SSres,startJ,startA,finJ,finA,a0,g0,c1a,c2a,c1g,c2g,beta0,h,alpha,f,initvec,version,orig_tmax);

% Set up vectors:
resJvec=resJvalvec;
resAvec=resAvalvec;
xvec=zeros(length(resJvec)*length(resAvec),1);
yvec=zeros(length(resJvec)*length(resAvec),1);
plotter1=zeros(length(resJvec)*length(resAvec),1);
plotter2=zeros(length(resJvec)*length(resAvec),1);
xcurveJ=zeros(length(resJvec)*length(resAvec),1);
ycurveJ=zeros(length(resAvec)*length(resJvec),1);
xcurveA=zeros(length(resJvec)*length(resAvec),1);
ycurveA=zeros(length(resAvec)*length(resJvec),1);
xnanvec=zeros(length(resAvec)*length(resJvec),1);
ynanvec=zeros(length(resAvec)*length(resJvec),1);

% Convert data into vectors to be plotted. Our phase plane will show the
% direction but not the magnitude of the fitness gradients (because the
% magnitudes vary so much that, if we showed the high magnitudes, we would
% not be able to see the low ones at all).
for i=1:length(resJvec)
    for j=1:length(resAvec)
        xvec(i+length(resJvec)*(j-1))=resJvec(i);
        yvec(i+length(resJvec)*(j-1))=resAvec(j);
        plotter1(i+length(resJvec)*(j-1))=fitgradJval(i,j);
        plotter2(i+length(resJvec)*(j-1))=fitgradAval(i,j);
        if isnan(fitgradJval(i,j)) && isnan(fitgradAval(i,j))
            xnanvec(i+length(resJvec)*(j-1))=i;
            ynanvec(i+length(resJvec)*(j-1))=j;
        end
        if j<length(resAvec)
            if fitgradJval(i,j)>0 && fitgradJval(i,j+1)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i,j+1)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        if j>1
            if fitgradJval(i,j)>0 && fitgradJval(i,j-1)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i,j-1)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        if i<length(resJvec)
            if fitgradJval(i,j)>0 && fitgradJval(i+1,j)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i+1,j)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        if i>1
            if fitgradJval(i,j)>0 && fitgradJval(i-1,j)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i-1,j)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        
        if j>1 && i>1
            if fitgradJval(i,j)>0 && fitgradJval(i-1,j-1)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i-1,j-1)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        if j<length(resAvec) && i<length(resJvec)
            if fitgradJval(i,j)>0 && fitgradJval(i+1,j+1)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i+1,j+1)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        if j<length(resAvec) && i>1
            if fitgradJval(i,j)>0 && fitgradJval(i-1,j+1)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i-1,j+1)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        if j>1 && i<length(resJvec)
            if fitgradJval(i,j)>0 && fitgradJval(i+1,j-1)<=0
                xcurveJ(i+length(resJvec)*(j-1))=i;
                ycurveJ(i+length(resJvec)*(j-1))=j;
            end
            if fitgradAval(i,j)>0 && fitgradAval(i+1,j-1)<=0
                xcurveA(i+length(resJvec)*(j-1))=i;
                ycurveA(i+length(resJvec)*(j-1))=j;
            end
        end
        
    end
end

% Scale the fitness gradient directions to be the same length:
xcurveJ=nonzeros(xcurveJ);
ycurveJ=nonzeros(ycurveJ);
xcurveA=nonzeros(xcurveA);
ycurveA=nonzeros(ycurveA);
for i=1:length(resJvec)*length(resAvec)
    scaler=sqrt(plotter1(i)^2 + plotter2(i)^2);
    if scaler~=0
        plotter1(i)=plotter1(i)/scaler;
        plotter2(i)=plotter2(i)/scaler;
    end  
end

% Refine the vectors so that an appropriate number of fitness gradients are
% plotted:
deletevec=zeros(length(resJvec)*length(resAvec),1);
for i=1:length(resJvec)
    for j=1:length(resAvec)
        if floor((i-(SSres/20))/(SSres/10))~=(i-(SSres/20))/(SSres/10)
            deletevec(i+length(resJvec)*(j-1))=i+length(resJvec)*(j-1);
        end
        if floor((j-(SSres/20))/(SSres/10))~=(j-(SSres/20))/(SSres/10)
            deletevec(i+length(resJvec)*(j-1))=i+length(resJvec)*(j-1);
        end
    end
end
deletevec=nonzeros(deletevec);
xvec(deletevec)=[];
yvec(deletevec)=[];
plotter1(deletevec)=[];
plotter2(deletevec)=[];


% Plot the phase plane:
q=quiver(xvec,yvec,plotter1,plotter2,'color','k');
q.AutoScaleFactor=0.2;
xlim([0,1])
ylim([0,1])
set(gca,'xtick',0:0.2:1,'xticklabel',resJvec(1):(resJvec(end)-resJvec(1))/5:resJvec(end));
set(gca,'ytick',0:0.2:1,'yticklabel',resAvec(1):(resAvec(end)-resAvec(1))/5:resAvec(end));
hold on
h1=plot(xcurveJ/SSres,ycurveJ/SSres,'o','color','r');
hold on
h2=plot(xcurveA/SSres,ycurveA/SSres,'x','color','b');
ax=gca;
ax.FontSize=14;
xlabel('Juvenile resistance, $r_J$','interpreter','latex','fontsize',16)
ylabel('Adult resistance, $r_A$','interpreter','latex','fontsize',16)

% Plot the singular strategies:
plot(0,0,'.','MarkerSize',50,'color',[0.4,0.1,0.5])
plot(1,1,'.','MarkerSize',50,'color',[0.4,0.1,0.5])
axis square

%% Section 2: Plot the boundary of the basin of attraction of each singular strategy

% Stipulate the precision (level of detail) you want:
lod=7;
resAfirststart=0.5;
resAsecondstart=0.8;
resJfirststart=0.8;
resJsecondstart=1;

% Each step in this 'while' loop zooms in to add greater precision:
while matsize<2^lod+2
    disp(matsize)
    resAnew=[];
    resJnew=[];
    
    % We seek to determine which basin of attraction each point in our
    % plane lies in. Other than on the first step, when no values have 
    % already been determined, we first assign values which are not near
    % the boundary with the same value as the points around them:
    if matsize~=1
        thematrixonly=thematrix;
        % Create vectors of the new values of resJ and resA which we are
        % going to consider as we zoom in. These are the mid-points of the
        % values which have been considered previously.
        resAnew=zeros(1,length(resAold)-1);
        resJnew=zeros(1,length(resJold)-1);
        for i=1:length(resAold)-1
            resAnew(i)=(resAold(i)+resAold(i+1))/2;
        end
        for i=1:length(resJold)-1
            resJnew(i)=(resJold(i)+resJold(i+1))/2;
        end
        % These are added to the values considered previously to get a new
        % set of parameter values.
        resAvec=sort([resAold,resAnew]);
        resJvec=sort([resJold,resJnew]);
        % Define useful quantities and set up vectors to be used later:
        numberofA=length(resAvec);
        numberofJ=length(resJvec);
        endpoint=zeros(numberofA,numberofJ);
        indicator=zeros(length(resAvec),length(resJvec));
        % If we have already found which basin of attraction a particular 
        % point lies in, then we can just use that previous result:
        for i=1:length(resAold)
            for j=1:length(resJold)
                endpoint(2*i-1,2*j-1)=thematrixonly(i,j);
                indicator(2*i-1,2*j-1)=1;
            end
        end
        % If the new parameters lie 'in the middle' of a square of values
        % which have already been considered, then the value here will be
        % the same as that in the corners, if all of the corners display
        % the same behaviour.
        for i=1:length(resAnew)
            for j=1:length(resJnew)
                if thematrixonly(i,j)==thematrixonly(i,j+1) && thematrixonly(i,j)==thematrixonly(i+1,j) && thematrixonly(i,j)==thematrixonly(i+1,j+1)
                    endpoint(2*i,2*j)=thematrixonly(i,j);
                    indicator(2*i,2*j)=1;
                end
            end
        end
        % If the new parameter values are 'on the edge' of a square of
        % values already considered, then the value here will be the same
        % as those on either side of it if all four corners of both squares
        % which neighbour it all have the same value.
        for i=1:length(resAnew)
            for j=1:length(resJold)
                if j==1
                    if thematrixonly(i,j)==thematrixonly(i,j+1) && thematrixonly(i,j)==thematrixonly(i+1,j) && thematrixonly(i,j)==thematrixonly(i+1,j+1)
                        endpoint(2*i,2*j-1)=thematrixonly(i,j);
                        indicator(2*i,2*j-1)=1;
                    end
                end
                if j==length(resJold)
                    if thematrixonly(i,j)==thematrixonly(i,j-1) && thematrixonly(i,j)==thematrixonly(i+1,j) && thematrixonly(i,j)==thematrixonly(i+1,j-1)
                        endpoint(2*i,2*j-1)=thematrixonly(i,j);
                        indicator(2*i,2*j-1)=1;
                    end
                end   
                if j>1 && j<length(resJold)
                    if thematrixonly(i,j)==thematrixonly(i,j+1) && thematrixonly(i,j)==thematrixonly(i,j-1) && thematrixonly(i,j)==thematrixonly(i+1,j) && thematrixonly(i,j)==thematrixonly(i+1,j+1) && thematrixonly(i,j)==thematrixonly(i+1,j-1)
                        endpoint(2*i,2*j-1)=thematrixonly(i,j);
                        indicator(2*i,2*j-1)=1;
                    end
                end   
            end
        end
        for i=1:length(resAold)
            for j=1:length(resJnew)
                if i==1
                    if thematrixonly(i,j)==thematrixonly(i,j+1) && thematrixonly(i,j)==thematrixonly(i+1,j) && thematrixonly(i,j)==thematrixonly(i+1,j+1)
                        endpoint(2*i-1,2*j)=thematrixonly(i,j);
                        indicator(2*i-1,2*j)=1;
                    end
                end
                if i==length(resAold)
                    if thematrixonly(i,j)==thematrixonly(i,j+1) && thematrixonly(i,j)==thematrixonly(i-1,j) && thematrixonly(i,j)==thematrixonly(i-1,j+1)
                        endpoint(2*i-1,2*j)=thematrixonly(i,j);
                        indicator(2*i-1,2*j)=1;
                    end
                end   
                if i>1 && i<length(resAold)
                    if thematrixonly(i,j)==thematrixonly(i-1,j) && thematrixonly(i,j)==thematrixonly(i,j+1) && thematrixonly(i,j)==thematrixonly(i+1,j) && thematrixonly(i,j)==thematrixonly(i-1,j+1) && thematrixonly(i,j)==thematrixonly(i+1,j+1)
                        endpoint(2*i-1,2*j)=thematrixonly(i,j);
                        indicator(2*i-1,2*j)=1;
                    end
                end   
            end
        end
    end
    
    % Define useful quantities and set up vectors to be used later:
    if matsize==1
        resAvec=[resAfirststart,resAsecondstart];
        resJvec=[resJfirststart,resJsecondstart];
        numberofA=length(resAvec);
        numberofJ=length(resJvec);
        endpoint=zeros(numberofA,numberofJ);
        indicator=zeros(length(resAvec),length(resJvec));
    end
    
    % For each combination of parameters, we run a simulation of the
    % evolution of the adult and juvenile resistance traits:
    for m=1:numberofA
        for n=1:numberofJ
            % If we have not already determined the basin of attraction for
            % this point, then we proceed:
            if indicator(m,n)==0
                resA=resAvec(m);
                resJ=resJvec(n);
        
                % Run the simulation to determine the end-point of the
                % trajectory from each starting value:
                [Jout,Aout]=sim_inout_function(resJ,resA,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
                if Jout<0.7 && Aout<0.7
                    endpoint(m,n)=1;
                elseif 0.9<Jout && 0.9<Aout 
                    endpoint(m,n)=2;
                else
                    disp("Simulation not finished at resJ=" +resJ +"and resA=" +resA)
                end
                
            end
        end
    end
    
    % Set up the matrix which we are going to use later:
    thematrix=endpoint;
    firstcolumn=transpose([0,resAvec]);
    matsize=length(thematrix)+1;
    % Record which parameter values we have already considered:
    resAold=resAvec;
    resJold=resJvec;
    % If matsize is smaller than the required value, then we need to zoom
    % in further. As such, we repeat the 'while' loop.
end
boa_matrix=thematrix;


% We now know where the basins of attraction are. We want to plot their
% boundary as a curve:
boa_curve=NaN(numberofA,1);
for m=1:numberofA
    switchvec=NaN(numberofJ,1);
    for n=1:numberofJ
        if n<numberofJ-1
            if boa_matrix(m,n)==1 && boa_matrix(m,n+1)==2
                switchvec(n)=1;
            end
        else
            if boa_matrix(m,n)==1
                switchvec(n)=1;
            end
        end
    end
    boundaryval=find(switchvec==1,1,'last');
    if ~isempty(boundaryval)
        boa_curve(m)=resJvec(boundaryval);
    else
        boa_curve(m)=NaN;
    end
    boa_curve(boa_curve>0.98)=NaN;
        
end

% Smooth the data and create the plot:
smooth_boa_curve=smoothdata(boa_curve,'movmean',10); 
hold on
grey=[0.6,0.6,0.7];
plot(smooth_boa_curve,resAvec,'--','color',grey,'linewidth',2)


%% Section 3: Plot some trajectories on the phase plane

% Pick a starting point and run a simulation:
resJstart=0.1;
resAstart=1;
[trajJ,trajA]=sim_traj_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
hold on
% Smooth the data and plot the trajectory:
smooth_trajJ=smoothdata(trajJ,'movmean',100);
smooth_trajA=smoothdata(trajA,'movmean',100);
smooth_trajJ(end+1)=0;
smooth_trajA(end+1)=0;
plot(smooth_trajJ,smooth_trajA,'linewidth',2,'color',[0.3,0.6,0.3])

% Now choose some more starting points:
resJstart=1;
resAstart=0.4;
[trajJ2,trajA2]=sim_traj_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
hold on
smooth_trajJ2=smoothdata(trajJ2,'movmean',100);
smooth_trajA2=smoothdata(trajA2,'movmean',100);
plot(smooth_trajJ2,smooth_trajA2,'linewidth',2,'color',[0.3,0.6,0.3])

resJstart=0.4;
resAstart=0.8;
[trajJ3,trajA3]=sim_traj_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
hold on
smooth_trajJ3=smoothdata(trajJ3,'movmean',100);
smooth_trajA3=smoothdata(trajA3,'movmean',100);
plot(smooth_trajJ3,smooth_trajA3,'linewidth',2,'color',[0.3,0.6,0.3])

resJstart=0.8;
resAstart=0.8;
[trajJ4,trajA4]=sim_traj_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
hold on
smooth_trajJ4=smoothdata(trajJ4,'movmean',100);
smooth_trajA4=smoothdata(trajA4,'movmean',100);
plot(smooth_trajJ4,smooth_trajA4,'linewidth',2,'color',[0.3,0.6,0.3])


resJstart=0.85;
resAstart=0.82;
[trajJ5,trajA5]=sim_traj_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
hold on
smooth_trajJ5=smoothdata(trajJ5,'movmean',100);
smooth_trajA5=smoothdata(trajA5,'movmean',100);
plot(smooth_trajJ5,smooth_trajA5,'linewidth',2,'color',[0.3,0.6,0.3])

resJstart=0.7;
resAstart=0.9;
[trajJ6,trajA6]=sim_traj_function(resJstart,resAstart,version,a0,g0,c1a,c2a,c1g,c2g,f,h,alpha,beta0);
hold on
smooth_trajJ6=smoothdata(trajJ6,'movmean',100);
smooth_trajA6=smoothdata(trajA6,'movmean',100);
plot(smooth_trajJ6,smooth_trajA6,'linewidth',2,'color',[0.3,0.6,0.3])

% Indicate the starting-points of the trajectories:
hold on
plot(smooth_trajJ(1),smooth_trajA(1),'o','markersize',12,'color',[0.3,0.6,0.3],'linewidth',4)
plot(smooth_trajJ2(1),smooth_trajA2(1),'o','markersize',12,'color',[0.3,0.6,0.3],'linewidth',4)
plot(smooth_trajJ3(1),smooth_trajA3(1),'o','markersize',12,'color',[0.3,0.6,0.3],'linewidth',4)
plot(smooth_trajJ4(1),smooth_trajA4(1),'o','markersize',12,'color',[0.3,0.6,0.3],'linewidth',4)
plot(smooth_trajJ5(1),smooth_trajA5(1),'o','markersize',12,'color',[0.3,0.6,0.3],'linewidth',4)
plot(smooth_trajJ6(1),smooth_trajA6(1),'o','markersize',12,'color',[0.3,0.6,0.3],'linewidth',4)

