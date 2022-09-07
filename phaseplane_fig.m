% This code draws six phase planes for the evolution of juvenile and adult 
% resistance traits for different combinations of trade-offs.

str{1}='J: maturation, A: reproduction';
str{2}='J: juvenile mortality, A: adult mortality';
str{3}='J: reproduction, A: adult mortality';
str{4}='J: maturation, A: adult mortality';
str{5}='J: juvenile mortality, A: reproduction';
str{6}='J: reproduction, A: reproduction';
ORDER = [1,5,6,4,2,3];
labs = {'A','B','C','D','E','F'};

% Define parameters:
a0=5;
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
f=0.1;
h=1;
alpha=0;
beta0=8;
initvec=[0.1,0.1,0.1,0.1];
orig_tmax=100;
matsize=1;
SSres=1000;

for k=1:6
    version=ORDER(k);

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
    
    % Scale the fitness gradient direction arrows to be the same length:
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
    
    % Plot the nullclines and selection gradient directions:
    subplot(2,3,k)
    q=quiver(xvec,yvec,plotter1,plotter2,'color','k');
    q.AutoScaleFactor=0.2;
    xlim([0,1])
    ylim([0,1])
    set(gca,'xtick',0:0.2:1,'xticklabel',resJvec(1):(resJvec(end)-resJvec(1))/5:resJvec(end));
    set(gca,'ytick',0:0.2:1,'yticklabel',resAvec(1):(resAvec(end)-resAvec(1))/5:resAvec(end));
    hold on
    h1=plot(xcurveJ/SSres,ycurveJ/SSres,'o','color','r');
    hold on
    h2=plot(xcurveA/SSres,ycurveA/SSres,'o','color','b');
    text(0.1,0.93,strcat('\bf{',labs{k},'}'),'interpreter','latex','fontsize',20);
    title(str(ORDER(k)),'interpreter','latex','fontsize',12);
    if k==5
        xlabel('Juvenile resistance, $r_J$','interpreter','latex','fontsize',16)
    elseif k==1
        ylabel('Adult resistance, $r_A$','interpreter','latex','fontsize',16)
    end
    axis square
end