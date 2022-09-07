% This code creates four sub-figures showing the effect of varying the 
% baseline transmissibility, beta0, for different sets of parameters.

% The red curve shows the juvenile resistance and the blue curve shows the
% adult resistance. The dashed, grey curve shows the proportion of the
% population which is infected.

% Stipulate the size of the figure and its position on the screen:
figure(100)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 16; ySize = 7;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

% Parameter values: 
a0=5;
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
alpha=0;

for i=1:4

    % We want this plot to show the different possible shapes of the
    % relationship between beta0 and the singular values of resistance. We
    % therefore choose the following parameters:
    if i==1
        version=2;
        f=0.3;
    elseif i==2
        version=6;
        f=0.5;
    elseif i==3
        version=2;
        f=0.5;
    elseif i==4
        version=6;
        f=0.3;
    end
    
    % The function 'beta0_fig_data' determines the singular strategies and
    % disease prevalence for a variety of values of the baseline
    % transmissibility.
    [vec,Jplotvec,Aplotvec,disprevvec1,totalpopvec1,disprevvec2,totalpopvec2]=beta0_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,f,alpha);
    upperx=max(log10(vec));
    diseasedensity1=disprevvec1.*totalpopvec1;
    diseasedensity2=disprevvec2.*totalpopvec2;
    Jplotvecpart1=Jplotvec(:,1);
    Jplotvecpart2=Jplotvec(:,2);
    Aplotvecpart1=Aplotvec(:,1);
    Aplotvecpart2=Aplotvec(:,2);
    vecpart1=vec;
    vecpart2=vec;
    
    % Plot the data:
    subplot(2,2,i)
    hold on
    
    % If the host population goes extinct:
    extinction_boundary1=find(totalpopvec1<1e-4,1);
    if ~isempty(extinction_boundary1)
        plot(log10(vec(extinction_boundary1))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log10(vec(extinction_boundary1))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart1(extinction_boundary1+1:end,:)=[];
        Aplotvecpart1(extinction_boundary1+1:end,:)=[];
        disprevvec1(extinction_boundary1+1:end,:)=[];
        vecpart1(extinction_boundary1+1:end)=[];
    end
    % If the second host population (in a bistable case) goes extinct:
    extinction_boundary2=find(totalpopvec2<1e-4,1);
    if ~isempty(extinction_boundary2)
        plot(log10(vec(extinction_boundary2))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log10(vec(extinction_boundary2))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart2(extinction_boundary2+1:end,:)=[];
        Aplotvecpart2(extinction_boundary2+1:end,:)=[];
        disprevvec2(extinction_boundary2+1:end,:)=[];
        vecpart2(extinction_boundary2+1:end)=[];
    end
    
    % The left-hand axis will show the level of resistance at the singular
    % strategies:
    yyaxis left
    plot(log10(vecpart1),Jplotvecpart1,'r-','linewidth',3)
    plot(log10(vecpart2),Jplotvecpart2,'r-','linewidth',3)
    plot(log10(vecpart1),Aplotvecpart1,'b--','linewidth',3)
    plot(log10(vecpart2),Aplotvecpart2,'b--','linewidth',3)
    ax=gca;
    ax.YColor='k';
    ylim([0,1])
    
    % The right-hand axis will show the total population density and the
    % density of infected hosts.
    yyaxis right
    plot(log10(vecpart1),totalpopvec1,':','linewidth',3,'color',[0.8,0.8,0.8])
    plot(log10(vecpart2),totalpopvec2,':','linewidth',3,'color',[0.8,0.8,0.8])
    plot(log10(vecpart1),diseasedensity1,'-','linewidth',3,'color',[0.95,0.95,0.95])
    plot(log10(vecpart2),diseasedensity2,'-','linewidth',3,'color',[0.95,0.95,0.95])
    set(gca,'xtick',0:1:6,'xticklabel',10.^(0:1:6));
    set(gca,'ytick',0:0.2:1);
    ax=gca;
    ax.YColor=[0.5,0.5,0.5];
    ylim([0,1])
    box on

    % We want the resistance curves to be on top:
    yyaxis right
    plot(log10(vecpart1),Jplotvecpart1,'r-','linewidth',3)
    plot(log10(vecpart2),Jplotvecpart2,'r-','linewidth',3)
    plot(log10(vecpart1),Aplotvecpart1,'b--','linewidth',3)
    plot(log10(vecpart2),Aplotvecpart2,'b--','linewidth',3)
    
    % Determine and label any regions where bistability occurs:
    bistability_boundaryJ=find(~isnan(Jplotvec(:,1).*Jplotvec(:,2)),1);
    bistability_boundaryA=find(~isnan(Aplotvec(:,1).*Aplotvec(:,2)),1);
    bistability_boundaryJ2=find(~isnan(Jplotvec(:,1).*Jplotvec(:,3)),1);
    bistability_boundaryA2=find(~isnan(Aplotvec(:,1).*Aplotvec(:,3)),1);
    if ~isempty(bistability_boundaryJ) && ~isempty(bistability_boundaryJ2)
        bistability_boundaryJ=min(bistability_boundaryJ,bistability_boundaryJ2);
    end
    if ~isempty(bistability_boundaryA) && ~isempty(bistability_boundaryA2)
        bistability_boundaryA=min(bistability_boundaryA,bistability_boundaryA2);
    end
    if isempty(bistability_boundaryJ) && ~isempty(bistability_boundaryJ2)
        bistability_boundaryJ=bistability_boundaryJ2;
    end
    if isempty(bistability_boundaryA) && ~isempty(bistability_boundaryA2)
        bistability_boundaryA=bistability_boundaryA2;
    end
    bistability_boundaryJ3=find(~isnan(Jplotvec(:,2).*Jplotvec(:,3)),1);
    bistability_boundaryA3=find(~isnan(Aplotvec(:,2).*Aplotvec(:,3)),1);
    if ~isempty(bistability_boundaryJ) && ~isempty(bistability_boundaryJ3)
        bistability_boundaryJ=min(bistability_boundaryJ,bistability_boundaryJ3);
    end
    if ~isempty(bistability_boundaryA) && ~isempty(bistability_boundaryA3)
        bistability_boundaryA=min(bistability_boundaryA,bistability_boundaryA3);
    end
    if isempty(bistability_boundaryJ) && ~isempty(bistability_boundaryJ3)
        bistability_boundaryJ=bistability_boundaryJ3;
    end
    if isempty(bistability_boundaryA) && ~isempty(bistability_boundaryA3)
        bistability_boundaryA=bistability_boundaryA3;
    end
    if ~isempty(bistability_boundaryJ) && ~isempty(bistability_boundaryA)
        bistability_boundary=min(bistability_boundaryJ,bistability_boundaryA);
    elseif ~isempty(bistability_boundaryJ)
        bistability_boundary=bistability_boundaryJ;
    elseif ~isempty(bistability_boundaryA)
        bistability_boundary=bistability_boundaryA;
    else
        bistability_boundary=[];
    end
    if ~isempty(bistability_boundary)
        plot(log10(vec(bistability_boundary))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log10(vec(bistability_boundary))+0.1,0.4,'bistability','fontsize',16);
        set(mytext,'Rotation',90)
    end
    xlim([0,upperx])
    
    % Label the plots:
    if i==1
        text(0.02,0.93,'A','fontsize',20);
        yyaxis left
        ylab=ylabel('Resistance','interpreter','latex','fontsize',18);
        temp = get(ylab,'position');
        temp(2) = temp(2) - 0.7;
        temp(1) = temp(1) ;
        set(ylab,'position',temp);
    elseif i==2
        text(0.02,0.93,'B','fontsize',20);
        yyaxis right
        ylab=ylabel('Population density','interpreter','latex','fontsize',18);
        temp = get(ylab,'position');
        temp(2) = temp(2) - 0.7;
        temp(1) = temp(1) ;
        set(ylab,'position',temp);
    elseif i==3
        text(0.02,0.93,'C','fontsize',20);
        xlab=xlabel('Baseline transmissibility, $\beta_0$','interpreter','latex','fontsize',16);
        temp = get(xlab,'position');
        temp(1) = temp(1) + 0.7;
        temp(2) = temp(2) ;
        set(xlab,'position',temp);
    elseif i==4
        text(0.02,0.93,'D','fontsize',20);
    end
    
end