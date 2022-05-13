function effect_of_f_figure

% This code creates six sub-figures showing the effect of varying the 
% sterility virulence, 1-f, for different combinations of trade-offs.

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

str{1}='J: maturation, A: reproduction';
str{2}='J: juvenile mortality, A: adult mortality';
str{3}='J: reproduction, A: adult mortality';
str{4}='J: maturation, A: adult mortality';
str{5}='J: juvenile mortality, A: reproduction';
str{6}='J: reproduction, A: reproduction';

% The plotting order is defined so that the top row contains the trade-offs 
% between adult resistance and reproduction and the bottom row contains the
% trade-offs between adult resistance and adult mortality.
ORDER = [1,5,6,4,2,3];
labs = {'A','B','C','D','E','F'};

for i=1:6
    version = ORDER(i);
    
    % The function 'f_fig_data' determines the singular strategies and
    % disease prevalence for a variety of values of the sterility
    % virulence, given certain parameters which are defined inside the
    % function.
    [vec,Jplotvec,Aplotvec,~,~,disprevvec,~,~]=f_fig_data(version);
    
    % Plot the data:
    subplot(2,3,i)
    hold on
    
    % The left-hand axis will show the level of resistance at the singular
    % strategies:
    yyaxis left
    plot(1-vec,Jplotvec,'r-','linewidth',3)
    plot(1-vec,Aplotvec,'b--','linewidth',3)
    ax=gca;
    ax.YColor='k';
    ylim([0,1])
    
    % The right-hand axis will show the proportion of hosts which are
    % infected:
    yyaxis right
    plot(1-vec,disprevvec,':','linewidth',3,'color',[0.8,0.8,0.8])
    set(gca,'xtick',0:0.2:1);
    set(gca,'ytick',0:0.2:1);
    ax=gca;
    ax.YColor=[0.5,0.5,0.5];
    ylim([0,1])
    box on

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
        plot(1-(vec(bistability_boundary))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(1-(vec(bistability_boundary))+0.2,0.4,'bistability');
        set(mytext,'Rotation',90)
    end
    
end

% Label the figure:
for i=1:6
    subplot(2,3,i)
    xlim([0,1])
    text(0.02,0.93,strcat('\bf{',labs{i},'}'),'interpreter','latex','fontsize',20);
    
    title(str(ORDER(i)),'interpreter','latex','fontsize',12);
    set(gca,'fontsize',10);
    if(i==4)
        yyaxis left
        ylab = ylabel('Resistance','interpreter','latex','fontsize',18);
        temp = get(ylab,'position');
        temp(2) = temp(2) + 0.5;
        temp(1) = temp(1) ;
        set(ylab,'position',temp);
    elseif (i==6)
        yyaxis right
        ylab=ylabel('Proportion infected','interpreter','latex','fontsize',18);
        temp = get(ylab,'position');
        temp(2) = temp(2) + 0.7;
        set(ylab,'position',temp);
    elseif(i==5)
        xlabel('Fecundity virulence, $1-f$','interpreter','latex','fontsize',16)
    end
    
end

end
