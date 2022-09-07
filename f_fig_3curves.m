% This code creates six sub-figures showing the effect of varying the 
% sterility virulence, 1-f, for different combinations of trade-offs. It 
% plots three curves (showing the effect for three different values of 
% another varying parameter). 

% The solid curve shows the juvenile resistance and the dashed curve shows 
% the adult resistance. Orange denotes a high value of the parameter being
% varied between the curves, grey denotes an intermediate value and purple
% denotes a low value. 

% Parameter values:
a0vec=[4,5,10];
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
beta0=8;
alpha=0;

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
    
    subplot(2,3,i)
    hold on
    
    % Intermediate parameter set:
    a0=a0vec(2);
    % 'f_fig_data' is a function which finds the singular strategies for
    % different values of sterility virulence. 
    [vec,Jplotvec,Aplotvec,~,~]=f_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,alpha,beta0);
    xplotvec=1-vec;
    % Make the plot:
    plot(xplotvec,Jplotvec,'-','linewidth',3,'color',[0.6,0.6,0.6])
    plot(xplotvec,Aplotvec,'--','linewidth',3,'color',[0.6,0.6,0.6])
    ylim([0,1])
    set(gca,'xtick',0:0.2:1);
    set(gca,'ytick',0:0.2:1);
    box on
    
    % High parameter set:
    a0=a0vec(3);
    [vec,Jplotvec,Aplotvec,~,~]=f_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,alpha,beta0);
    xplotvec=1-vec;
    plot(xplotvec,Jplotvec,'r-','linewidth',3)
    plot(xplotvec,Aplotvec,'r--','linewidth',3)
    
    % Low parameter set:
    a0=a0vec(1);
    [vec,Jplotvec,Aplotvec,~,~]=f_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,alpha,beta0);
    xplotvec=1-vec;
    plot(xplotvec,Jplotvec,'b-','linewidth',3)
    plot(xplotvec,Aplotvec,'b--','linewidth',3)
    
end

% Label the figure:
for i=1:6
    subplot(2,3,i)
    xlim([0,1])
    text(0.02,0.93,strcat('\bf{',labs{i},'}'),'interpreter','latex','fontsize',20);
    
    title(str(ORDER(i)),'interpreter','latex','fontsize',12);
    set(gca,'fontsize',10);
    if(i==4)
        ylab = ylabel('Resistance','interpreter','latex','fontsize',18);
        temp = get(ylab,'position');
        temp(2) = temp(2) + 0.5;
        temp(1) = temp(1) ;
        set(ylab,'position',temp);
    elseif(i==5)
        xlabel('Sterility virulence, $1-f$','interpreter','latex','fontsize',16)
    elseif i==1
        xlabel('Effect of baseline reproduction, $a_0$','interpreter','latex','fontsize',16)
    end
    
end
