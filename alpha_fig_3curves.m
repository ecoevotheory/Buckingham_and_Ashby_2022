% This code creates six sub-figures showing the effect of varying the 
% mortality virulence, alpha, for different combinations of 
% trade-offs. It plots three curves (showing the effect for three different
% values of another varying parameter). 

% The solid curve shows the juvenile resistance and the dashed curve shows 
% the adult resistance. Orange denotes a high value of the parameter being
% varied between the curves, grey denotes an intermediate value and purple
% denotes a low value. 

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

% Parameter values:
beta0=8;
g0=1;
c1a=0.5;
c2a=-3;
c1g=0.5;
c2g=-3;
f=1;
a0vec=[4,5,10];

for i=1:6
    version = ORDER(i);
    
    subplot(2,3,i)
    hold on
    
    % Find singular strategies given an intermediate value of a0:
    a0=a0vec(2);
    [vec,Jplotvec,Aplotvec,~,totalpopvec1,~,totalpopvec2]=alpha_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,f,beta0);
    upperx=max(log2(vec));
    lowerx=min(log2(vec));
    Jplotvecpart1=Jplotvec(:,1);
    Jplotvecpart2=Jplotvec(:,2);
    Aplotvecpart1=Aplotvec(:,1);
    Aplotvecpart2=Aplotvec(:,2);
    vecpart1=vec;
    vecpart2=vec;
    % If the host population goes extinct:
    extinction_boundary1=find(totalpopvec1<1e-4,1);
    if ~isempty(extinction_boundary1)
        plot(log2(vec(extinction_boundary1))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log2(vec(extinction_boundary1))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart1(extinction_boundary1+1:end,:)=[];
        Aplotvecpart1(extinction_boundary1+1:end,:)=[];
        vecpart1(extinction_boundary1+1:end)=[];
    end
    % If the second host population (in a bistable case) goes extinct:
    extinction_boundary2=find(totalpopvec2<1e-4,1);
    if ~isempty(extinction_boundary2)
        plot(log2(vec(extinction_boundary2))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log2(vec(extinction_boundary2))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart2(extinction_boundary2+1:end,:)=[];
        Aplotvecpart2(extinction_boundary2+1:end,:)=[];
        vecpart2(extinction_boundary2+1:end)=[];
    end
    % Make the plot:
    plot(log2(vecpart1),Jplotvecpart1,'-','linewidth',3,'color',[0.6,0.6,0.6])
    plot(log2(vecpart2),Jplotvecpart2,'-','linewidth',3,'color',[0.6,0.6,0.6])
    plot(log2(vecpart1),Aplotvecpart1,'--','linewidth',3,'color',[0.6,0.6,0.6])
    plot(log2(vecpart2),Aplotvecpart2,'--','linewidth',3,'color',[0.6,0.6,0.6])
    set(gca,'xtick',-3:1:4,'xticklabel',{'1/8','1/4','1/2','1','2','4','8','16'});
    set(gca,'ytick',0:0.2:1);
    box on
    ylim([0,1])
    
    
    % Find singular strategies given a high value of a0:
    a0=a0vec(3);
    [vec,Jplotvec,Aplotvec,~,totalpopvec1,~,totalpopvec2]=alpha_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,f,beta0);
    Jplotvecpart1=Jplotvec(:,1);
    Jplotvecpart2=Jplotvec(:,2);
    Aplotvecpart1=Aplotvec(:,1);
    Aplotvecpart2=Aplotvec(:,2);
    vecpart1=vec;
    vecpart2=vec;
    % If the host population goes extinct:
    extinction_boundary1=find(totalpopvec1<1e-4,1);
    if ~isempty(extinction_boundary1)
        plot(log2(vec(extinction_boundary1))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log2(vec(extinction_boundary1))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart1(extinction_boundary1+1:end,:)=[];
        Aplotvecpart1(extinction_boundary1+1:end,:)=[];
        vecpart1(extinction_boundary1+1:end)=[];
    end
    % If the second host population (in a bistable case) goes extinct:
    extinction_boundary2=find(totalpopvec2<1e-4,1);
    if ~isempty(extinction_boundary2)
        plot(log2(vec(extinction_boundary2))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log2(vec(extinction_boundary2))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart2(extinction_boundary2+1:end,:)=[];
        Aplotvecpart2(extinction_boundary2+1:end,:)=[];
        vecpart2(extinction_boundary2+1:end)=[];
    end
    % Make the plot:
    plot(log2(vecpart1),Jplotvecpart1,'r-','linewidth',3)
    plot(log2(vecpart2),Jplotvecpart2,'r-','linewidth',3)
    plot(log2(vecpart1),Aplotvecpart1,'r--','linewidth',3)
    plot(log2(vecpart2),Aplotvecpart2,'r--','linewidth',3)
    
    % Find singular strategies given a low value of a0:
    a0=a0vec(1);
    [vec,Jplotvec,Aplotvec,~,totalpopvec1,~,totalpopvec2]=alpha_fig_data(version,a0,g0,c1a,c2a,c1g,c2g,f,beta0);
    Jplotvecpart1=Jplotvec(:,1);
    Jplotvecpart2=Jplotvec(:,2);
    Aplotvecpart1=Aplotvec(:,1);
    Aplotvecpart2=Aplotvec(:,2);
    vecpart1=vec;
    vecpart2=vec;
    % If the host population goes extinct:
    extinction_boundary1=find(totalpopvec1<1e-4,1);
    if ~isempty(extinction_boundary1)
        plot(log2(vec(extinction_boundary1))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log2(vec(extinction_boundary1))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart1(extinction_boundary1+1:end,:)=[];
        Aplotvecpart1(extinction_boundary1+1:end,:)=[];
        vecpart1(extinction_boundary1+1:end)=[];
    end
    % If the second host population (in a bistable case) goes extinct:
    extinction_boundary2=find(totalpopvec2<1e-4,1);
    if ~isempty(extinction_boundary2)
        plot(log2(vec(extinction_boundary2))+zeros(101,1),0:0.01:1,'k--','linewidth',1.5,'color',[0.5,0.5,0.6]);
        mytext=text(log2(vec(extinction_boundary2))+0.2,0.4,'host extinction');
        set(mytext,'Rotation',90)
        Jplotvecpart2(extinction_boundary2+1:end,:)=[];
        Aplotvecpart2(extinction_boundary2+1:end,:)=[];
        vecpart2(extinction_boundary2+1:end)=[];
    end
    % Make the plot:
    plot(log2(vecpart1),Jplotvecpart1,'b-','linewidth',3)
    plot(log2(vecpart2),Jplotvecpart2,'b-','linewidth',3)
    plot(log2(vecpart1),Aplotvecpart1,'b--','linewidth',3)
    plot(log2(vecpart2),Aplotvecpart2,'b--','linewidth',3)
    
end

% Label the figure:
for i=1:6
    subplot(2,3,i)
    xlim([lowerx,upperx])
    text(-3.2,0.93,strcat('\bf{',labs{i},'}'),'interpreter','latex','fontsize',20);
    
    title(str(ORDER(i)),'interpreter','latex','fontsize',12);
    set(gca,'fontsize',10);
    if(i==4)
        ylab = ylabel('Resistance','interpreter','latex','fontsize',18);
        temp = get(ylab,'position');
        temp(2) = temp(2) + 0.5;
        temp(1) = temp(1) ;
        set(ylab,'position',temp);
    elseif(i==5)
        xlabel('Mortality virulence, $\alpha$','interpreter','latex','fontsize',16)
    elseif i==1
        xlabel('Effect of baseline reproduction, $a_0$','interpreter','latex','fontsize',16)
    end
    
end
