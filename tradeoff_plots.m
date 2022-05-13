% This code sketches different trade-off functions and loads the model 
% schematic.

%% Trade-off between adult resistance and birth rate

% Set up variables and parameters:
syms a(resA)
a0=5;

% Here are some example parameter values in the case of a strong, high
% curvature trade-off:
c1a=0.5;
c2a=3;

% Define the trade-off function:
a(resA)=a0*(1-c1a*(1-exp(c2a*resA))/(1-exp(c2a)));

% Create vectors of the values of the trade-off function:
xvec=zeros(100,1);
avec1=zeros(100,1);
for i=1:100
    xvec(i)=i/100;
    avec1(i)=a(xvec(i));
end

% Now do the same for a weak, low curvature trade-off:
c1a=0.25;
c2a=3;
a(resA)=a0*(1-c1a*(1-exp(c2a*resA))/(1-exp(c2a)));
avec2=zeros(100,1);
for i=1:100
    avec2(i)=a(xvec(i));
end

% Now do the same for a strong, high curvature trade-off:
c1a=0.5;
c2a=10;
a(resA)=a0*(1-c1a*(1-exp(c2a*resA))/(1-exp(c2a)));
avec3=zeros(100,1);
for i=1:100
    avec3(i)=a(xvec(i));
end

% Now do the same for a weak, high curvature trade-off:
c1a=0.25;
c2a=10;
a(resA)=a0*(1-c1a*(1-exp(c2a*resA))/(1-exp(c2a)));
avec4=zeros(100,1);
for i=1:100
    avec4(i)=a(xvec(i));
end

% Plot the trade-offs:
firstcolour=1/255*[215,48,39];
secondcolour=1/255*[69,117,180];
subplot(3,2,2)
plot(xvec,avec1,'Color',firstcolour,'Linewidth',1.5)
hold on
plot(xvec,avec2,'Color',secondcolour,'Linewidth',1.5)
hold on
plot(xvec,avec3,'Color',firstcolour,'LineStyle','--','Linewidth',1.5) 
hold on
plot(xvec,avec4,'Color',secondcolour,'LineStyle','--','Linewidth',1.5)
title('(b) (i)','position',[0.05,5.1],'FontSize',14)
xlabel('Adult resistance, $r_A$','interpreter','latex','fontsize',16)
ylabel('Birth rate, $a(r_A)$','interpreter','latex','fontsize',16)
ylim([2.5,5])

%%  Trade-off between adult resistance and adult natural death rate

% Set up variables and parameters:
syms bA(resA)

% Here are some example parameter values in the case of a strong, high
% curvature trade-off:
c1a=0.5;
c2a=3;

% Define the trade-off function:
bA(resA)=1+c1a*(1-exp(c2a*resA))/(1-exp(c2a));

% Create vectors of the values of the trade-off function:
xvec=zeros(100,1);
bAvec1=zeros(100,1);
for i=1:100
    xvec(i)=i/100;
    bAvec1(i)=bA(xvec(i));
end

% Now do the same for a weak, low curvature trade-off:
c1a=0.25;
c2a=3;
bA(resA)=1+c1a*(1-exp(c2a*resA))/(1-exp(c2a));
bAvec2=zeros(100,1);
for i=1:100
    bAvec2(i)=bA(xvec(i));
end

% Now do the same for a strong, high curvature trade-off:
c1a=0.5;
c2a=10;
bA(resA)=1+c1a*(1-exp(c2a*resA))/(1-exp(c2a));
bAvec3=zeros(100,1);
for i=1:100
    bAvec3(i)=bA(xvec(i));
end

% Now do the same for a weak, high curvature trade-off:
c1a=0.25;
c2a=10;
bA(resA)=1+c1a*(1-exp(c2a*resA))/(1-exp(c2a));
bAvec4=zeros(100,1);
for i=1:100
    bAvec4(i)=bA(xvec(i));
end

% Plot the trade-offs:
firstcolour=1/255*[215,48,39];
secondcolour=1/255*[69,117,180];
thirdcolour=1/255*[217,95,2];
fourthcolour=1/255*[27,158,119];
subplot(3,2,4)
plot(xvec,bAvec1,'Color',firstcolour,'Linewidth',1.5)
hold on
plot(xvec,bAvec2,'Color',secondcolour,'Linewidth',1.5)
hold on
plot(xvec,bAvec3,'Color',firstcolour,'LineStyle','--','Linewidth',1.5) 
hold on
plot(xvec,bAvec4,'Color',secondcolour,'LineStyle','--','Linewidth',1.5)
title('(b) (ii)','position',[0.05,1.52],'FontSize',14)
xlabel('Adult resistance, $r_A$','interpreter','latex','fontsize',16)
ylabel({'Adult natural'; 'death rate, $b_A(r_A)$'},'interpreter','latex','fontsize',16)
ylim([1,1.5])

%% Both resistance traits with reproduction

% Set up variables and parameters:
syms a(resJ,resA)
a0=5;

% Here are some example parameter values in the case of a strong, high
% curvature trade-off:
c1a=0.5;
c2a=3;
c1g=0.5;
c2g=3;

% Define the trade-off function:
a(resJ,resA)=a0*(1-c1a*(1-exp(c2a*resA))/(1-exp(c2a)))*(1-c1g*(1-exp(c2g*resJ))/(1-exp(c2g)));

% Create vectors of the values of the trade-off function:
resJvec=zeros(100,1);
resAvec=zeros(100,1);
avec=zeros(100,100);
for i=1:100
    for j=1:100
        resJvec(i)=i/100;
        resAvec(j)=j/100;
        avec(i,j)=a(resJvec(i),resAvec(j));
    end
end

% Plot the trade-off:
subplot(3,2,6)
imagesc(resJvec,resAvec,avec)
colormap(flipud(hot))
c=colorbar;
set(gca,'YDir','normal')
xlabel('Juvenile resistance, $r_J$','interpreter','latex','fontsize',16)
ylabel('Adult resistance, $r_A$','interpreter','latex','fontsize',16)
xlim([0,1])
ylim([0,1])
title('(b) (iii)','position',[0.07,1.05],'FontSize',14)
c.Label.Interpreter='latex';
c.Label.String='Birth rate, $$a(r_J,r_A)$$';
c.Label.Position(1)=2;
set(c,'ylim',[1,5])

%% Include model schematic image

s=subplot(3,2,[1,3,5]);
A=imread('Model Schematic.jpg');
B=imresize(A,2);
image(B)
title('(a)','FontSize',14)
set(gca,'xtick',[])
set(gca,'ytick',[])
s.Position=s.Position+[-0.01,0.05,0.01,-0.05];
