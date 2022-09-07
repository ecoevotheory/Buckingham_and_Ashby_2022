function [resJssvalue,resAssvalue] = singstrat_finder_analytical(resJlowerbound,resJupperbound,resAlowerbound,resAupperbound,c1g,c2g,c1a,c2a,a0,g0,beta0,h,f,alpha,version)

% This function finds singular strategies analytically.

% Set up variables, parameters and vectors:
syms SJ SA IJ IA resJ resA resJm resAn bJ bA
eps=1e-5;
resJssvalue=[];
resAssvalue=[];
betaJ=beta0*(1-resJ);
betaA=beta0*(1-resA);

% Define the trade-offs and invasion fitness for each version:
if version==1
    a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)));
    g=g0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
    am=subs(a,resA,resAn);
    gm=subs(g,resJ,resJm);
    bJ=1;
    bA=1;
    
    lambdaJstar=beta0*(1-resJ)*(h*IJ+IA);
    lambdaJstarm=beta0*(1-resJm)*(h*IJ+IA);
    lambdaAstar=beta0*(1-resA)*(h*IJ+IA);
    lambdaAstarn=beta0*(1-resAn)*(h*IJ+IA);
    
    Am=(bA+alpha)*(bJ+gm+alpha)+f*(bJ+gm+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+gm+alpha)*(bA+lambdaAstar)*(bJ+gm+lambdaJstarm);
    An=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bA+lambdaAstarn);
    Bn=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstarn)*(bJ+g+lambdaJstar);
    
    wJ=((gm*a*(1-(SJ+SA+IJ+IA))*Am)/Bm) -1;
    wA=((g*am*(1-(SJ+SA+IJ+IA))*An)/Bn) -1;

elseif version==2
    bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
    bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
    bAm=subs(bA,resA,resAn);
    bJm=subs(bJ,resJ,resJm);
    a=a0;
    g=g0;
    
    lambdaJstar=beta0*(1-resJ)*(h*IJ+IA);
    lambdaJstarm=beta0*(1-resJm)*(h*IJ+IA);
    lambdaAstar=beta0*(1-resA)*(h*IJ+IA);
    lambdaAstarn=beta0*(1-resAn)*(h*IJ+IA);
    
    Am=(bA+alpha)*(bJm+g+alpha)+f*(bJm+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJm+g+alpha)*(bA+lambdaAstar)*(bJm+g+lambdaJstarm);
    An=(bAm+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bAm+lambdaAstarn);
    Bn=(bAm+alpha)*(bJ+g+alpha)*(bAm+lambdaAstarn)*(bJ+g+lambdaJstar);
    
    wJ=((g*a*(1-(SJ+SA+IJ+IA))*Am)/Bm) -1;
    wA=((g*a*(1-(SJ+SA+IJ+IA))*An)/Bn) -1;

elseif version==3
    bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
    a=a0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
    bAm=subs(bA,resA,resAn);
    am=subs(a,resJ,resJm);
    g=g0;
    bJ=1;
    
    lambdaJstar=beta0*(1-resJ)*(h*IJ+IA);
    lambdaJstarm=beta0*(1-resJm)*(h*IJ+IA);
    lambdaAstar=beta0*(1-resA)*(h*IJ+IA);
    lambdaAstarn=beta0*(1-resAn)*(h*IJ+IA);
    
    Am=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstar)*(bJ+g+lambdaJstarm);
    An=(bAm+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bAm+lambdaAstarn);
    Bn=(bAm+alpha)*(bJ+g+alpha)*(bAm+lambdaAstarn)*(bJ+g+lambdaJstar);
    
    wJ=((g*am*(1-(SJ+SA+IJ+IA))*Am)/Bm) -1;
    wA=((g*a*(1-(SJ+SA+IJ+IA))*An)/Bn) -1;

elseif version==4
    g=g0*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
    bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
    bAm=subs(bA,resJ,resJm);
    gm=subs(g,resA,resAn);
    a=a0;
    bJ=1;
    
    lambdaJstar=beta0*(1-resJ)*(h*IJ+IA);
    lambdaJstarm=beta0*(1-resJm)*(h*IJ+IA);
    lambdaAstar=beta0*(1-resA)*(h*IJ+IA);
    lambdaAstarn=beta0*(1-resAn)*(h*IJ+IA);
    
    Am=(bA+alpha)*(bJ+gm+alpha)+f*(bJ+gm+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+gm+alpha)*(bA+lambdaAstar)*(bJ+gm+lambdaJstarm);
    An=(bAm+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bAm+lambdaAstarn);
    Bn=(bAm+alpha)*(bJ+g+alpha)*(bAm+lambdaAstarn)*(bJ+g+lambdaJstar);
    
    wJ=((g*a*(1-(SJ+SA+IJ+IA))*Am)/Bm) -1;
    wA=((gm*a*(1-(SJ+SA+IJ+IA))*An)/Bn) -1;

elseif version==5
    bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
    a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)));
    am=subs(a,resA,resAn);
    bJm=subs(bJ,resJ,resJm);
    bA=1;
    g=g0;
    
    lambdaJstar=beta0*(1-resJ)*(h*IJ+IA);
    lambdaJstarm=beta0*(1-resJm)*(h*IJ+IA);
    lambdaAstar=beta0*(1-resA)*(h*IJ+IA);
    lambdaAstarn=beta0*(1-resAn)*(h*IJ+IA);
    
    Am=(bA+alpha)*(bJm+g+alpha)+f*(bJm+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJm+g+alpha)*(bA+lambdaAstar)*(bJm+g+lambdaJstarm);
    An=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bA+lambdaAstarn);
    Bn=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstarn)*(bJ+g+lambdaJstar);
    
    wJ=((g*a*(1-(SJ+SA+IJ+IA))*Am)/Bm) -1;
    wA=((g*am*(1-(SJ+SA+IJ+IA))*An)/Bn) -1;

elseif version==6
    a=a0*(1-c1a*(1-exp(-c2a*resA))/(1-exp(-c2a)))*(1-c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g)));
    am=subs(a,resJ,resJm);
    an=subs(a,resA,resAn);
    g=g0;
    bJ=1;
    bA=1;
    
    lambdaJstar=beta0*(1-resJ)*(h*IJ+IA);
    lambdaJstarm=beta0*(1-resJm)*(h*IJ+IA);
    lambdaAstar=beta0*(1-resA)*(h*IJ+IA);
    lambdaAstarn=beta0*(1-resAn)*(h*IJ+IA);
    
    Am=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstar)*(bJ+g+lambdaJstarm);
    An=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bA+lambdaAstarn);
    Bn=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstarn)*(bJ+g+lambdaJstar);
    
    wJ=((g*am*(1-(SJ+SA+IJ+IA))*Am)/Bm) -1;
    wA=((g*an*(1-(SJ+SA+IJ+IA))*An)/Bn) -1;

end
            
% State the ecological system:
eq1=a*(1-SJ-SA-IJ-IA)*(SA+f*IA) -(bJ+g+betaJ*(h*IJ+IA))*SJ;
eq2=g*SJ-(bA+betaA*(h*IJ+IA))*SA;
eq3=betaJ*(h*IJ+IA)*SJ-(bJ+g+alpha)*IJ;
eq4=g*IJ+betaA*(h*IJ+IA)*SA-(bA+alpha)*IA;

% Determine the fitness gradient:
fitgradJworking(resA,resJm,IJ,IA,SJ,SA)=diff(wJ,resJm);
fitgradAworking(resJ,resAn,IJ,IA,SJ,SA)=diff(wA,resAn);
fitgradJagain(resA,resJ,IJ,IA,SJ,SA)=subs(fitgradJworking,resJm,resJ);
fitgradAagain(resA,resJ,IJ,IA,SJ,SA)=subs(fitgradAworking,resAn,resA);

% Find the roots of the fitness gradients and the ecological system:
overallsol=vpasolve([eq1,eq2,eq3,eq4,fitgradJagain,fitgradAagain],[resJ,resA,SJ,SA,IJ,IA],[resJlowerbound resJupperbound;resAlowerbound resAupperbound; eps inf; eps inf; eps inf; eps inf]);

% This gives the singular strategy:
if ~isempty(overallsol.resJ)
    if ~isempty(overallsol.resA)
        resJssvalue=overallsol.resJ;
        resAssvalue=overallsol.resA;
    end
end