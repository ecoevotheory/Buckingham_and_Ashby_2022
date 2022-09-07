% This code determines analytical expressions for the fitness gradient and
% the second derivatives of the invasion fitness, for different
% combinations of trade-offs. These expressions are used directly 
% throughout the code for this model.

% State which version of the model you want the expressions to hold for:
version=1;

% For each versions, we produce the expressions:
if version==1

    % Define variables and trade-offs:
    syms SJval SAval IJval IAval resJ resA resJm resAn a alpha g g0 h a0 f c1a c1g c2a c2g beta0 resJval resAval
    a=a0*(1-(c1a*(1-exp(-c2a*resA))/(1-exp(-c2a))));
    an=subs(a,resA,resAn);
    g=g0*(1-(c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g))));
    gm=subs(g,resJ,resJm);
    bJ=1;
    bA=1;
    
    % Set up parameters and functions:
    lambdaJstar=beta0*(1-resJ)*(h*IJval+IAval);
    lambdaJstarm=beta0*(1-resJm)*(h*IJval+IAval);
    lambdaAstar=beta0*(1-resA)*(h*IJval+IAval);
    lambdaAstarn=beta0*(1-resAn)*(h*IJval+IAval);
    Am=(bA+alpha)*(bJ+gm+alpha)+f*(bJ+gm+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+gm+alpha)*(bA+lambdaAstar)*(bJ+gm+lambdaJstarm);
    An=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bA+lambdaAstarn);
    Bn=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstarn)*(bJ+g+lambdaJstar);
    % Define the invasion fitness:
    wJ=((gm*a*(1-(SJval+SAval+IJval+IAval))*Am)/Bm) -1;
    wA=((g*an*(1-(SJval+SAval+IJval+IAval))*An)/Bn) -1;
    % Determine the fitness gradient:
    fitgradJworking=diff(wJ,resJm);
    fitgradAworking=diff(wA,resAn);
    fitgradJfinal=subs(fitgradJworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    fitgradAfinal=subs(fitgradAworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine the second derivatives with respect to the mutant trait
    % only:
    deriv2Jworking=diff(fitgradJworking,resJm);
    deriv2J=subs(deriv2Jworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    deriv2Aworking=diff(fitgradAworking,resAn);
    deriv2A=subs(deriv2Aworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    
    % Determine the cross-derivatives at (resJval,resAval) and at points
    % nearby:
    syms SJupJ SAupJ IJupJ IAupJ SJdownJ SAdownJ IJdownJ IAdownJ eps SJupA SAupA IJupA IAupA SJdownA SAdownA IJdownA IAdownA
    cderivJ1working=(subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivJ2working=(subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivJ1=subs(cderivJ1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivJ2=subs(cderivJ2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA1working=(subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivA2working=(subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivA1=subs(cderivA1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA2=subs(cderivA2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);

elseif version==2

    % Define variables and trade-offs:
    syms SJval SAval IJval IAval resJ resA resJm resAn g0 alpha h a0 f c1a c1g c2a c2g beta0 bJ bA bJm bAm resJval resAval
    a=a0;
    g=g0;
    bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
    bJm=subs(bJ,resJ,resJm);
    bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
    bAm=subs(bA,resA,resAn);
    
    % Set up parameters and functions:
    lambdaJstar=beta0*(1-resJ)*(h*IJval+IAval);
    lambdaJstarm=beta0*(1-resJm)*(h*IJval+IAval);
    lambdaAstar=beta0*(1-resA)*(h*IJval+IAval);
    lambdaAstarn=beta0*(1-resAn)*(h*IJval+IAval);
    Am=(bA+alpha)*(bJm+g+alpha)+f*(bJm+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJm+g+alpha)*(bA+lambdaAstar)*(bJm+g+lambdaJstarm);
    An=(bAm+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bAm+lambdaAstarn);
    Bn=(bAm+alpha)*(bJ+g+alpha)*(bAm+lambdaAstarn)*(bJ+g+lambdaJstar);
    % Define the invasion fitness:
    wJ=((g*a*(1-(SJval+SAval+IJval+IAval))*Am)/Bm) -1;
    wA=((g*a*(1-(SJval+SAval+IJval+IAval))*An)/Bn) -1;
    % Determine the fitness gradient:
    fitgradJworking=diff(wJ,resJm);
    fitgradAworking=diff(wA,resAn);
    fitgradJfinal=subs(fitgradJworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    fitgradAfinal=subs(fitgradAworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine second derivatives with respect to mutant traits only:
    deriv2Jworking=diff(fitgradJworking,resJm);
    deriv2J=subs(deriv2Jworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    deriv2Aworking=diff(fitgradAworking,resAn);
    deriv2A=subs(deriv2Aworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine cross-derivatives at (resJval, resAval) and at points
    % nearby:
    syms SJupJ SAupJ IJupJ IAupJ SJdownJ SAdownJ IJdownJ IAdownJ eps SJupA SAupA IJupA IAupA SJdownA SAdownA IJdownA IAdownA
    cderivJ1working=(subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivJ2working=(subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivJ1=subs(cderivJ1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivJ2=subs(cderivJ2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA1working=(subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivA2working=(subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivA1=subs(cderivA1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA2=subs(cderivA2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);

elseif version==3

    % Define variables and trade-offs:
    syms SJval SAval IJval IAval resJ resA resJm resAn alpha a g0 h a0 f c1a c1g c2a c2g beta0 bA bAm resJval resAval
    g=g0;
    bJ=1;
    bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
    bAm=subs(bA,resA,resAn);
    a=a0*(1-(c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g))));
    am=subs(a,resJ,resJm);
    
    % Set up parameters and functions:
    lambdaJstar=beta0*(1-resJ)*(h*IJval+IAval);
    lambdaJstarm=beta0*(1-resJm)*(h*IJval+IAval);
    lambdaAstar=beta0*(1-resA)*(h*IJval+IAval);
    lambdaAstarn=beta0*(1-resAn)*(h*IJval+IAval);
    Am=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstar)*(bJ+g+lambdaJstarm);
    An=(bAm+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bAm+lambdaAstarn);
    Bn=(bAm+alpha)*(bJ+g+alpha)*(bAm+lambdaAstarn)*(bJ+g+lambdaJstar);
    % Define the invasion fitness:
    wJ=((g*am*(1-(SJval+SAval+IJval+IAval))*Am)/Bm) -1;
    wA=((g*a*(1-(SJval+SAval+IJval+IAval))*An)/Bn) -1;
    % Determine the fitness gradient:
    fitgradJworking=diff(wJ,resJm);
    fitgradAworking=diff(wA,resAn);
    fitgradJfinal=subs(fitgradJworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    fitgradAfinal=subs(fitgradAworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine the second derivatives with respect to the mutant traits
    % only:
    deriv2Jworking=diff(fitgradJworking,resJm);
    deriv2J=subs(deriv2Jworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    deriv2Aworking=diff(fitgradAworking,resAn);
    deriv2A=subs(deriv2Aworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine the cross-derivatives at (resJval, resAval) and at points
    % nearby:
    syms SJupJ SAupJ IJupJ IAupJ SJdownJ SAdownJ IJdownJ IAdownJ eps SJupA SAupA IJupA IAupA SJdownA SAdownA IJdownA IAdownA
    cderivJ1working=(subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivJ2working=(subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivJ1=subs(cderivJ1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivJ2=subs(cderivJ2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA1working=(subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivA2working=(subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivA1=subs(cderivA1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA2=subs(cderivA2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);

elseif version==4

    % Define variables and trade-offs:
    syms SJval SAval IJval IAval resJ resA resJm resAn g alpha g0 h a0 f c1a c1g c2a c2g beta0 bA bAm resJval resAval
    g=g0*(1-(c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g))));
    gm=subs(g,resJ,resJm);
    bA=1+c1a*(1-exp(-c2a*resA))/(1-exp(-c2a));
    bAm=subs(bA,resA,resAn);
    a=a0;
    bJ=1;
    
    % Set up parameters and functions:
    lambdaJstar=beta0*(1-resJ)*(h*IJval+IAval);
    lambdaJstarm=beta0*(1-resJm)*(h*IJval+IAval);
    lambdaAstar=beta0*(1-resA)*(h*IJval+IAval);
    lambdaAstarn=beta0*(1-resAn)*(h*IJval+IAval);
    Am=(bA+alpha)*(bJ+gm+alpha)+f*(bJ+gm+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+gm+alpha)*(bA+lambdaAstar)*(bJ+gm+lambdaJstarm);
    An=(bAm+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bAm+lambdaAstarn);
    Bn=(bAm+alpha)*(bJ+g+alpha)*(bAm+lambdaAstarn)*(bJ+g+lambdaJstar);
    % Define the invasion fitness:
    wJ=((gm*a*(1-(SJval+SAval+IJval+IAval))*Am)/Bm) -1;
    wA=((g*a*(1-(SJval+SAval+IJval+IAval))*An)/Bn) -1;
    % Determine the fitness gradient:
    fitgradJworking=diff(wJ,resJm);
    fitgradAworking=diff(wA,resAn);
    fitgradJfinal=subs(fitgradJworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    fitgradAfinal=subs(fitgradAworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine the second derivatives with respect to the mutant trait
    % only:
    deriv2Jworking=diff(fitgradJworking,resJm);
    deriv2J=subs(deriv2Jworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    deriv2Aworking=diff(fitgradAworking,resAn);
    deriv2A=subs(deriv2Aworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine cross-derivatives at (resJval, resAval) and at points
    % nearby:
    syms SJupJ SAupJ IJupJ IAupJ SJdownJ SAdownJ IJdownJ IAdownJ eps SJupA SAupA IJupA IAupA SJdownA SAdownA IJdownA IAdownA
    cderivJ1working=(subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivJ2working=(subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivJ1=subs(cderivJ1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivJ2=subs(cderivJ2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA1working=(subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivA2working=(subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivA1=subs(cderivA1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA2=subs(cderivA2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);

elseif version==5

    % Define variables and trade-offs:
    syms SJval SAval IJval IAval resJ resA resJm resAn alpha a g0 h a0 f c1a c1g c2a c2g beta0 bJ bJm resJval resAval
    a=a0*(1-(c1a*(1-exp(-c2a*resA))/(1-exp(-c2a))));
    an=subs(a,resA,resAn);
    bJ=1+c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g));
    bJm=subs(bJ,resJ,resJm);
    bA=1;
    g=g0;
    
    % Set up parameters and functions:
    lambdaJstar=beta0*(1-resJ)*(h*IJval+IAval);
    lambdaJstarm=beta0*(1-resJm)*(h*IJval+IAval);
    lambdaAstar=beta0*(1-resA)*(h*IJval+IAval);
    lambdaAstarn=beta0*(1-resAn)*(h*IJval+IAval);
    Am=(bA+alpha)*(bJm+g+alpha)+f*(bJm+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJm+g+alpha)*(bA+lambdaAstar)*(bJm+g+lambdaJstarm);
    An=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bA+lambdaAstarn);
    Bn=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstarn)*(bJ+g+lambdaJstar);
    % Define the invasion fitness:
    wJ=((g*a*(1-(SJval+SAval+IJval+IAval))*Am)/Bm) -1;
    wA=((g*an*(1-(SJval+SAval+IJval+IAval))*An)/Bn) -1;
    % Determine the fitness gradient:
    fitgradJworking=diff(wJ,resJm);
    fitgradAworking=diff(wA,resAn);
    fitgradJfinal=subs(fitgradJworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    fitgradAfinal=subs(fitgradAworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine second derivatives with respect to the mutant trait only:
    deriv2Jworking=diff(fitgradJworking,resJm);
    deriv2J=subs(deriv2Jworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    deriv2Aworking=diff(fitgradAworking,resAn);
    deriv2A=subs(deriv2Aworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine cross-derivatives at (resJval, resAval) and at points
    % nearby:
    syms SJupJ SAupJ IJupJ IAupJ SJdownJ SAdownJ IJdownJ IAdownJ eps SJupA SAupA IJupA IAupA SJdownA SAdownA IJdownA IAdownA
    cderivJ1working=(subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivJ2working=(subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivJ1=subs(cderivJ1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivJ2=subs(cderivJ2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA1working=(subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivA2working=(subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivA1=subs(cderivA1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA2=subs(cderivA2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);

elseif version==6

    % Define variables and trade-offs:
    syms SJval SAval IJval IAval resJ resA resJm resAn a g0 alpha h a0 f c1a c1g c2a c2g beta0 resJval resAval
    a=a0*(1-(c1a*(1-exp(-c2a*resA))/(1-exp(-c2a))))*(1-(c1g*(1-exp(-c2g*resJ))/(1-exp(-c2g))));
    an=subs(a,resA,resAn);
    am=subs(a,resJ,resJm);
    bJ=1;
    bA=1;
    g=g0;
    
    % Set up parameters and functions:
    lambdaJstar=beta0*(1-resJ)*(h*IJval+IAval);
    lambdaJstarm=beta0*(1-resJm)*(h*IJval+IAval);
    lambdaAstar=beta0*(1-resA)*(h*IJval+IAval);
    lambdaAstarn=beta0*(1-resAn)*(h*IJval+IAval);
    Am=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstar+f*lambdaJstarm*(bA+lambdaAstar);
    Bm=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstar)*(bJ+g+lambdaJstarm);
    An=(bA+alpha)*(bJ+g+alpha)+f*(bJ+g+alpha)*lambdaAstarn+f*lambdaJstar*(bA+lambdaAstarn);
    Bn=(bA+alpha)*(bJ+g+alpha)*(bA+lambdaAstarn)*(bJ+g+lambdaJstar);
    % Define the invasion fitness:
    wJ=((g*am*(1-(SJval+SAval+IJval+IAval))*Am)/Bm) -1;
    wA=((g*an*(1-(SJval+SAval+IJval+IAval))*An)/Bn) -1;
    % Determine the fitness gradient:
    fitgradJworking=diff(wJ,resJm);
    fitgradAworking=diff(wA,resAn);
    fitgradJfinal=subs(fitgradJworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    fitgradAfinal=subs(fitgradAworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine second derivatives with respect to the mutant trait only:
    deriv2Jworking=diff(fitgradJworking,resJm);
    deriv2J=subs(deriv2Jworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    deriv2Aworking=diff(fitgradAworking,resAn);
    deriv2A=subs(deriv2Aworking,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    % Determine cross-derivatives at (resJval, resAval) and at points
    % nearby:
    syms SJupJ SAupJ IJupJ IAupJ SJdownJ SAdownJ IJdownJ IAdownJ eps SJupA SAupA IJupA IAupA SJdownA SAdownA IJdownA IAdownA
    cderivJ1working=(subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradJworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivJ2working=(subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradJworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivJ1=subs(cderivJ1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivJ2=subs(cderivJ2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA1working=(subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA+eps,SJupA,SAupA,IJupA,IAupA])-subs(fitgradAworking,[resA,SJval,SAval,IJval,IAval],[resA-eps,SJdownA,SAdownA,IJdownA,IAdownA]))/(2*eps);
    cderivA2working=(subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ+eps,SJupJ,SAupJ,IJupJ,IAupJ])-subs(fitgradAworking,[resJ,SJval,SAval,IJval,IAval],[resJ-eps,SJdownJ,SAdownJ,IJdownJ,IAdownJ]))/(2*eps);
    cderivA1=subs(cderivA1working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);
    cderivA2=subs(cderivA2working,[resJm,resJ,resAn,resA],[resJval,resJval,resAval,resAval]);

end

% Display the output expressions:
disp("fitgradJ=")
disp(fitgradJfinal)
disp("fitgradA=")
disp(fitgradAfinal)
disp("deriv2J=")
disp(deriv2J)
disp("deriv2A=")
disp(deriv2A)
disp("cderivJ1=")
disp(cderivJ1)
disp("cderivJ2=")
disp(cderivJ2)
disp("cderivA1=")
disp(cderivA1)
disp("cderivA2=")
disp(cderivA2)
