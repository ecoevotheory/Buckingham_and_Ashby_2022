function [fitgradJval,fitgradAval,resJvalvec,resAvalvec,R0counter]=fitness_gradients(SSres,startJ,startA,finJ,finA,a0,g0,c1a,c2a,c1g,c2g,beta0,h,alpha,f,initvec,version,orig_tmax)

% This function determines the two fitness gradients at different values of
% juvenile and adult resistance.

% Set up vectors to use later:
fitgradJval=zeros(SSres);
fitgradAval=zeros(SSres);
resJvalvec = linspace(startJ,finJ,SSres);
resAvalvec = linspace(startA,finA,SSres);
R0counter=0;

% For each value of adult and juvenile resistance we determine the fitness
% gradients:
for j=1:SSres
    for k=1:SSres
        resJval=resJvalvec(j);
        resAval=resAvalvec(k);
        
        % The trade-offs will be different depending on which version
        % of the model we are using:
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
        
        % Check that the host population is viable in the absence
        % of disease
        hostviability=aval*gval>bAval*(bJval+gval);
        
        if(hostviability)
            
            % Check that disease is viable
            R0=beta0*(aval*gval-bAval*(bJval+gval))*((bAval*(h*bAval+gval+h*alpha)*(1-resJval))+(gval*(bJval+gval+alpha)*(1-resAval)))/(aval*gval*(bAval+alpha)*(bJval+gval)*(bJval+gval+alpha));
            
            if(R0>1)
                % Find the endemic equilibrium:
                [SJval,SAval,IJval,IAval,problem_marker]=endemic_equilibrium_function(resJval,resAval,aval,gval,bJval,bAval,beta0,h,f,alpha,initvec,orig_tmax);
                
                % These are the expressions for the fitness gradients
                % generated by the code called 'fitgrad_functions'.
                if problem_marker==0
                    if version==1
                        fitgradJval(j,k)=(a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) + (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (a0*c1g*c2g*g0^2*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2);
                        fitgradAval(j,k)=(a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0^2*f*(IAval + IJval*h)^2*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1));
                    elseif version==2
                        fitgradJval(j,k)=(a0*g0*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) + (c1g*c2g*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) - (beta0*c1g*c2g*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (a0*g0*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2);
                        fitgradAval(j,k)=(a0*g0*(beta0*f*(IAval + IJval*h)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (a0*g0*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (a0*c1a*c2a*g0*exp(-c2a*resAval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1));
                    elseif version==3
                        fitgradJval(j,k)=(a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (a0*beta0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^2*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (a0*beta0*f*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1));
                        fitgradAval(j,k)=- (a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha + g0 + 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha + g0 + 1))/(exp(-c2a) - 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)) - (a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + g0 + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2);
                    elseif version==4
                        fitgradJval(j,k)=(a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (c1g*c2g*g0*exp(-c2g*resJval)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*g0*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) - (c1g*c2g*g0*exp(-c2g*resJval))/(exp(-c2g) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)^2*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*c1g*c2g*g0*exp(-c2g*resJval)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) - (a0*c1g*c2g*g0^2*exp(-c2g*resJval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)^2);
                        fitgradAval(j,k)=(a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + (c1a*c2a*exp(-c2a*resAval)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))/(exp(-c2a) - 1) - beta0*f*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(IAval + IJval*h)*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (a0*g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1a*c2a*exp(-c2a*resAval))/(exp(-c2a) - 1))*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/(((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)^2*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1)) + (a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(resJval - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1) - (alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1) + beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - beta0*(IAval + IJval*h)*(resAval - 1) + 1)*(alpha + (c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) + 1)^2*(g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + beta0*(IAval + IJval*h)*(resJval - 1) - 1)*(alpha - g0*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1) + 1));
                    elseif version==5
                        fitgradJval(j,k)=- (a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h) - (c1g*c2g*exp(-c2g*resJval)*(alpha + 1))/(exp(-c2g) - 1) + (beta0*c1g*c2g*f*exp(-c2g*resJval)*(IAval + IJval*h)*(resAval - 1))/(exp(-c2g) - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0*(IAval + IJval*h) + (c1g*c2g*exp(-c2g*resJval))/(exp(-c2g) - 1))*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)^2);
                        fitgradAval(j,k)=(a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) - (a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*(beta0^2*f*(IAval + IJval*h)^2*(resJval - 1) - beta0*f*(IAval + IJval*h)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)) + (a0*c1a*c2a*g0*exp(-c2a*resAval)*((alpha + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1)*(alpha + g0 + (c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) + 1));
                    elseif version==6
                        fitgradJval(j,k)=(a0*beta0*f*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*(IAval + IJval + SAval + SJval - 1))/((alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) + (a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)^2*(alpha + g0 + 1)) - (a0*c1g*c2g*g0*exp(-c2g*resJval)*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2g) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1));
                        fitgradAval(j,k)=- (a0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(beta0*f*(IAval + IJval*h)*(alpha + g0 + 1) - beta0^2*f*(IAval + IJval*h)^2*(resJval - 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) - (a0*beta0*g0*((c1a*(exp(-c2a*resAval) - 1))/(exp(-c2a) - 1) - 1)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*(IAval + IJval*h)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((beta0*(IAval + IJval*h)*(resAval - 1) - 1)^2*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1)) - (a0*c1a*c2a*g0*exp(-c2a*resAval)*((c1g*(exp(-c2g*resJval) - 1))/(exp(-c2g) - 1) - 1)*((alpha + 1)*(alpha + g0 + 1) + beta0*f*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(IAval + IJval*h)*(resJval - 1) - beta0*f*(IAval + IJval*h)*(resAval - 1)*(alpha + g0 + 1))*(IAval + IJval + SAval + SJval - 1))/((exp(-c2a) - 1)*(beta0*(IAval + IJval*h)*(resAval - 1) - 1)*(alpha + 1)*(g0 - beta0*(IAval + IJval*h)*(resJval - 1) + 1)*(alpha + g0 + 1));
                    end
                    
                else
                    fitgradJval(j,k)=NaN;
                    fitgradAval(j,k)=NaN;
                end
                
            else % If there is no disease then resistance should always decrease
                fitgradJval(j,k)=-1;
                fitgradAval(j,k)=-1;
                R0counter=R0counter+1;
            end
        else % If the host population is not viable then there are no hosts to exhibit the trait and hence no fitness gradients
            fitgradJval(j,k)=NaN;
            fitgradAval(j,k)=NaN;
        end
    end
end