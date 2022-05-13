function [SJval,SAval,IJval,IAval,problem_marker]=endemic_equilibrium_function(resJval,resAval,aval,gval,bJval,bAval,beta0,h,f,alpha,initvec,orig_tmax)

% This function determines the endemic equilibrium of the ecological
% system.

tol = 1e-3;
betaJ=beta0*(1-resJval);
betaA=beta0*(1-resAval);
problem_marker=0;

% 'ode_fast' is an ODE solver.
[t,y] = ode_fast(orig_tmax,aval,h,f,gval,bJval,bAval,betaJ,betaA,alpha,initvec);

% Check that the system has actually reached equilibrium:
t1=find(t>t(end)*0.9,1);

% If the system has not reached equilibrium, then we need to extend the
% time:
if(any(range(y(t1:end,:))>tol))
    [t,y] = ode_fast(orig_tmax*9,aval,h,f,gval,bJval,bAval,betaJ,betaA,alpha,y(end,:));
    t=t+orig_tmax;
    
    % Check again and extend again...
    t2=find(t>t(end)*0.9,1);
    if(any(range(y(t2:end,:))>tol)) 
        [t,y] = ode_fast(orig_tmax*90,aval,h,f,gval,bJval,bAval,betaJ,betaA,alpha,y(end,:));
        t=t+orig_tmax*10;
        t3=find(t>t(end)*0.9,1);
        
        % Check again and extend again...
        if(any(range(y(t3:end,:))>tol))
            [t,y] = ode_fast(orig_tmax*900,aval,h,f,gval,bJval,bAval,betaJ,betaA,alpha,y(end,:));
            t=t+orig_tmax*100;
            t4=find(t>t(end)*0.9,1);
            
            % Display an error if the system has still not reached
            % equilibrium:
            if(any(range(y(t4:end,:))>50*tol))
                problem_marker=1;
            end
        end
    end
end

% The output of this function is the equilibrium state of the system:
SJval=y(end,1);
SAval=y(end,2);
IJval=y(end,3);
IAval=y(end,4);

end