function [resJss,resAss]=singstrats_at_0or1(fitgradJval,fitgradAval,resJvalvec,resAvalvec,SSres,maxsingstrats,R0counter)

% This function determines whether there are singular strategies when the
% adult or juvenile resistance is at zero or one (and hence where there may
% not be sign changes in both fitness gradients).

finished=0;
resJss=-10+zeros(1,maxsingstrats);
resAss=-10+zeros(1,maxsingstrats);
foundsofar=1;

for j=1:SSres-1
    % Determine fitness gradients for each value of resA when resJ is zero:
    if fitgradJval(1,j)<0 && fitgradAval(1,j)>0 && fitgradJval(1,j+1)<0 && fitgradAval(1,j+1)<0 && resJvalvec(1)==0
        newSS=(resAvalvec(j)+resAvalvec(j+1))/2;
        resJss(foundsofar)=0;
        resAss(foundsofar)=newSS;
        foundsofar=foundsofar+1;
        finished=1;
    % Determine fitness gradients for each value of resA when resJ is one:
    elseif fitgradJval(end,j)>0 && fitgradAval(end,j)>0 && fitgradJval(end,j+1)>0 && fitgradAval(end,j+1)<0 && resJvalvec(end)==1
        newSS=(resAvalvec(j)+resAvalvec(j+1))/2;
        resJss(foundsofar)=1;
        resAss(foundsofar)=newSS;
        foundsofar=foundsofar+1;
        finished=1;
    % Determine fitness gradients for each value of resJ when resA is zero:
    elseif fitgradJval(j,1)>0 && fitgradJval(j+1,1)<0 && fitgradAval(j,1)<0 && fitgradAval(j+1,1)<0 && resAvalvec(1)==0
        newSS=(resJvalvec(j)+resJvalvec(j+1))/2;
        resJss(foundsofar)=newSS;
        resAss(foundsofar)=0;
        foundsofar=foundsofar+1;
        finished=1;
    % Determine fitness gradients for each value of resA when resJ is zero:
    elseif fitgradJval(j,end)>0 && fitgradJval(j+1,end)<0 && fitgradAval(j,end)>0 && fitgradAval(j+1,end)>0 && resAvalvec(end)==1
        newSS=(resJvalvec(j)+resJvalvec(j+1))/2;
        resJss(foundsofar)=newSS;
        resAss(foundsofar)=1;
        foundsofar=foundsofar+1;
        finished=1;
    end
end

% Other cases:
if finished==0
    if R0counter==SSres^2 % The disease is never viable
        resJss=[0,-10+zeros(1,maxsingstrats-1)];
        resAss=[0,-10+zeros(1,maxsingstrats-1)];
    elseif sum(isnan(fitgradJval),'all')==SSres^2 % hosts are never viable
        resJss=-10+zeros(1,maxsingstrats);
        resAss=-10+zeros(1,maxsingstrats);
    % Singular strategies at the corners:
    elseif sum(fitgradJval<0,'all')==0 && sum(fitgradAval<0,'all')==0 && resJvalvec(end)==1 && resAvalvec(end)==1
        resJss=[1,-10+zeros(1,maxsingstrats-1)];
        resAss=[1,-10+zeros(1,maxsingstrats-1)];
    elseif sum(fitgradJval<0,'all')==0 && sum(fitgradAval>0,'all')==0 && resJvalvec(end)==1 && resAvalvec(1)==0
        resJss=[1,-10+zeros(1,maxsingstrats-1)];
        resAss=[0,-10+zeros(1,maxsingstrats-1)];
    elseif sum(fitgradJval>0,'all')==0 && sum(fitgradAval<0,'all')==0 && resJvalvec(1)==0 && resAvalvec(end)==1
        resJss=[0,-10+zeros(1,maxsingstrats-1)];
        resAss=[1,-10+zeros(1,maxsingstrats-1)];
    elseif sum(fitgradJval>0,'all')==0 && sum(fitgradAval>0,'all')==0 && resJvalvec(1)==0 && resAvalvec(1)==0
        resJss=[0,-10+zeros(1,maxsingstrats-1)];
        resAss=[0,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradJval(1,1)<0 && fitgradAval(1,1)<0 && resJvalvec(1)==0 && resAvalvec(1)==0
        resJss=[0,-10+zeros(1,maxsingstrats-1)];
        resAss=[0,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradJval(1,end)<0 && fitgradAval(1,end)>0 && resJvalvec(1)==0 && resAvalvec(end)==1
        resJss=[0,-10+zeros(1,maxsingstrats-1)];
        resAss=[1,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradJval(end,1)>0 && fitgradAval(end,1)<0 && resJvalvec(end)==1 && resAvalvec(1)==0
        resJss=[1,-10+zeros(1,maxsingstrats-1)];
        resAss=[0,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradJval(end,end)>0 && fitgradAval(end,end)>0 && resJvalvec(end)==1 && resAvalvec(end)==1
        resJss=[1,-10+zeros(1,maxsingstrats-1)];
        resAss=[1,-10+zeros(1,maxsingstrats-1)];
    else % No singular strategies
        resJss=-10+zeros(1,maxsingstrats);
        resAss=-10+zeros(1,maxsingstrats);
    end
end

end