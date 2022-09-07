function [markerJ,markerA,numberofsingstrats]=fitgrad_signchange_function(fitgradJval,fitgradAval)

% This function finds where the fitness gradients change sign.

n = length(fitgradJval);
markerJ = zeros(n,n,2);
markerA = zeros(n,n,2);
for k=1:4
    
    % First, we consider sign changes in fitgradJval and then we repeat for
    % fitgradAval:
    if(k<3)
        input = sign(fitgradJval);
    else
        input = sign(fitgradAval);
    end
    
    % We surround the matrix of fitness gradient signs with either 1's or
    % -1's
    INPUT = zeros(n+2);
    if(rem(k,2)==0)
        input(isnan(input))=-1;
        INPUT(1,:)=-1;
        INPUT(n+2,:)=-1;
        INPUT(:,1)=-1;
        INPUT(:,n+2)=-1;
    else
        input(isnan(input))=1;
        INPUT(1,:)=1;
        INPUT(n+2,:)=1;
        INPUT(:,1)=1;
        INPUT(:,n+2)=1;
    end
    
    % We determine whether the signs are the same above, below, to the left
    % of and to the right of each entry:
    INPUT(2:(n+1),2:(n+1)) = input;
    UP = [INPUT(2:(n+2),:);zeros(1,n+2)];
    DOWN = [zeros(1,n+2);INPUT(1:(n+1),:)];
    LEFT = [INPUT(:,2:(n+2)),zeros(n+2,1)];
    RIGHT = [zeros(n+2,1),INPUT(:,1:(n+1))];
    A = INPUT==UP;
    B = INPUT==DOWN;
    C = INPUT==LEFT;
    D = INPUT==RIGHT;
    
    % If they are not all the same then we have detected a sign change:
    OUTPUT = (A&B&C&D);
    output = 1-OUTPUT(2:(end-1),2:(end-1));
    
    % The places where these sign changes occur give our output:
    if(k<3)
        markerJ(:,:,k) = output;
    else
        markerA(:,:,k-2) = output;
    end
end

% This method also wrongly detects sign changes at the edge of the fitness
% gradient matrix. For this reason, we look for sign changes both when the
% edge is assumed to be positive and when it is assumed to be negative, and
% require that both give a sign change:
markerJ = double(markerJ(:,:,1)&markerJ(:,:,2));
markerA = double(markerA(:,:,1)&markerA(:,:,2));

% We do not want changes from positive or negative to NaN (when the host is
% not viable) to count as sign changes, so we remove them:
markerJ(isnan(fitgradJval)|isnan(fitgradAval)) = 0;
markerA(isnan(fitgradJval)|isnan(fitgradAval)) = 0;

% We determine where both fitgradJval and fitgradAval change sign:
list = find(markerA&markerJ);
list(isnan(fitgradJval(list)) | isnan(fitgradAval(list)))=[];
[singstratpointerJ,~] = ind2sub(size(markerA),list);
% This gives us a list of locations of singular strategies:
singstratpointerJ=unique(singstratpointerJ);

% We then know how many potential singular strategies have been found:
numberofsingstrats=length(singstratpointerJ);

end
