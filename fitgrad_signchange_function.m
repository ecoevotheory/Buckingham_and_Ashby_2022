function [markerJ,markerA,singstratpointerJ,singstratpointerA]=fitgrad_signchange_function(fitgradJval,fitgradAval)

% This function finds where the fitness gradients change sign.

n = length(fitgradJval);
markerJ = zeros(n,n,2);
markerA = zeros(n,n,2);
for k=1:4
    
    if(k<3)
        input = sign(fitgradJval);
    else
        input = sign(fitgradAval);
    end
    
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
    
    INPUT(2:(n+1),2:(n+1)) = input;
    UP = [INPUT(2:(n+2),:);zeros(1,n+2)];
    DOWN = [zeros(1,n+2);INPUT(1:(n+1),:)];
    LEFT = [INPUT(:,2:(n+2)),zeros(n+2,1)];
    RIGHT = [zeros(n+2,1),INPUT(:,1:(n+1))];
    
    A = INPUT==UP;
    B = INPUT==DOWN;
    C = INPUT==LEFT;
    D = INPUT==RIGHT;
    
    OUTPUT = (A&B&C&D);
    output = 1-OUTPUT(2:(end-1),2:(end-1));
    
    if(k<3)
        markerJ(:,:,k) = output;
    else
        markerA(:,:,k-2) = output;
    end
end

markerJ = double(markerJ(:,:,1)&markerJ(:,:,2));
markerA = double(markerA(:,:,1)&markerA(:,:,2));

markerJ(isnan(fitgradJval)|isnan(fitgradAval)) = 0;
markerA(isnan(fitgradJval)|isnan(fitgradAval)) = 0;

list = find(markerA&markerJ);
list(isnan(fitgradJval(list)) | isnan(fitgradAval(list)))=[];

[singstratpointerJ,singstratpointerA] = ind2sub(size(markerA),list);
singstratpointerJ=unique(singstratpointerJ);
singstratpointerA=unique(singstratpointerA);

end
