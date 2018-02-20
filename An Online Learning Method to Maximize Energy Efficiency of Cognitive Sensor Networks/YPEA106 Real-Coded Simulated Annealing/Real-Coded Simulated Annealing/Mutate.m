%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA106
% Project Title: Real-Coded Simulated Annealing in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function y=Mutate(x,mu,sigma,VarMin,VarMax)

    A=(rand(size(x))<=mu);
    J=find(A==1);

    y=x;
    y(J)=x(J)+sigma*randn(size(J));

    % Clipping
    y=max(y,VarMin);
    y=min(y,VarMax);
    
end