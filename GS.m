function [x, k] = GS(a, b, N, TOL, XO)
%GS Iterative method 
%   input£ºA: Coefficient matrix
%          B: Inhomogeneous terms
%          tol: Threshold, default 10^-5
%          n: The default maximum number of iterations, the default is 10^5
%          Xi:Default 0 vector
%   output£ºX0:output
Xstar = ones(length(a), 1);
if nargin == 2
    TOL = 10e-5;
    N = 10^5;
    XO = zeros(length(a), 1);
end

if nargin == 3
    TOL = 10e-5;
    XO = zeros(length(a), 1);
end

if nargin == 4
    XO = zeros(length(a), 1);
end

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
n = length(a);
p = 0; % let p be 0. We will use this variable for summation functions
while k <= N % While k is less than or equal to our max iterations N
    for i = 1:n % For entries 1 to n:
        for j = 1:i-1 % For entries 1 to i-1:
            s = s + a(i,j).*x(j); % compute the first summation 
        end % end for loop
        for j = i+1:n % For entries i+1 to n:
            p = p + a(i,j).*XO(j); % compute second summation
        end % end for loop
        % calculate the values in our approximated x vector
        x(i) = (1./a(i,i)).*(-s - p + b(i)); 
        s = 0; % reset s back to zero
        p = 0; % reset p back to zero
    end % end for loop
    % If our change between iterations is less than our tolerance:
    if norm(x' - Xstar) < TOL 
        x; % print the approximation
        break; % break while loop
    end % end if loop
    % Otherwise, our tolerance has not been reached, so increment 
    % iterative count by one to go through the loop again
    k = k+1; 
    for i = 1:n % For entries 1 to n:
        XO(i) = x(i)'; % set our x vector equal to XO for the next iteration
    end % end for loop
    
end % end while loop
x = x';


end

