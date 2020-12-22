function [x, k] = J(a, b, TOL, N, XO)
%J Iterative method Input A, B, and output X0 through Jacobian iteration method
%   input£ºA: Coefficient matrix
%          B: Inhomogeneous terms
%          tol: Threshold, default 10^-5
%          n: The default maximum number of iterations, the default is 10^5
%          Xi:Default 0 vector
%   output£ºX0:output
if nargin == 2
    TOL = 10e-5;
    N = 10^5;
    XO = zeros(length(a), 1);
end

if nargin == 3
    N = 10^5;
    XO = zeros(length(a), 1);
end

if nargin == 4
    XO = zeros(length(A), 1);
end

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
n = length(a);
while k <= N % While k is less than or equal to our max iterations N:
    for i = 1:n % For all entries 1 to n: 
        for j = 1:n % For all entries 1 to n:
   % If j is not equal to i, since we exclude this term in our approximations:
            if j ~= i 
                s = s + a(i,j).*XO(j); % compute s
            end % end if loop
        end % end for loop
        % calculate the approximate ith iteration solution for vector x
        x(i) = (1./a(i,i)).*(-s + b(i)); 
        s = 0; % reset s back to zero
    end % end for loop
    % If our change between iterations is less than our tolerance:
    if norm( x' - XO, 2 ) < TOL 
        x; % print x
        break; % break the while loop
    end % end if loop
    k = k+1; % Otherwise, our tolerance has not been reached, 
             % so increment iterative count by one to go through the loop again
    for i = 1:n % For all entries 1 to n:
        XO(i) = x(i)'; % set our x vector equal to XO for the next iteration
    end % end for loop
end % end while loop
end

