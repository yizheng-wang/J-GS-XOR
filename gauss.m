function [ x ] = gauss( A, b, order )
%Automatic Gaussian elimination and Gaussian elimination 
%with manual selection of main element, solve Ax=b
%    input£º A - matrix to be solved
%            b - Inhomogeneous terms
%    output: order-the selection of pivots, where 0 is automatic, 
%                   1 is column pivot elimination, -1 is the smallest pivot,
%                   and the rest 2, 3, 4... are the order of the pivot
%                   The order is the second largest, the third largest,
%                   the fourth largest..., -2, -3, -4... are 
%                   in the order of the pivot size,followed by the  
%                   second smallest, third smallest, fourth smallest ...
%                   The default order parameter is 0
if nargin < 3
    order = 0;
end
[rownum, colnum] = size(A);
Ab = [A, b]; % consist augmented matrix
% the process of elimination
for j = 1 : colnum-1
            % selection of the main element
    if order ~= 0
            needorder = Ab(j:colnum, j);
            [~, I] = sort(needorder, 'descend');
    end
        
    if order >0
            getindex = I(order);
            Ab([getindex+j-1, j],:) = Ab([j, getindex+j-1],:);
    end
        
    if order <0
            getindex = I(end+order+1);
            Ab([getindex+j-1, j],:) = Ab([j, getindex+j-1],:);            
    end
    
    for i = j+1 : rownum
        l = Ab(i, j)/Ab(j, j);
        Ab(i, :) = Ab(i, :)- l * Ab(j, :);
    end
end
% back propogation process
x = zeros(colnum, 1);
for i = colnum : -1 : 1
    sum = 0 ;
    for m = i + 1 :colnum
        sum = sum + Ab(i, m) * x(m);
    end
    x(i) = (Ab(i, colnum + 1) - sum)/Ab(i, i);
end

end

