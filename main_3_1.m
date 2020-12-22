%% n = 10
A = diag(ones(1, 10)*6)+diag(ones(1, 9), 1)+diag(ones(1, 9)*8, -1);
cond(A)
b = ones(10,1)*15;
b(1) = 7;
b(10) = 14;
[x] = gauss(A,b);
%% n = 20
A = diag(ones(1, 20)*6)+diag(ones(1, 19), 1)+diag(ones(1, 19)*8, -1);
cond(A)
b = ones(20,1)*15;
b(1) = 7;
b(20) = 14;
[x] = gauss(A,b);
%% 3.1.4
for n = 1:40

H = hilb(n);
xex = ones(n, 1);
b = H * xex;
xauto = gauss(H,b);
xmax = gauss(H,b,1);
error_xauto(n) = norm(xex-xauto);
error_xmax(n) = norm(xex-xmax);
end
plot(1:40, error_xauto,'r')
hold on 
plot(1:40, error_xmax)
legend('auto','max')
xlabel('Hilbert order')
ylabel('Absolute error under two norm')
title(['The relationship between the order of Hilbert matrix and the absolute'...
    'error of column principal element method and Gaussian elimination method'])