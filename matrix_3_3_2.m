%% The second question of experiment 3.3
clear
p = 0;
for n = 10 : 20
% creat hilbert matrix of order n
p = p + 1;
H = hilb(n);% creat hilbert matrix of order n
Xstar = ones(n, 1); % specify exact solution for obtaining b 
b = H * Xstar; % creat b of AX=b 
xgauss = gauss(H, b);
errorgauss(p) = norm(xgauss - Xstar);

[xgs, kgs(p)] = GS(H, b, 10000, errorgauss(p));  
yxgs(p) = norm(xgs-Xstar);

s = 0; % s stands for different W

for W = 0.1:0.1:1.9 % Calculate the relationship between the number of 
%                       iterations of different W and the error
    s = s +1;
    [xxor, kxor(s, p)] = SOR(H, b, W, 30000, 0.1);  
    yxxor(s, p) = norm(xxor - Xstar);
end
end
%% Plot all comparisons of W
x = [10:20];
for i = 1:19
plot(x, kxor(i, :), 'Color', rand(1,3))
hold on 
end
legend('w = 0.1','w = 0.2','w = 0.3','w = 0.4','w = 0.5','w = 0.6',...
'w = 0.7','w = 0.8','w = 0.9','w = 1.0','w = 1.1','w = 1.2','w = 1.3',...
'w = 1.4','w = 1.5','w = 1.6','w = 1.7','w = 1.8','w = 1.9')
xlabel('the order of matrix')
ylabel('the number of iteration needed')
%% Plot all comparisons of W=0.2 and W=1
x = [300:1000:300000 ];
m = 0; % Calculate the relationship between GS method iteration times and error
for i = 300:1000:300000 
        m = m + 1;  % m represents different number of iterations
        [xxor1(:, m), kxor1(m)] = SOR(H, b, 1, i,10e-30);  
        yxxor1(m) = norm(xxor1(:, m)-ones(n, 1));
        [xxor2(:, m), kxor2(m)] = SOR(H, b, 0.2, i,10e-30);  
        yxxor2(m) = norm(xxor2(:, m)-ones(n, 1));
end

plot(x, yxxor1, x, yxxor2)
legend('w = 1.0','w = 0.2')
xlabel('the number of iteration')
ylabel('the error under two norm')


