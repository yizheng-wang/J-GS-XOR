%%The first question of experiment 3.3
clear
n = 6 ; % creat hilbert matrix of order n
H = hilb(n);% creat hilbert matrix of order n
Xstar = ones(n, 1); % specify exact solution for obtaining b 
b = H * Xstar; % creat b of AX=b 
%%
xgauss = gauss(H, b);
errorgauss(p) = norm(xgauss - Xstar);
m=0;
for i = 300:1:1290 % Calculate the relationship between GS method iteration 
                    %  times and error
     m = m + 1; 
     [xgs(:, m), kgs(m)] = GS(H, b, i,10e-30);  
     yxgs(m) = norm(xgs(:, m)-ones(n, 1));
end

s = 0; % s stands for different W
for W = 0:0.1:2 % Calculate the relationship between the number of iterations
                %    of different W and the error
 %Calculate the relationship between GS method iteration times and error            
    m = 0; 
    s = s +1;
    for i = 300:1:1290 
        m = m + 1;  % m represents different number of iterations
        [xxor(:, m, s), kxor(m,s)] = XOR(H, b, W, i,10e-30);  
        yxxor(m, s) = norm(xxor(:, m, s)-ones(n, 1));
    end
end
%% plot all comparisons of W
x = [300:1:1290];
for i = 2:20
plot(x, yxxor(:, i), 'Color', rand(1,3))
hold on 
end
legend('w = 0.1','w = 0.2','w = 0.3','w = 0.4','w = 0.5','w = 0.6',...
'w = 0.7','w = 0.8','w = 0.9','w = 1.0','w = 1.1','w = 1.2','w = 1.3',...
'w = 1.4','w = 1.5','w = 1.6','w = 1.7','w = 1.8','w = 1.9')
xlabel('the number of iteration')
ylabel('the error under two norm')
%% graphing all comparisons of W=0.2 and W=1
x = [300:1000:300000 ];
m = 0; % Calculate the relationship between GS method iteration times and error
for i = 300:1000:300000 
        m = m + 1;  % m represents different number of iterations
        [xxor1(:, m), kxor1(m)] = XOR(H, b, 1, i,10e-30);  
        yxxor1(m) = norm(xxor1(:, m)-ones(n, 1));
        [xxor2(:, m), kxor2(m)] = XOR(H, b, 0.2, i,10e-30);  
        yxxor2(m) = norm(xxor2(:, m)-ones(n, 1));
end

plot(x, yxxor1, x, yxxor2)
legend('w = 1.0','w = 0.2')
xlabel('the number of iteration')
ylabel('error under two norm')


