%%RUNDOWN
%%Problem 2
load west0479;
A = west0479;

spy(A);

lam = eig(full(A));
plot(real(lam),imag(lam),'r+');
title('Real Eigenvalue of west0479')
print(gcf, '-dpng', '-r280', 'p2.png');
%%Problem 3

v0_3 = rand(size(A, 2), 1); 
j  = 30; 

[V, H, iter] = myarnoldi(A, v0_3, j, 1);

res1 = norm(A * V(:,1:j) - V * H ) 
res2 = norm(eye(j+1) - V' * V)


%%Problem 4

Hj = H(1:j, 1:j); 
ritz = eig(Hj); 
plot(real(ritz), imag(ritz), 'bo'); 
hold on; 
plot(real(lam),imag(lam),'r+');
hold off; 
lgd = legend('ritz-30', 'exact');
set(lgd,'FontSize', 16); 
print(gcf, '-dpng', '-r280', 'p4.png');

%%Problem 5
j = 10; 
[V1, H1, iter] = myarnoldi(A, v0_3, j, 1);
H1j = H1(1:j, 1:j); 
ritz1 = eig(H1j); 
res1 = norm(A * V1(:,1:j) - V1 * H1 ) 
res2 = norm(eye(j+1) - V1' * V1)
j = 20;
[V2, H2, iter] = myarnoldi(A, v0_3, j, 1);
H2j = H2(1:j, 1:j); 
ritz2 = eig(H2j); 
res1 = norm(A * V2(:,1:j) - V2 * H2 ) 
res2 = norm(eye(j+1) - V2' * V2)

plot(real(ritz1), imag(ritz1), 'y*'); 
hold on;
plot(real(ritz2), imag(ritz2), 'g*'); 
hold on;
plot(real(ritz), imag(ritz), 'bo'); 
hold on;
plot(real(lam),imag(lam),'r+');
hold off; 
lgd = legend('ritz-10','ritz-20','ritz-30', 'exact');
set(lgd,'FontSize', 16); 
print(gcf, '-dpng', '-r280', 'p5.png');

%%Problem 6
j  = 30; 
res = [];
res_r = [];

for j = 0:10:100;
    [Vi, Hi, iter] = arnoldi(A, v0_3, j, 0);
    res =[res norm(eye(j+1) - Vi' * Vi)];
    
    [Vr, Hr, iter] = arnoldi(A, v0_3, j, 1);
    res_r =[res_r norm(eye(j+1) - Vr' * Vr)];
end
    
plot(0:10:100, log10(res), 'bx-');
hold on;
plot(0:10:100, log10(res_r), 'ro-');
hold off; 
legend('without','reorthogonization','Location','northwest');
title('Implementations of Arnoldi method','FontSize', 16);
xlabel('iteration step');
ylabel('loss of orthogonality');
print(gcf, '-dpng', '-r280', 'p62.png');


    
j = 30;
[V3, H3, iter] = arnoldi(A, v0_3, j, 0);
%res1 = norm(A * V(:,1:j) - V * H ); 
res2_3 = norm(eye(j+1) - V3' * V3)
H3j = H3(1:j, 1:j); 
ritz3 = eig(H3j); 

j = 60;
[V6, H6, iter] = arnoldi(A, v0_3, j, 0);
%res1 = norm(A * V(:,1:j) - V * H ); 
res2_6 = norm(eye(j+1) - V6' * V6)
H6j = H6(1:j, 1:j); 
ritz6 = eig(H6j); 

plot(real(ritz3), imag(ritz3), 'bo'); 
hold on;
plot(real(ritz6), imag(ritz6), 'g*'); 
hold on;
plot(real(lam),imag(lam),'r+');
hold off; 
legend('ritz-30','ritz-60', 'exact');
print(gcf, '-dpng', '-r280', 'p6.png');
%%Problem 7

tau = 10; 
j = 30; 
[V71, H71, iter] = myarnolditau(A, v0_3, j, 0, tau);
H71j = H71(1:j, 1:j); 
ritz71 = eig(H71j).^(-1) + tau;



tau = -10; 
[V72, H72, iter] = myarnolditau(A, v0_3, j, 0, tau);
H72j = H72(1:j, 1:j); 
ritz72 = eig(H72j).^(-1) + tau; 

plot(real(lam),imag(lam),'r+');
hold on;
plot(real(ritz71), imag(ritz71), 'bo'); 
hold on;
plot(real(ritz72), imag(ritz72), 'g*'); 
hold off; 

legend('exact','tau = 10','tau = -10');
print(gcf, '-dpng', '-r280', 'p7.png');