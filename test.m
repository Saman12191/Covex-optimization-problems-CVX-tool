
clc; clear all;
global A p n alpha beta eta b
n=100;
p=30; 
alpha=.2;% between 0.005 and .5
beta=.5;% between 0 and 1
eta=10^-6;
x=rand(n,1);
A=rand(p,n);
while rank(A)~=p
    A=rand(p,n);
end
x0=x;
b=A*x0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=eye(n);
Z=I-A'*(A*A')^(-1)*A;
%%%%%%%%%% (a)
disp('a:')
figure(1);
[k1,F1,xstar1,fstar1,gfstar1]=newton_method(x0);
K=1:k1;
hold on
plot(K,log(F1-fstar1),'k','linewidth',3);
disp(['feasibility error=',num2str(norm(A*xstar1-b))]);
%disp(['The values of reduced gradient norm at  x^* is: ', num2str(norm(Z*gfstar1))]);
%%%%%%%% (b)
[k2,F2,xstar2,fstar2,gfstar2]=infeasible_start_newton(x);
K=1:k2;
plot(K,log(F2-fstar2),'g--','linewidth',3);
disp(['feasibility error=',num2str(norm(A*xstar2-b))]);
%disp(['The values of reduced gradient norm at  x^* is: ', num2str(norm(Z*gfstar2))]);
%%%%%%%%%%%%%% (c)
x=ones(n,1);
[k3,F3,xstar3,fstar3,gfstar3]=infeasible_start_newton(x);
K=1:k3;
plot(K,log(F3-fstar3),'m','linewidth',3);
disp(['feasibility error=',num2str(norm(A*xstar3-b))]);
% %%%%%%%%%%%%%%%%% cvx
% cvx_setup
% cvx_begin
% variable x(n)
% minimize(sum(x.*log(x)))
% subject to
% A*x == b;of
% cvx_end
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title({'Plot of log(f-p^*) versus iteration number: newton method'; ['n=',num2str(n)];['\alpha=',num2str(alpha)];['\beta=',num2str(beta)]})
axis normal
xlabel('k=iteration number');
ylabel('log(f(x_k)-p^*)');
legend('newton method with x0=x^','infeasible newton method with x0=x^','infeasible newton method with x0=1')
hold off
