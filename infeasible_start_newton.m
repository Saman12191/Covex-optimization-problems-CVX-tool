function [k,F,xstar,fstar,gfstar]=infeasible_start_newton(x)
global A  b alpha beta eta p n
X(:,1)=x;%starting point
f=feval(@fun,x);
F(1)=f;
gf=feval(@gradf,x);
H=feval(@hesf,x);
nu=zeros(p,1);
Q=[H  A';A zeros(p,p)] ;
r = [gf+A'*nu; A*x-b];
lambda=norm(r);% stopping criteria: norm(r)<eta
k=1;
while lambda>eta
    e=-Q\r;
    Dx = e(1:n) ;
    Dnu = e(n+1:n+p);
    t=1;
    while checkfeasible(x+t*Dx)==0
        t=beta*t;
    end
    while norm([feval(@gradf,x+t*Dx)+A'*(nu+t*Dnu); A*(x+t*Dx)-b]) > (1-alpha*t)*norm(r)
    t=beta*t;
    end
       x = x + t*Dx; nu = nu + t*Dnu;
       lambda=norm(r);
       gf=feval(@gradf,x);
       k=k+1;
       X(:,k)=x;%starting point
       f=feval(@fun,x);
       F(k)=f;
       gf=feval(@gradf,x);
       H=feval(@hesf,x);
       Q=[H  A';A zeros(p,p)] ;
       r = [gf+A'*nu; A*x-b];
       lambda=norm(r);% stopping criteria: norm(r)<eta
end

T(k)=0;
K=1:k;
%figure(3)
% plot(K,F, ' m ','LineWidth',1,'MarkerSize',2);%arghavani
% title({'Plot of objective function  versus iteration number: newton method'; ['n=',num2str(n)];['\alpha=',num2str(alpha)];['\beta=',num2str(beta)]})
% axis normal
% xlabel('k=iteration number');
% ylabel('f(x_k)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(5)
% plot(K,T, ' b ','LineWidth',1,'MarkerSize',2)
% title({'Plot of  step length versus iteration number: newton method'; ['n=',num2str(n)];['\alpha=',num2str(alpha)];['\beta=',num2str(beta)]})
% axis normal
% xlabel('k=iteration number');
% ylabel('t_k: step length)')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(4)
% plot(K,log(F-F(k)))
% title({'Plot of log(f-p^*) versus iteration number: newton method'; ['n=',num2str(n)];['\alpha=',num2str(alpha)];['\beta=',num2str(beta)]})
% axis normal
% xlabel('k=iteration number');
% ylabel('log(f(x_k)-p^*)')
k
xstar=x;
fstar=f;
gfstar=gf;