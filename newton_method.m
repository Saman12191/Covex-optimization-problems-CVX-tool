function [k,F,xstar,fstar,gfstar]=newton_method(x)
global A  alpha beta eta p n
X(:,1)=x;%starting point
f=feval(@fun,x);
F(1)=f;
gf=feval(@gradf,x);
H=feval(@hesf,x);
Q=[H  A';A zeros(p,p)] ;
v=[gf; zeros(p,1)];
e=-Q\v;
d=e(1:n);%L'd=p  ==>  LL'd=-gf ==> Hd=-gf ==> d=-H^-1*gf is newton direction
lambda=sqrt(d'*H*d);
k=1;
while lambda^2/2>eta
    t=1;
    while checkfeasible(x+t*d)==0
        t=beta*t;
    end
    f1=feval(@fun,x+t*d);
    while f1>f+alpha*t*gf'*d
        t=beta*t;
        f1=feval(@fun,x+t*d);
    end
    x=x+t*d;
    k=k+1;
    f=f1;
    F(k)=f1;
    X(:,k)=x;
    T(k-1)=t;
    gf=feval(@gradf,x);
    H=feval(@hesf,x);
    Q=[H A';A zeros(p,p)] ;
    v=[gf; zeros(p,1)];
    e=-Q\v;
    d=e(1:n);
    lambda=sqrt(d'*H*d);

    %norm(gf)
end
k
xstar=x;
fstar=f;
gfstar=gf;