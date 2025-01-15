function feas=checkfeasible(x)
global  n
feas=1;
for i=1:n
    if x(i)<0
        feas=0;
        break
    end
end
