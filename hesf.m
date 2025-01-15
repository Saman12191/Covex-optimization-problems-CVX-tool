function H=hesf(x)
global A p n
for i=1:n
H(i,i)=1/x(i);
end