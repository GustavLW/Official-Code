clc
hold off
w = [0.030 0.100 0.005];
s = [0.030 0.150 0.035];
n = 100;
r = 1;
x = linspace(0,1,n);
yw = birth_rate(x,w);
ys = birth_rate(x,s);
plot(x,yw,'b')
hold on
plot(x,ys,'r')
plot([0 1],[0 0],'k')
grid on
axis([0 1/r sqrt(r)*max(min(yw(1:n/r)),min(ys(1:n/r))) 1.1*max(max(yw(1:n/r)),max(ys(1:n/r)))])
function y = birth_rate(x,p)
    y = x.*((p(1) + p(2).*x).*(1-x) - p(3));
end