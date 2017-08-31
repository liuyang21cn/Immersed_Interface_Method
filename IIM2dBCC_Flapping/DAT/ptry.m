m=100;
n=2*m+1;
x=[-m:m]/m;
y=x;
p=zeros(n,n);
ratio=1.25;

for j=1:n
    for i=1:n
        r2=x(i)*x(i)+y(j)*y(j);
        if r2<0.25
            p(i,j)=ratio*(x(i)+y(j))+r2;
        end
    end
end

figure(1)
contour(x,y,p,60);
axis equal;
            
