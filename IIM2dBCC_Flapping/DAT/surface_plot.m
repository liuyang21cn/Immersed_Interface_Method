clear all

theta=180/3.1415926;

ll=1;
if ll==1
    load surface1.dat
    f=surface1;
elseif ll==2
    load surface2.dat
    f=surface2;
elseif ll==3
    load surface3.dat
    f=surface3;
end

figure(1)
plot(f(:,1),f(:,2))
xlabel('xs')
ylabel('ys')
axis equal


% figure(3)
% plot(theta*f(:,1),f(:,5),theta*f(:,1),f(:,7),'o')
% xlabel('alpha')
% ylabel('vs')
% 
% figure(4)
% plot(theta*f(:,1),f(:,8),'-',theta*f(:,1),f(:,9),'--')
% xlabel('alpha')
% ylabel('ft-,fn--')
%  
% load falfa.dat;
% ip=1+3;
% 
% figure(5)
% plot(theta*falfa(:,1),falfa(:,ip),'-')
% xlabel('alpha')
% ylabel('falpha')


