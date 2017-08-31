close all
clear all

load fbyfs.dat
load fbypp.dat

figure,plot(fbyfs(:,1),fbyfs(:,2));title('Drag without corrector');
figure,plot(fbyfs(:,1),fbyfs(:,3));title('Lift without corrector');
figure,plot(fbypp(:,1),fbypp(:,2));title('Drag with corrector');
figure,plot(fbypp(:,1),fbypp(:,3));title('Lift with corrector');