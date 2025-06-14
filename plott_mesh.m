clc;
clear all;
close all;


nel = length(v)/4;

figure;
hold on;
for i=1:nel
    
    xx = [v(4*(i-1)+1,1),v(4*(i-1)+2,1),v(4*(i-1)+3,1),v(4*(i-1)+4,1),v(4*(i-1)+1,1)];
    yy = [v(4*(i-1)+1,2),v(4*(i-1)+2,2),v(4*(i-1)+3,2),v(4*(i-1)+4,2),v(4*(i-1)+1,2)];
    plot(xx,yy,'-r');
    pause;
    
end