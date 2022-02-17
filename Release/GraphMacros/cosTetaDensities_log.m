clear
A2 = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\cos tetta density from non-central point(simplified).txt');
A1 = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\cos tetta density from non-central point.txt');
x1 = A1(:,1);
y1 = A1(:,2);

x2 = A2(:,1);
y2 = A2(:,2);
semilogy(x1,y1)
hold on
semilogy(x2,y2)
legend('cosTetta density','cosTetta density simplified')

parametres = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\parametres.txt')
initCond = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\initial conditions.txt')
title("v= " + parametres(2) + ", D = " + parametres(3) + ", tau = " + parametres(4) + ", R = " + initCond(1) + ", r = " + initCond(2) + "  logaxis")