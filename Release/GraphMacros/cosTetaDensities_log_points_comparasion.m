clear
A1 = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\cos tetta density from non-central point.txt');
x1 = A1(:,1);
y1 = A1(:,2);
plot(x1,y1)
parametres = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\parametres.txt')
initCond = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\initial conditions.txt')

v = 7
k = v*initCond(1)/(2*parametres(3))
x2 = 0:0.01:2
y2 = k/(1 - exp(-2*k))*exp(-k*x2)

hold on
plot(x2,y2)

v = 6
k = v*initCond(1)/(2*parametres(3))
x2 = 0:0.01:2
y2 = k/(1 - exp(-2*k))*exp(-k*x2)

hold on
plot(x2,y2)

v = 5
k = v*initCond(1)/(2*parametres(3))
x2 = 0:0.01:2
y2 = k/(1 - exp(-2*k))*exp(-k*x2)

hold on
plot(x2,y2)

v = 4
k = v*initCond(1)/(2*parametres(3))
x2 = 0:0.01:2
y2 = k/(1 - exp(-2*k))*exp(-k*x2)

hold on
plot(x2,y2)
legend('1-cosTetta density from non central point without drift velocity','1-cosTetta density from central point, v = 7','1-cosTetta density from central point, v = 6','1-cosTetta density from central point, v = 5','1-cosTetta density from central point, v = 4')