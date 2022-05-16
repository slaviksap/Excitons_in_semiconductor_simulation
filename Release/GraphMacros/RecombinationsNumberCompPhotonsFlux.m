hold on
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\photonsIntensityAtBoundaryOverTimeOneDoubleRecombination.txt')
A3 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\photonsIntensityAtBoundaryOverTimeMultyDoubleRecombinations.txt')



t2 = A2(:,1)
f2 = A2(:,2)

t3 = A3(:,1)
f3 = A3(:,2)

xlab = xlabel('t');
ylab = ylabel('Photons flux at boundary');

set(xlab,'fontsize',16)
set(ylab,'fontsize',16)

p2 = plot(t2,f2)
hold on
p3 = plot(t3,f3)
set(p1,'linewidth', 2)
set(p2,'linewidth', 2)
set(p3,'linewidth', 2)

legend('single double recombination possibility','multiple recombinations');
set(legend,'fontsize',16)