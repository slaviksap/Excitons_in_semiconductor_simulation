clear

A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonsIntensityAtBoundaryOverTime.txt')
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\photonsIntensityAtBoundaryOverTime.txt')

t1 = A1(:,1)
t2 = A2(:,1)
exc = A1(:,2)
phot = A2(:,2)

xlabel('t');
ylabel('intensity');
plot(t1,exc)

hold on

plot(t2,phot)

legend('exciton intensity','photon intensity')