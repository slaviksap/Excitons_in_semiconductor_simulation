hold on
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonsIntensityAtBoundaryOverTime.txt')

t1 = A1(:,1)
exc = A1(:,2)

xlabel('t');
ylabel('intensity');
plot(t1,exc)
