hold on
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\photonsIntensityAtBoundaryOverTime.txt')

t2 = A2(:,1)
phot = A2(:,2)

xlabel('t');
ylabel('intensity');
plot(t2,phot)
