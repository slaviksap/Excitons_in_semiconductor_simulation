A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\photonsIntensityAtBoundaryOverTime.txt')
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\photonsIntensityAtBoundaryOverTimeWithCylinderDislocation.txt')

t1 = A1(:,1)
phot1 = A1(:,2)
t2 = A2(:,1)
phot2 = A2(:,2)


xlabel('t');
ylabel('intensity');

plot(t1,phot1)
hold on
plot(t2,phot2)

legend('Photon intensiny','Photon intensity with cylindrical dislocation')