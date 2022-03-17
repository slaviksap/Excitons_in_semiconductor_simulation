hold on
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonsNumberAtBoundaryOverTimeWithMultipleDislocations.txt')

t2 = A2(:,1)
phot2 = A2(:,2)


xlabel('t');
ylabel('Relative number of particles');

plot(t2,phot2)

