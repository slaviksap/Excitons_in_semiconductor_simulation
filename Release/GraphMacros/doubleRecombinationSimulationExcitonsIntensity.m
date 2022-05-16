hold on
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonsIntensityAtBoundaryOverTime.txt')

t1 = A1(:,1)
exc = A1(:,2)

xlabel('t');
ylabel('intensity');
p = plot(t1,exc)
set(p,'linewidth', 2)
legend('small spheres simulation r = H/20','max spheres simulation','max spheres with correction')
set(legend,'fontsize',18)
