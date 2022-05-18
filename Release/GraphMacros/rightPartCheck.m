hold on
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\equationSystemCheckWithIntegration.txt')

t = A1(:,1)
f = A1(:,2)

xlabel('t');
ylabel('F');

e = t/2
p1 = plot(t,f)
set(p1,'linewidth', 2)

hold on
p2 = plot(t,e)
set(p2,'linewidth', 2)

legend('simulated right part','1-t\2');
set(legend,'fontsize',16)