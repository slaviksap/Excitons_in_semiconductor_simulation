hold on
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\equationSystemCheckWithIntegration.txt')

t = A1(:,1)
f = A1(:,2)

xlabel('t');
ylabel('F');

e = exp(-2*t)
plot(t,f)
hold on
plot (t,e)
