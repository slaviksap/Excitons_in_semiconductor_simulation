A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\checkSim.txt')
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\checkAn.txt')

t1 = A1(:,1)
sim = A1(:,2)
t2 = A2(:,1)
an = A2(:,2)

plot(t1,sim)
hold on
plot(t2,an)

legend('sim','an');