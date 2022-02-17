hold on
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonGrid.txt')
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonGrid2.txt')

z = A1(:,1)
t = A1(:,2)
f = A1(:,3)

z2 = A2(:,1)
t2 = A2(:,2)
f2 = A2(:,3)


xlabel('z');
ylabel('t');
zlabel('\Phi');


plot3(z,t,f)
hold on
plot3(z2,t2,f2)
