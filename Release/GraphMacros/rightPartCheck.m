hold on
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\equationRightPart.txt')

z = A1(:,1)
t = A1(:,2)
f = A1(:,3)

xlabel('z');
ylabel('t');
zlabel('\Phi');

e = 0.05*0.05*exp(-0.05*t)
plot(z,t,f)
