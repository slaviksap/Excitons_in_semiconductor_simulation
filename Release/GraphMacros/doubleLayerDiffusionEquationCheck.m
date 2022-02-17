clear
A1 = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\analitGraph.txt')
A2 = load('C:\\Projects\\Excitons in semiconductor simulation\\Release\\calcGraph.txt')

x = A1(:,1)
y1 = A1(:,2)
y2 = A2(:,2)

plot(x,y1)
hold on
plot(x,y2)

legend('analytical solution','simulated solution ')