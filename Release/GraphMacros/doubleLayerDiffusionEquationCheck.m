clear
A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\analitGraph.txt')
A2 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\calcGraph.txt')

x = A1(:,1)
y1 = A1(:,2)
y2 = A2(:,2)

plot(x,y1,"linewidth", 3)
hold on
plot(x,y2,"linewidth", 2)

legend('analytical solution','simulated solution ')