A1 = load('C:\\Projects\\Excitons_in_semiconductor_simulation\\Release\\excitonsIntensityAtBoundaryHitMap.txt')
x = A1(:,1)
y = A1(:,2)
z = A1(:,3)

n = 200;
[X, Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z = griddata(x,y,z,X,Y);
%// Remove the NaNs for imshow:
Z(isnan(Z)) = 0;
m = min(Z(Z~=0));
M = max(Z(Z~=0));
imshow((Z-m)/(M-m));
