ga = Grid_Axis();
Dnons = ga.get_derivative_operator(2);
Dstag = ga.get_derivative_operator(2,1);

xn = ga.node;
xm = ga.mid;

p4 = 4*xn.^2 + 2*xn;
p3 = 8*xm + 2*ones(size(p4));
plot(Dstag*p4 - p3)


g = FD_Grid(); g

gxn = g.grid_axis{1}.node;
gxm = g.grid_axis{1}.mid;
gyn = g.grid_axis{2}.node;
gym = g.grid_axis{2}.mid;

[X,Y] = meshgrid(gxn,gyn);

Z = ones(32,32);%2*X+Y;
Z = X.^2 + Y.^2;
Z = 2*X+Y;
surf(Z)

surf(reshape(g.Lap*Z(:),32,32))

g.grid_axis{1}.set('N',64)
g.rebuild();

gxn = g.grid_axis{1}.node;
gxm = g.grid_axis{1}.mid;
gyn = g.grid_axis{2}.node;
gym = g.grid_axis{2}.mid;

[X,Y] = meshgrid(gxn,gyn);

Z = ones(32,32);%2*X+Y;
Z = X.^2 + Y.^2;
Z = 2*X'+Y';
surf(Z)

Z1 = g.Grad*Z(:);
figure
surf(reshape(Z1(1:64*32),64,32))
figure
surf(reshape(Z1(64*32+1:end),64,32))