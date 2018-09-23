epsilon_0 = 8.85e-12;        % vacuum permittivity in F/m
k = 4*pi*epsilon_0;          % electrostatic const.
N = 100;

xrange = linspace(-1.5,1.5,N);
yrange = linspace(-0.5,2.5,N);
[Y,X] = meshgrid(yrange,xrange);
V1 = [];
V2 = [];

f1 = @(x,y) (integral(@(xp) (k*2*xp./sqrt((xp-x).^2+y.^2)),0,1));

f2 = @(x,y) (integral(@(xp) (k*xp.^2./sqrt((xp-x).^2+y.^2)),0,1) + ...
    integral(@(yp) (k*yp./sqrt(x.^2+(yp-y).^2)),1,2));

for i=1:N
    for j=1:N
        x = xrange(i);
        y = yrange(j);
        V1(i,j) = f1(x,y);
        V2(i,j) = f2(x,y);
    end
end

%%--------------- 1st Plot --------------------------------
[E1y,E1x] = gradient(-V1);
E1 = sqrt(E1x.^2 + E1y.^2);     % Magnitude of the field
E1xn = E1x./E1;                 % Normalization
E1yn = E1y./E1;                 % Normalization

figure
contour(X,Y,V1,50)              % Plot contour with 50 levels
xlabel('x')
ylabel('y')
title('Electric potential and field of a line charge')


Xs = X(1:5:end,1:5:end);        % Slice matrix to plot fewer arrows
Ys = Y(1:5:end,1:5:end);
E1xns = E1xn(1:5:end,1:5:end);
E1yns = E1yn(1:5:end,1:5:end);

hold on
quiver(Xs,Ys,E1xns,E1yns,0.5)   % Scale arrow lenght as half of its original
hold off


%%--------------- 2nd Plot --------------------------------
[E2y,E2x] = gradient(-V2);
E2 = sqrt(E2x.^2 + E2y.^2);
E2xn = E2x./E2;
E2yn = E2y./E2;

figure
contour(X,Y,V2,50)
xlabel('x')
ylabel('y')
title('Electric potential and field of a L-shaped charge')

E2xns = E2xn(1:5:end,1:5:end);
E2yns = E2yn(1:5:end,1:5:end);

hold on
quiver(Xs,Ys,E2xns,E2yns,0.5)
hold off


%%--------------- 3rd Plot --------------------------------
xrange = linspace(-3,3,N);       %Plot in a different region
yrange = linspace(-3,3,N);
[Y,X] = meshgrid(yrange,xrange);

V3 = [];

for i=1:N
    for j=1:N
        x = xrange(i);
        y = yrange(j);
        
        f3_rect = @(xp,yp) (k*xp./sqrt((xp-x).^2+(yp-y).^2));           %define function in rect coordinates
        f3_polar = @(r,theta) f3_rect(r.*cos(theta),r.*sin(theta)).*r;  %convert to polar coordinates
        V3(i,j) = integral2(f3_polar,0,2,0,2*pi);                       %Evaluate double integral
    end
end

[E3y,E3x] = gradient(-V3);
E3 = sqrt(E3x.^2 + E3y.^2);
E3xn = E3x./E3;
E3yn = E3y./E3;

figure
contour(X,Y,V3,50)
xlabel('x')
ylabel('y')
title('Electric potential and field of a charged disc')

Xs = X(1:5:end,1:5:end);
Ys = Y(1:5:end,1:5:end);
E3xns = E3xn(1:5:end,1:5:end);
E3yns = E3yn(1:5:end,1:5:end);

hold on
quiver(Xs,Ys,E3xns,E3yns,0.5)
hold off