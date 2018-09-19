load('derdata.mat');

figure
plot(X,Y,'-o','MarkerSize',5)
title('Y v.s. X')

%Calculate the first derivative of (X,Y)
[fdy, fdx] = Der(Y,X,'fd');
[bdy, bdx] = Der(Y,X,'bd');
[cdy, cdx] = Der(Y,X,'cd');
[hody, hodx] = Der(Y,X,'hod');

figure
plot(cdx,cdy,'k','LineWidth',1)
title('1st derivative')

figure
plot(cdx,fdy-cdy,'b','LineWidth',1)
title('Compare forwardiff and centraldiff')
figure
plot(cdx,bdy-cdy,'g','LineWidth',1)
title('Compare backdiff and centraldiff')
figure
plot(cdx,hody-cdy,'r','LineWidth',1)
title('Compare higherorderdiff and centraldiff')

%Calculate the second derivative of (X,Y)
[cddy, cddx] = Der(cdy,cdx,'cd');

figure
plot(cddx,cddy,'k','LineWidth',1)
title('2nd derivative')


