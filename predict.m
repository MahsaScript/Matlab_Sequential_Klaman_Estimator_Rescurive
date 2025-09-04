function yp=predict(tb,yb,n,f_unk,B_un,mass);

% functions of the generalized state equations of a Linear System

g(1:n,1)=yb(n+1:2*n);
[stiff,damp,xkp,xcp]=kcm(n,yb);
g(n+1:2*n,1)=B_un*f_unk+inv(mass)*(-stiff*yb(1:10)-damp*yb(10+1:2*10));
g(2*n+1:33,1)=zeros(13,1);

yp =g;

   