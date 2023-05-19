function [u_h,err] =solveNeumann(n,w,plot)
h = 1/(n+1);
u_h = 1;
%n=n+2 needed for the generation of T 
n = n+2;
X1 = 0:h:1;
Y1 = 0:h:1;
[xg,yg] = meshgrid(X1,Y1);
u = @(x,y) sin(2*pi*w*(x*x+y))+y;
lphu = @(x,y) 1/(h*h)*(-4*u(x,y)+u(x+h,y)+u(x-h,y)+u(x,y-h)+u(x,y+h));
fg = lphu(xg,yg);
fh = fg(:);
T = sparse(1:n,1:n,4*ones(1,n),n,n)+sparse(2:n,1:n-1,-1*ones(1,n-1),n,n)+sparse(1:n-1,2:n,-1*ones(1,n-1),n,n);
T(1,2)=-2;
T(n,n-1) = -2;
TT = -1*(T-4*eye(n,n));
%disp(full(T));
%disp(full(TT));
%create d for projection
dg = ones(n,n);
dg(:,1) = dg(:,1).*1/2;
dg(:,end) = dg(:,end)*1/2;
dg(1,:) = dg(1,:)*1/2;
dg(end,:) = dg(end,:)*1/2;

% todo nicer component wise mult with ones(1,n)*1/2 thus corners -> 1/4
d = dg(:);
Ah = 1/(u_h*u_h) * ((kron(eye(n,n),T))+kron(TT,eye(n,n)));
Ah = [Ah;ones(1,n*n)];
orthProject = fh-dot(d,fh)/norm(d,"inf")^2 * d;
orthProject = [orthProject;0];
u_h =  mldivide(Ah,orthProject);
%disp(u_h);
surf(reshape(u_h,n,n));
%disp(Ah);
