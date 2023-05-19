function [u_h,err] =solveNeumann(n,w,plot)
h = 1/(n+1);
%n=n+2 needed for the generation of T 
n = n+2;
X1 = 0:h:1;
Y1 = 0:h:1;
[xg,yg] = meshgrid(X1,Y1);
u = @(x,y) sin(2*pi*w*(x*x+y))+y;
lphu = @(x,y) 1/(h*h)*(-4*u(x,y)+u(x+h,y)+u(x-h,y)+u(x,y-h)+u(x,y+h));
%lphu = @(x,y) 1/h *(4*pi*w*(cos(2*pi*w*(x^2+y))-pi*w*(4*x^2+1)*sin(2*pi*w*(x^2+y))));
fg = lphu(xg,yg);
fh = fg(:);
%T = sparse(1:n,1:n,4*ones(1,n),n,n)+sparse(2:n,1:n-1,-1*ones(1,n-1),n,n)+sparse(1:n-1,2:n,-1*ones(1,n-1),n,n);
%T(1,2)=-2;
%T(n,n-1) = -2;
T = sparse(1:n,1:n,4*ones(1,n),n,n)+sparse(2:n,1:n-1,[-1*ones(1,n-2),-2],n,n)+sparse(1:n-1,2:n,[-2,-1*ones(1,n-2)],n,n);

TT = 1*(T-4*speye(n,n));
%create d for projection
dg = ones(n,n);
dg(:,1) = dg(:,1)*1/2;
dg(:,end) = dg(:,end)*1/2;
dg(1,:) = dg(1,:)*1/2;
dg(end,:) = dg(end,:)*1/2;

% todo nicer component wise mult with ones(1,n)*1/2 thus corners -> 1/4
d = dg(:);
Ah = 1/(h*h) * ((kron(speye(n,n),T))+kron(TT,speye(n,n)));
Ah = [Ah;ones(1,n*n)];
orthProjector =  @(x) x-dot(d,x)/norm(d,"inf")^2 * d;
orthProject=  orthProjector(fh);
orthProject = [orthProject;0];
u_h =  mldivide(Ah,orthProject);
%now we do error estimation
uGrid= u(xg,yg);
uGridVector= uGrid(:);
errV = orthProjector(uGridVector)-u_h;
err = norm(errV,"inf");
%disp(u_h);
%disp(Ah);

if plot ~= "false"
    hold on
    tiledlayout(1,2)
    nexttile
    surf(reshape(errV,n,n));
    title("error")
    nexttile
    surf(reshape(u_h,n,n));
    title("Solution")
end
