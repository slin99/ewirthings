function [u_h,err] =solveNeumann(n,w,plot)
h = 1/(n+1);
%n=n+2 needed for the generation of T 
n = n+2;
X1 = 0:h:1;
Y1 = 0:h:1;
omega = 2*pi*w;
grad_u = @(x, y) omega * [2 * x * cos(omega * (x* x + y)); cos(omega*(x * x + y)) + 1];
g = @(n1, n2, x, y) dot(grad_u(x, y), [n1, n2]);
[xg,yg] = meshgrid(X1,Y1);

u = @(x,y) sin(2.*pi.*w.*(x.*x+y))+y;
%lphu = @(x,y) 1/(h*h).*(-4.*u(x,y)+u(x+h,y)+u(x-h,y)+u(x,y-h)+u(x,y+h));
lphu = @(x,y) 1/h^2 * (4*pi*w.*(cos(2*pi*w.*(x.^2+y))-pi*w.*(4*x.^2+1).*sin(2*pi*w.*(x.^2+y))));
fg =lphu(xg,yg);
left = arrayfun(@(y) g(-1,0,0,y),yg);
right = arrayfun(@(y) g(1,0,1,y),yg);
top = arrayfun(@(x) g(0,1,x,1),xg);
bot = arrayfun(@(x) g(0,-1,x,0),xg);
lt_edge = g(-1/2,1/2,xg(1),yg(1))+fg(1,1);
rt_edge = g(1/2,1/2,xg(end),yg(1))+fg(1,end);
lb_edge = g(-1/2,-1/2,xg(1),yg(end))+fg(end,1);
rb_edge = g(1/2,-1/2,xg(end),yg(end))+fg(end,end);
fg(1,:) = left(1,:)+fg(1,:);
fg(end,:) = fg(end,:)+right(end,:);
fg(:,1) = fg(:,1)+top(:,1);
fg(:,end) = fg(:,end)+bot(:,end);
fg(1,1) = lt_edge(1,1);
fg(1,end) = lb_edge(1,end);
fg(end,end) = rb_edge(end,end);
fg(end,1) = rt_edge(end,1);
fh = fg(:);

%T = sparse(1:n,1:n,4*ones(1,n),n,n)+sparse(2:n,1:n-1,-1*ones(1,n-1),n,n)+sparse(1:n-1,2:n,-1*ones(1,n-1),n,n);
%T(1,2)=-2;
%T(n,n-1) = -2;
T = sparse(1:n,1:n,4.*ones(1,n),n,n)+sparse(2:n,1:n-1,[-1.*ones(1,n-2),-2],n,n)+sparse(1:n-1,2:n,[-2,-1.*ones(1,n-2)],n,n);

TT = 1*(T-4*speye(n,n));
%create d for projection
dg = ones(n,n);
dg(:,1) = dg(:,1)*1/2;
dg(:,end) = dg(:,end)*1/2;
dg(1,:) = dg(1,:)*1/2;
dg(end,:) = dg(end,:)*1/2;
d = dg(:);
Ah = ((kron(speye(n,n),T))+kron(TT,speye(n,n)));

Ah = [Ah;ones(1,n*n)];
orthProjector =  @(x) x-dot(d,x)/norm(d,"inf")^2 .* d;
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
disp(err)
if plot ~= "false"
    hold on
    tiledlayout(1,3)
    nexttile
    surf(X1,Y1,reshape(errV,n,n));
    title("error")
    nexttile
    %set(surf(X1,Y1,reshape(u_h,n,n)),"linestyle","none");
    surf(reshape(u_h,n,n))
    title("Solution")
    
end
