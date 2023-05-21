function [u_h,err] =solveNeumann(n,w,plot)

function res = stencil(u,xg,yg,h)
    u_eval = u(xg,yg);
    res = u_eval;
 %   grad_u = @(x, y) omega * [2 * x * cos(omega * (x* x + y)); cos(omega*(x * x + y)) + 1];
%g = @(n1, n2, x, y) dot(grad_u(x, y), [n1, n2]);
    %cool tenary operator
    tenary = @(varargin)varargin{end-varargin{1}};
    % does not say no the corners
    isEdge = @(x,y) x==length(u_eval) ||x == 1 || y == length(u_eval) || y == 1;
    isCorner = @(x,y) (x== length(u_eval) || x== 1) && ( y==1 || y == length(u_eval));

    %fuck efficiency
    for i= 1:length(u_eval) 
        for j= 1:length(u_eval)
            r = 0;
            if isCorner(i,j)
                r =  1/h^2 .*(4*u_eval(i,j)-2*u_eval(tenary(i-1==0,i+1,i-1),j)-2*u_eval(i,tenary(j-1==0,j+1,j-1)));
            %edge case should handle corner case 
            elseif isEdge(i,j)
                 r = 1/h^2  .*(4*u_eval(i,j)-u_eval(tenary(i-1==0,i+1,i-1),j)-u_eval(i,tenary(j-1==0,j+1,j-1))-u_eval(tenary(i==length(u_eval),i-1,i+1),j)-u_eval(i,tenary(j==length(u_eval),j-1,j+1)));
            else
                %we need this case since we need to filp the signs
             r =  1/h^2 .*(-4.*u_eval(i,j)+u_eval(i-1,j)+u_eval(i,j-1)+u_eval(i+1,j)+u_eval(i,j+1));
            end




            res(i,j)=r;
        end
    end
end
h = 1/(n+1);
%n=n+2 needed for the generation of T 
n = n+2;
X1 = 0:h:1;
Y1 = 0:h:1;
omega = 2*pi*w;

[xg,yg] = meshgrid(X1,Y1);

u = @(x,y) sin(2.*pi.*w.*(x.*x+y))+y;

fg = stencil(u,xg,yg,h);
%disp(fg);
%fg = fg';
fg = flip(fg)';
fh = fg(:);
%T = sparse(1:n,1:n,4*ones(1,n),n,n)+sparse(2:n,1:n-1,-1*ones(1,n-1),n,n)+sparse(1:n-1,2:n,-1*ones(1,n-1),n,n);
%T(1,2)=-2;
%T(n,n-1) = -2;

T = sparse(1:n,1:n,4.*ones(1,n),n,n)+sparse(2:n,1:n-1,[-1.*ones(1,n-2),-2],n,n)+sparse(1:n-1,2:n,[-2,-1.*ones(1,n-2)],n,n);

TT = (T-4*speye(n,n));
%create d for projection
dg = ones(n,n);
dg(:,1) = dg(:,1)*1/2;
dg(:,end) = dg(:,end)*1/2;
dg(1,:) = dg(1,:)*1/2;
dg(end,:) = dg(end,:)*1/2;
d = flip(dg)';
d = d(:);
%calc Ah
Ah =1/h^2 * ((kron(speye(n,n),T))+kron(TT,speye(n,n)));
Bh = 1/h^2 * gallery("neumann",n^2);
disp("diff");
%disp([size(Ah),size(Bh)]);
%disp(full(Ah));
%disp(full(Bh));
disp(norm(Ah-Bh,"inf"));
% mittelwert freiheit
Ah = [Ah;ones(1,n*n)];
%disp(full(T))
%calc projector 
orthProjector =  @(x) x-(dot(d,x)/norm(d,"inf")^2) .* d;
orthProject=  orthProjector(fh);
%mittelwert freiheit
orthProject = [orthProject;0];
u_h =  mldivide(Ah,orthProject);
%now we do error estimation
uGrid= u(xg,yg);
uGrid = flip(uGrid)';
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
    surf(X1,Y1,reshape(uGridVector,n,n));
    title("error")
    nexttile
    %set(surf(X1,Y1,reshape(u_h,n,n)),"linestyle","none");
    surf(flip(reshape(u_h,n,n)'))
    title("Solution")
    
end
end
