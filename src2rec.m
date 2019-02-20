function [path, itr]=src2rec(src_cord, rec_cord, pad, ttbl, nitr, fdOrder)
    srcx1=src_cord(1); srcx2=src_cord(2);
    recx1=rec_cord(1)+pad; 
    recx2=rec_cord(2)+pad;
    [n1, n2]=size(ttbl);
    dtdz=diff(ttbl, fdOrder, 1);
    dtdz=[dtdz; repmat(dtdz(n1-fdOrder,:), fdOrder, 1)];
    dtdx=diff(ttbl, fdOrder, 2);
    dtdx=[dtdx repmat(dtdx(:,n2-fdOrder), 1, fdOrder)];
    %[nz1, nx1]=size(dtdz);
    %[X, Z]=meshgrid(0:nx1-1, 0:nz1-1);
    %[Xq, Zq]=meshgrid(0:0.5:nx1-1,0:0.5:nz1-1);
    x=[recx1, recx2];
    path=zeros(nitr, 2);
    for itr=1:nitr-1
        path(itr, :)=x;
        x(1)=(x(1)>n1) * n1 + (~(x(1)>n1))*x(1);
        x(2)=(x(2)>n2) * n2 + (~(x(2)>n2))*x(2);
        ix1=int32(x(1)); ix2=int32(x(2));
        if x(2)-srcx2<2 
            break; 
        end
        %r=ttbl(ix1, ix2);
        g=[dtdz(ix1, ix2) dtdx(ix1, ix2)];
        %g=[interp2(X,Z,dtdz,x(2),x(1)), interp2(X,Z,dtdx,x(2),x(1))];
        %dx=g'*r;
        dx=g;
        %dr=g*dx;
        alpha=-2;%% fixed step length
        while(1)
            x0=x+alpha*dx/norm(dx);
            if(x0(1)<1 || x0(2)<1)
                alpha=alpha*0.5;
            else
                break;
            end
        end
        x=x+alpha*dx/norm(dx);
    end
    path(itr+1, :)=[srcx1 srcx2];
end