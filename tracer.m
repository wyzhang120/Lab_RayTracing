% A ray tracor function
% inputs:
%   vel    (2D) - velocity model
%   dx, dz (#)  - sampling interval in x & z directions
%   dir (2VEC)  - take off direction 
%   
% outputs:
%   x, y (VEC) - positions along the ray
%   s    (VEC) - distances from the source point
%   t    (VEC) - travel-time distances from the source
function [x,z,s,t]=tracer(vel, dx, dz, dir, srcpos, scount)
    x=zeros(scount,1);
    z=x; s=x; t=x;
    tau=0;
    % Computing partial derivatives of velocity fields
    [N1,N2]=size(vel);
    vx=zeros(N1,N2); % x-vel gradient
    vz=zeros(N1,N2); % z-vel gradient
    for m=1:N2-1
        vx(:,m)=(-vel(:,m)+vel(:,m+1))/dx;
    end
    vx(:,N2)=vx(:,N2-1); % extrapolating the gradient for edge
    vx=-vx./vel.^2;
    for m=1:N1-1
        vz(m,:)=(-vel(m,:)+vel(m+1,:))/dz;
    end
    vz(N1,:)=vz(N1-1,:); % extrapolating the gradient for edge
    vz=-vz./vel.^2;
    p=dir/norm(dir)/vel(srcpos(1),srcpos(2));
    curpos=srcpos.*[dz;dx];
    
    for is=1:scount
        ipos=int32(curpos.*[1/dz; 1/dx]+1);
        iz=ipos(1);
        ix=ipos(2);
        if(ix>0 && ix<= N2 && iz>0 && iz<= N1)
            dpz=vz(iz,ix);
            dpx=vx(iz,ix);
            vs=vel(iz,ix);
            dxzds=vs*p;
            tau=tau+1/vs;
            
            curpos=curpos+dxzds;
            p=p+[dpz;dpx];
            
            z(is)=curpos(1)/dz;
            x(is)=curpos(2)/dx;
            s(is)=is;
            t(is)=tau;
        else
            z(is:scount)=z(is-1);
            x(is:scount)=x(is-1);
            s(is:scount)=s(is-1);
            t(is:scount)=t(is-1);
            break;
        end
    end