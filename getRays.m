function [rays_z, rays_x]=getRays(sln, nz, nx,  ns, nr, s_z, r_x,  h)
rays_z=zeros(length(ns), length(nr), nz); % geom: src, rec, seg
rays_x=zeros(length(ns), length(nr), nz);
dx=0.1; % spacing for calculating derivatives
%%
parfor is_x=1:length(ns)
    s_x=ns(is_x);
    temp_rays_z=zeros(length(nr), nz);
    temp_rays_x=zeros(length(nr), nz);
    for ir_z=1:length(nr)
        r_z=nr(ir_z);
        ray_z=s_z:sign(r_z-s_z):r_z;
        ray_x=s_x+(ray_z-s_z)*(r_x-s_x)/(length(ray_z)-1);
        g=zeros(1, length(ray_z));
        dir=g*0;
        H=zeros(1, length(g));
        
        for itr=1:2000%min(length(ray_z)^2, 1000)
            for iseg=2:length(ray_z)-1
                % estimating the first derivative with respect to current x
                % position of segment end
                fx=         lineIntegral(sln, ray_x(iseg-1), ray_z(iseg-1),ray_x(iseg), ray_z(iseg), nx, h);
                fx=fx+      lineIntegral(sln, ray_x(iseg+1), ray_z(iseg+1),ray_x(iseg), ray_z(iseg), nx, h);
                % change in travel time due to +dx shifting
                fxp1=-fx+   lineIntegral(sln, ray_x(iseg-1), ray_z(iseg-1),ray_x(iseg)+dx, ray_z(iseg), nx,h);
                fxp1=fxp1+  lineIntegral(sln, ray_x(iseg+1), ray_z(iseg+1),ray_x(iseg)+dx, ray_z(iseg), nx,h);
                %fxp1=fxp1+10000000*(ray_x(iseg)+dx>=nx);
                % change in travel time due to -dx shifting
                fxm1=-fx+   lineIntegral(sln, ray_x(iseg-1), ray_z(iseg-1),ray_x(iseg)-dx, ray_z(iseg), nx,h);
                fxm1=fxm1+  lineIntegral(sln, ray_x(iseg+1), ray_z(iseg+1),ray_x(iseg)-dx, ray_z(iseg), nx,h);
                %fxm1=fxm1+10000000*(ray_x(iseg)+dx<=1);
                g(iseg)=(fxp1-fxm1)/(2*dx);
                H(iseg)=(fxp1+fxm1-2*fx)/(12*dx);
            end
            if(norm(g)<1e-5); break; end;
            dir(2:length(g)-1)=1./H(2:length(g)-1)'.*g(2:length(g)-1)';
            ray_x=ray_x+dir;
            %ray_x=ray_x+(r_x-ray_x).*(ray_x<r_x);
            if (itr<80 && mod(itr,10)==0)
                tempray_x=ray_x;
                for pass=1:5
                 for iseg=2:length(ray_z)-1
                     tempray_x(iseg)=(ray_x(iseg-1)+ray_x(iseg+1)+ray_x(iseg))/3.0;
                 end
                end
                 ray_x=tempray_x;
            end
        end
        nseg=length(ray_z);
        temp_rays_z( ir_z, 1:nseg)=ray_z;
        temp_rays_x( ir_z, 1:nseg)=ray_x;
    end
    rays_z(is_x, :, :)=temp_rays_z(:,:);
    rays_x(is_x, :, :)=temp_rays_x(:,:);
end

end