%%
% asymptotic Green's function
function [tt, A]=asympgreens(sln, n1, n2, h, srcx1, srcx2)
    
    tt=sln*0;
    %tt(srcx1, srcx2)=0;
    tt(srcx1+1, srcx2)=h*(sln(srcx1+1, srcx2)+sln(srcx1, srcx2))/2;
    tt(srcx1-1, srcx2)=h*(sln(srcx1-1, srcx2)+sln(srcx1, srcx2))/2;
    tt(srcx1, srcx2+1)=h*(sln(srcx1, srcx2+1)+sln(srcx1, srcx2))/2;
    tt(srcx1, srcx2-1)=h*(sln(srcx1, srcx2-1)+sln(srcx1, srcx2))/2;
    K=tt*0;
    %K(srcx1-1:srcx1+1, srcx2)=1;
    %K(srcx1, srcx2-1:srcx2+1)=1;
    %%
    for itr=1:3*sqrt(n1^2+n1^2)
        tt2=tt;
        if(mod(itr,5)>0)
            for i2=1:n2
            for i1=1:n1
                if(abs(i1-srcx1)>0 || abs(i2-srcx2)>0)
                    if(i1==1)
                        tmin1=tt(i1+1,i2);
                    elseif(i1==n1)
                        tmin1=tt(i1-1,i2);
                    else
                        tmin1=min([tt(i1-1,i2) tt(i1+1,i2)]);
                    end
                    if(i2==1)
                        tmin2=tt(i1,i2+1);
                    elseif(i2==n2)
                        tmin2=tt(i1,i2-1);
                    else
                        tmin2=min([tt(i1,i2-1) tt(i1,i2+1)]);
                    end
                    a=2;
                    b=-2*(tmin1+tmin2);
                    c=tmin1^2+tmin2^2-h^2*sln(i1, i2)^2;
                    d=b^2-4*a*c;            
                    if(d>=0)
                        tt2(i1, i2)=(-b+sqrt(d))/(2*a);
                        if(abs(tt2(i1,i2)-tt(i1,i2))<.0001)
                            K(i1,i2)=1;
                        end
                    end

                end
            end
            end
        else
            for i2=1:n2
            for i1=1:n1
                if(abs(i1-srcx1)>1 || abs(i2-srcx2)>1)
                    if(i1==1)
                        if(i2<n2)
                            tmin1=tt(i1+1,i2+1);
                        else
                            tmin1=tt(i1+1,i2-1);%rvrslogic
                        end
                    elseif(i1==n1)
                        if(i2>1)
                            tmin1=tt(i1-1,i2-1);
                        else
                            tmin1=tt(i1-1,i2+11);%rvrslogic
                        end
                    else
                        if(i2==1)
                            tmin1=tt(i1+1,i2+1);
                        elseif(i2==n2)
                            tmin1=tt(i1-1,i2-1);
                        else
                            tmin1=min([tt(i1-1,i2-1) tt(i1+1,i2+1)]);
                        end
                    end
                    
                    if(i2==1)
                        if(i1>1)
                            tmin2=tt(i1-1,i2+1);
                        else
                            tmin2=tt(i1+1,i2+1);
                        end
                    elseif(i2==n2)
                        if(i1<n1)
                            tmin2=tt(i1+1,i2-1);
                        else
                            tmin2=tt(i1-1,i2-1);
                        end
                    else
                        if(i1==1)
                            tmin2=tt(i1+1,i2-1) ;
                        elseif (i1==n1)
                            tmin2=tt(i1-1,i2-1) ;  
                        else
                            tmin2=min([tt(i1+1,i2-1) tt(i1-1,i2+1)]);
                        end
                    end
                    
                    
                    a=2;
                    b=-2*(tmin1+tmin2);
                    c=tmin1^2+tmin2^2-(sqrt(2)*h)^2*sln(i1, i2)^2;
                    d=b^2-4*a*c;            
                    if(d>=0)
                        tt2(i1, i2)=(-b+sqrt(d))/(2*a);
                    end

                end
            end
            end
        end
        tt=tt2;
        if(mod(itr,25)==0)
            imagesc(tt); title(itr); colorbar; 
            %caxis([0 1]);
            pause(.01);
        end
    end
    %% smoothing for distance computation
    ttbl=tt;
    nitr=5;
    filt=[sqrt(2) 1 sqrt(2); 1 2 1; sqrt(2) 1 sqrt(2)];
    filt=filt/sum(sum(filt));
    for itr=1:nitr
        tmp=ttbl;
        for i2=2:n2-1
            for i1=2:n1-1
                ttbl(i1, i2)=sum(sum(filt.*tmp(i1-1:i1+1, i2-1:i2+1)));
            end
            ttbl(1, i2)=(2*tmp(1,i2)+tmp(2,i2))/3;
            ttbl(n1, i2)=(2*tmp(n1,i2)+tmp(n1-1,i2))/3;
        end
    end
    imagesc(ttbl);
    %%
    dtdz=diff(ttbl, 1, 1);
    dtdz=[dtdz; dtdz(n1-1,:)];
    dtdx=diff(ttbl, 1, 2);
    dtdx=[dtdx dtdx(:,n2-1)];
    nitr=500;
    A=ttbl*0;
    for recx2=1:n2;
    for recx1=1:n1;
        %[recx1 recx2]
        x=[recx1; recx2];
        path=zeros(nitr, 2);
        for itr=1:nitr
            path(itr, :)=x;
            ix1=int32(x(1));
            ix2=int32(x(2));
            if(abs(ix1-srcx1)<2 && abs(ix2-srcx2)<2); break; end;
            r=ttbl(ix1, ix2);
            g=[dtdz(ix1, ix2) dtdx(ix1, ix2)];
            dx=g'*r;
            %dr=g*dx;
            alpha=-2;%%-(r*dr)/(dr*dr)
            while(1)
                x0=x+alpha*dx/norm(dx);
                if(x0(1)<1 || x0(2)<1 || x0(1)>n1 || x0(2)>n2); 
                    alpha=alpha*0.5; 
                else
                    break;
                end;
            end
            x=x+alpha*dx/norm(dx);
        end
        %
        if(mod(recx1, 10)==0)
            imagesc(tt);
            hold on;
            plot(path(1:itr-1,2), path(1:itr-1,1), 'k');
            hold off;
            pause(.01)
        end
        %% distance calculations
        dist=0;
        for is=1:itr-2
            dist=dist+norm(path(is,:)-path(is+1,:));
        end
        A(recx1, recx2)=(dist+2)*10;
    end
    end
    %% smoothing calculated results
    aa=A;
    nitr=3;
    filt=[sqrt(2) 1 sqrt(2); 1 2 1; sqrt(2) 1 sqrt(2)];
    filt=filt/sum(sum(filt));
    for itr=1:nitr
        tmp=aa;
        for i2=2:n2-1
            for i1=2:n1-1
                aa(i1, i2)=sum(sum(filt.*tmp(i1-1:i1+1, i2-1:i2+1)));
            end
            aa(1, i2)=(2*tmp(1,i2)+tmp(2,i2))/3;
            aa(n1, i2)=(2*tmp(n1,i2)+tmp(n1-1,i2))/3;
        end
        aa(srcx1, srcx2)=tmp(srcx1, srcx2);
    end
    imagesc(aa);
    A=aa;
end