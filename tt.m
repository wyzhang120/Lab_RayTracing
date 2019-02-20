function tt=tt(sln, n1, n2, h, srcx1, srcx2)
    % sln: float, array, slowness
    % n1: int, nz, size in fast/vertical dimension
    % n2: int, nx, size in slow/horizontal dimension    
    % srcx1: int, index of source coord in dimension1 / z
    % srcx2: int, index of source coord in dimension2 / x
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
    end
end