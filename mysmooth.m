function ttbl=mysmooth(ttt1, nitr)
ttbl=ttt1;
[n1, n2]=size(ttt1);
filt=[1/sqrt(2) 1 1/sqrt(2); 1 2 1; 1/sqrt(2) 1 1/sqrt(2)];
filt=filt/sum(sum(filt));
if nitr==0
    ttbl=ttt1;
else
    for itr=1:nitr
        tmp=ttbl;
        for i2=2:n2-1
            for i1=2:n1-1
                ttbl(i1, i2)=sum(sum(filt.*tmp(i1-1:i1+1, i2-1:i2+1)));
            end
            ttbl(1, i2)=(2*tmp(1,i2)+tmp(2,i2))/3;
            ttbl(n1, i2)=(2*tmp(n1,i2)+tmp(n1-1,i2))/3;
        end

        if(itr<nitr-2)
            %ttbl(srcx1, srcx2) = ttt2(srcx1, srcx2);
        end
    end
end
end