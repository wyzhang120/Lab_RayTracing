function f=lineIntegral(fld, s_x, s_z, r_x,  r_z, nx, h)
if(r_x<1 || s_x<1|| r_x>nx || s_x>nx); f=10000; return; end;
% linear integration of a field
% assumes uniform spacing horizontally and vertically
p=[(r_z-s_z); (r_x-s_x)]; p=p/norm(p);
cur_pos=[s_z; s_x];
%hold on;
rem=[mod(-s_z,sign(r_z-s_z));mod(-s_x,sign(r_x-s_x))]; % remainder to the edge of the cell
t=0;
while(norm(cur_pos-[r_z; r_x])>1e-4)
%for q=1:10
    if(p(1)~=0); ratio1=rem(1)/p(1); else ratio1=1000; end;
    if(p(2)~=0); ratio2= rem(2)/p(2); else ratio2=1000; end;
    if(ratio1<ratio2)
        l=ratio1;
        rem(1)=sign(r_z-s_z);
        rem(2)=rem(2)-l*p(2);
    else
        l=ratio2;
        rem(2)=sign(r_x-s_x);
        rem(1)=rem(1)-l*p(1);
    end
    lpos=cur_pos;
    cur_pos=cur_pos+l*p;
    t=t+l*h*fld(floor((cur_pos(1)+lpos(1))/2),floor((cur_pos(2)+lpos(2))/2));
    %scatter(cur_pos(2), cur_pos(1));%pause(.001);
end
%scatter(cur_pos(2), cur_pos(1));pause(.1);
%hold off;
f=t;
end