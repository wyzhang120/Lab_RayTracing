path2model = 'E:/Geophysics/Project/Crosswell/FWI_2arr/vp22_elp';
%path2model = 'C:\DFiles\Geophysics\Project\Figs_Crosswell/vp22_elp';
ncpu=2;
nz=106; nx=301; dx=1; dz=1;
fid = fopen(path2model, 'r');
vp = fread(fid, [nz, nx], 'float32')/1000;
% vp2=vp;
% mask2=vp==4.06;
% mask1=vp==5.8;
% vp2(mask2)=5.8;
% vp2(mask1)=4.06;
% vp=vp2;
%% 
srcz=50; % source z position
srcx=0;% source x position
h=dx;
% pad vp model for time table calculation
fdOrder=1;
pad=2*fdOrder;
vpsmooth = mysmooth(vp, 5);
vppad = padarray(vp, [pad pad], 'replicate','both' );
n1=nz+2*pad;
n2=nx+2*pad;
srcz=pad+srcz; 
srcx=pad+srcx;
sln=1./vppad;
ttemp=tt(sln, n1, n2, h, srcz, srcx);
ttbl=mysmooth(ttemp,2);
%ttbl=ttemp;
%%
recx=srcx+nx-1-pad;% receiver x-position
recz=0:5:nz-1;
nitr=500; nrec=size(recz, 2);
path=zeros(nrec, nitr, 2);
itr=zeros([1, nrec]);
for irec=1:nrec
    [path(irec,:,:), itr(irec)]=src2rec([srcz, srcx], [recz(irec), recx], ...
                                        pad, ttbl, nitr, fdOrder);
end
% parpool(ncpu);
% parfor irec=1:nrec
%     [path(irec,:,:), itr(irec,:,:)]=src2rec([srcz, srcx], [recz, recx], pad, ttbl, nitr);
% end
%%
imagesc(h*(0:nx-1), h*(0:nz-1), vp,[3,6]);
colormap(flipud(jet)); axis tight; axis equal;
colorbar;
hold on;
% plot source and receiver
plot(h*(srcx- pad), h*(srcz-pad), 'r*');
plot(h*(recx)*ones(nrec), h*(recz), 'g*');
% plot ray path
for irec=1:nrec
    plot(h*(path(irec, 1:itr(irec)+1, 2)-pad), h*(path(irec, 1:itr(irec)+1,1)-pad),...
    'color', [0 0 0], 'LineWidth', 1);
end
hold off;
grid;
xlabel('x[m]');
ylabel('z[m]');
title('Two-point ray tracing results overlayed on the velocity model');