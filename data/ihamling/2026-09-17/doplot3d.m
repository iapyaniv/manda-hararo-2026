%%% script to plot fault geometry in 3D?

% saveview = [   -0.7370    3.0025    0.0000   -1.1327
%    -0.2860   -1.3856    2.9849   -0.6566
%    -0.6124   -2.9666   -1.3942   30.5886
%          0         0         0    1.0000];

%view(saveview);
plot3d = load('plot_3d_dyke_only.dat');
faults = plot3d;
%faults(faults(:,4)<1,:) = [];
op = load('dyke_only_apr08_ex.txt');
%caxis([0 3]);
d2r=pi/180;
hold off
%caxis([0 3]);
[rp,cp] = size(faults);

%tmp = (faults(:,9)==2);
%tmp = tmp + (faults(:,4)==90);
%for i = [rp:-1:1]
%    if (tmp(i) == 2)
%        faults(i,:)=[];
%    end
%end
[rp,cp] = size(faults);


for p = 1:rp
xcent = faults(p,1);
ycent = faults(p,2);
length = faults(p,7)*1000;
strike = faults(p,3)*d2r;
top = faults(p,8)*1000;
bot = faults(p,9)*1000;
slip = op(p,1);
dip = faults(p,4)*d2r;
%dip = 45*d2r;

x_surf1 = xcent - 0.5*length*sin(strike);
x_surf2 = xcent + 0.5*length*sin(strike);
y_surf1 = ycent - 0.5*length*cos(strike);
y_surf2 = ycent + 0.5*length*cos(strike);

x_top1 = x_surf1 + top*cos(strike)/tan(dip);
x_top2 = x_surf2 + top*cos(strike)/tan(dip);
y_top1 = y_surf1 - top*sin(strike)/tan(dip);
y_top2 = y_surf2 - top*sin(strike)/tan(dip);


x_bot1 = x_surf1 + bot*cos(strike)/tan(dip);
x_bot2 = x_surf2 + bot*cos(strike)/tan(dip);
y_bot1 = y_surf1 - bot*sin(strike)/tan(dip);
y_bot2 = y_surf2 - bot*sin(strike)/tan(dip);

X(1:4,p) = [x_top1 x_top2 x_bot2 x_bot1]'/1000;
Y(1:4,p) = [y_top1 y_top2 y_bot2 y_bot1]'/1000;
Z(1:4,p) = [-top  -top  -bot  -bot]'/1000;

%subplot(4,1,1:3)
fill3(X(:,p),Y(:,p),Z(:,p),slip)
%caxis([0 .05])
%view(saveview);
axis equal
hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%caxis([0 max(plot3d(:,6))])
colorbar('vert');
make_slipcol
colormap(col);
caxis([0 3])

%% add mogi sphere

%1. make unit sphere
[xs ys zs] = sphere(10);

xs1=xs*1.6+ 677.6;
ys1=ys*1.6+ 1360;
zs1=zs*1.6-8;

%xs2=xs*2+660.7; 
%ys2=ys*2+1394.1 ;
%zs2=zs*2-5;

surf(xs1,ys1,zs1);
%surf(xs2,ys2,zs2);

%load azel.mat
%view(az,el)
