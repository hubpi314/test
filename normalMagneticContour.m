function h = normalMagneticContour(Q,Rq,center,Rs,Ex,Ey,N,num_div)
%UNTITLED 等磁場線図の描画
%   Detailed explanation goes here
% 

num_grid_div = 5;
num_ch = size(Rs,2);
num_grid = num_ch*num_grid_div^2;
num_dip = size(Rq,2);

Rg = reshape(sensorSurfaceGrid(Rs,Ex,Ey,num_grid_div),[3,num_grid]);

% 電流源からセンサ表面へのリードフィールド行列
L = computeMEGLeadField(Rq,Rg,center);
b = reshape(L,[3*num_grid,3*num_dip])*Q(:);
% センサ面に対する磁場の法線成分を計算．
b = sum(bsxfun(@times,N,reshape(b,[3,num_ch,num_grid_div^2])),1);
b = reshape(b,[num_grid,1]);

% チャネル位置を2次元平面にマッピング．
dist = sqrt(sum(Rs.*Rs,1));
theta = acos(Rs(3,:)./dist);
Rs_2d = bsxfun(@times,Rs,theta./(dist.*sin(theta)));

% グリッドを2次元平面にマッピング．
dist = sqrt(sum(Rg.*Rg,1));
theta = acos(Rg(3,:)./dist);
Rg_2d = bsxfun(@times,Rg,theta./(dist.*sin(theta)));

lwr = min(Rg_2d,[],2);
upr = max(Rg_2d,[],2);

ind_x = linspace(lwr(1),upr(1),num_div);
ind_y = linspace(lwr(2),upr(2),num_div);
[X,Y] = meshgrid(ind_x,ind_y);

b2 = griddata(Rg_2d(1,:)',Rg_2d(2,:)',b,X(:),Y(:));

% 等磁場線図の描画
h = figure;
contourf(ind_x,ind_y,reshape(b2,[num_div,num_div]));
% チャネル位置の表示
line('XData',Rs_2d(1,:),'YData',Rs_2d(2,:),'LineStyle','none','Marker','o');
axis equal;
end

