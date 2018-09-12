% Minimum Norm Estimates (MNE)のテストスクリプト

close all;
clear variables;

% addpath('../meglib','../fifffile');
filename = '../hoge';
ind_stim = 1;
ind_t = 118;
badch = [28,99,100,103,111,120];
num_ch = 122;

ind_ch = 1:num_ch;
% ind_ch = [...
% %          14, 20, 18,...
% %          13, 19, 17,...
%          12, 10,  8,...
%          11,  9,  7,...
%           6,  4,  2,...
%           5,  3,  1,...
%         116,114,112,...
%         115,113,111
%         ];


% Fiffファイルの読み込み
[S, status] = read_fif(filename);
Rch = S.meas_info.ch_info.r0(:,1:num_ch); % センサの位置ベクトル
Ex = S.meas_info.ch_info.ex(:,1:num_ch); % 
Ey = S.meas_info.ch_info.ey(:,1:num_ch);
N = S.meas_info.ch_info.ez(:,1:num_ch); % 各センサの法線ベクトル

% HPIコイルの位置ベクトル
Rhpi = zeros(3,3);
Rhpi(:,1)=S.meas_info.hpi_result.dig_point(1).r;
Rhpi(:,2)=S.meas_info.hpi_result.dig_point(2).r;
Rhpi(:,3)=S.meas_info.hpi_result.dig_point(3).r;

% Device座標系からHead座標系に変換する回転行列と平行移動
Rot = S.meas_info.hpi_result.coord_trans.rot;
mv = S.meas_info.hpi_result.coord_trans.move;


% 加算波形の読み込み
X = S.processed_data.evoked(ind_stim).aspect(1).epoch(1:122,:);
cal = S.meas_info.ch_info.cal(1:122)*1e13;
X = bsxfun(@times,X,cal');
ind_base =  -S.processed_data.evoked(ind_stim).first_samp;


% filter setting
sf = S.meas_info.sfreq;
wc = [2 30]*2/sf; % 2-30
wt = [2 2]*2/sf; % 2,20
fftsize = 1024;

clear S;


% センサ位置とその座標系を取得
ind_odd = 1:2:num_ch;% 奇数チャネルのインデックス
Rs = Rch(:,ind_odd);
Ex_s = Ex(:,ind_odd);
Ey_s = Ey(:,ind_odd);
N_s = N(:,ind_odd);


% チャネルの状態を示す論理ベクトル
sta_ch = false(num_ch,1);
sta_ch(ind_ch) = true;
sta_ch(badch) = false;
num_ch = nnz(sta_ch);

% センサを示すパッチ
% [Xs, Ys, Zs] = MEGSensorPatch(Rch,Ex);

% 使用するチャネルの限定
Rch = Rch(:,sta_ch);
Ex = Ex(:,sta_ch);
N = N(:,sta_ch);
X = X(sta_ch,:);


% フィルタリング
Y = neuromagfilter2(X,wc,wt,'bandpass',fftsize,2);
% ベースライン処理
Y = bsxfun(@minus,Y,mean(Y(:,1:ind_base),2));

% ノイズの分散共分散行列の推定
C = Y(:,1:ind_base)*Y(:,1:ind_base)'/(ind_base-1);
% 分散共分散行列の正則化
ep = 0.05;
C = C + ep*trace(C)*eye(num_ch)/num_ch;

% Neuromag122の固定パラメータ
dx = 8.1e-3; % センサの中心からコイルの中心までの距離は8.1mm
w = 0.5./[dx -dx]; % センサにおける空間微分を近似するための重み


% コイルの中心の位置ベクトル，サイズは[3, 2, num_ch]
Rc = reshape(bsxfun(@plus,Rch(:),Ex(:)*[dx -dx]),[3,num_ch*2]);


% Head座標系への変換
Rc = permute(reshape(bsxfun(@plus,Rot*Rc,mv),...
    [3,num_ch,2]),[1,3,2]);
N = Rot*N;
Rhpi = bsxfun(@plus,Rot*Rhpi,mv);

% 導体球モデルのパラメータ
center = [0;0;0.04];% Head座標系における中心
% 導体球モデルの中心からHPIコイルまでの距離
dist = sqrt(sum(bsxfun(@minus,Rhpi,center).^2,1));


% MNEの正則化パラメータ
lamda = 7;


% Head座標系でMNEの解を計算．
[M,Rq,L,snr,ev] = mneForSphereModel(Rc,N,w,center,0.9*mean(dist),C,lamda);
q = M*Y(:,ind_t);

num_dip = size(Rq,2);
fprintf('Number of grids: %d\n',num_dip);

% 誤差の計算
y_p = L*q;
err = Y(:,ind_t) - y_p;
var_err = err'*C^-1*err;
fprintf('Power SNR: %f\n',snr);
fprintf('Loss: %f\n',var_err);
figure;
plot(ev);

% Device座標系に変換
Q = Rot'*reshape(q,[3,size(Rq,2)]);


% Device座標系に座標変換
Rq = Rot'*bsxfun(@minus,Rq,mv);
center = Rot'*(center-mv);


% 等磁界線図を計算
num_div = 64; % 2次元グリッドの分割数
normalMagneticContour(Q,Rq,center,Rs,Ex_s,Ey_s,N_s,num_div);


% 電流分布のプロット
figure;
quiver3(Rq(1,:),Rq(2,:),Rq(3,:),Q(1,:),Q(2,:),Q(3,:));
line('XData',Rs(1,:),'YData',Rs(2,:),'ZData',Rs(3,:),...
    'LineStyle','none','Marker','o');
% patch(Xs,Ys,Zs,'green','FaceAlpha',0);
line('XData',center(1),'YData',center(2,:),'ZData',center(3,:),...
    'LineStyle','none','Marker','o','Color','r');
