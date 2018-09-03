function [M,Rq,L,snr,ev] = mneForSphereModel(Rc,N,w,center,r,C,lamda)
% Minimum Norm Estimates (MNE)法によるMEG逆問題のソルバー
% 入力引数
% Rc: センサー位置の座標．サイズは[3,num_pt,num_ch]
% N: 各列にセンサ位置の法線ベクトルを持つ行列
% w: センサ出力をシミュレートするための重み
% center: 導体球モデルの中心座標
% r : 導体球モデルの半径
% C: ノイズの分散共分散行列
% lamda: 正則化パラメータ

num_ch = size(C,1);
num_pt = length(w);

% 導体球の表面および内部上のグリッドをおよそ1cm間隔で作成
Rq = megSphereModel(center,r,1e-2);
num_dip = size(Rq,2);

% fprintf('Computation of lead field matrix\n');

% リードフィールド行列の計算
L = reshape(...
    computeMEGLeadField(Rq,reshape(Rc,[3,num_pt*num_ch]),center),...
    [3,num_pt,num_ch,3,num_dip]);
% センサ面の法線方向の磁場成分を取り出すようにする．
L = squeeze(sum(bsxfun(@times,L,reshape(N,[3,1,num_ch])),1));
% センサで計測される信号を取り出すようにする．
L = squeeze(sum(bsxfun(@times,L,reshape(w,[2,1])),1));
% ここまででLのサイズは[num_ch,3,num_dip]になる．
% これを用いて電流の3次元ベクトルから計測された信号への順問題が解ける．

% Lを整形
L = reshape(L,[num_ch,3*num_dip]);

% 各ダイポール位置における半径方向に垂直な平面の基底を得る．
Ed = dipoleAxesOnShpere(Rq,center);
% ダイポールごとの2次元座標系から元の3次元座標系への変換行列
H = zeros(3*num_dip,2*num_dip);
for i=1:num_dip
    ri = 3*(i-1)+1;
    ci = 2*(i-1)+1;
    H(ri:ri+2,ci:ci+1) = Ed(:,1:2,i);
end

% fprintf('Computation of MNE\n');
[M,snr,ev] = mne(L*H,num_dip,C,lamda,0.5);
M = H*M;
end

