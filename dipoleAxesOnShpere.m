function Ed = dipoleAxesOnShpere(X,c)
% 中心の位置ベクトルがcの球の表面に接する座標座標系の単位ベクトルを返す．
% X: 各列に各点の3次元座標が並べられた行列
% c: 球の中心
num_dip = size(X,2);
Rd = bsxfun(@minus,X,c);

Ed = zeros(3,3,num_dip);
Ed(1:2,1,:) = bsxfun(@times,Rd(1:2,:),Rd(3,:));
Ed(1,2,:) = -Rd(2,:);
Ed(2,2,:) = Rd(1,:);
tmp = sqrt(sum(Rd(1:2,:).^2,1)); % xy平面成分の大きさ
lv = tmp<1e-6;
Ed(3,1,:) = -tmp;
Ed(1:2,1,:) = bsxfun(@rdivide,Ed(1:2,1,:),reshape(tmp,[1 1 num_dip]));
Ed(1:2,2,:) = bsxfun(@rdivide,Ed(1:2,2,:),reshape(tmp,[1 1 num_dip]));
tmp = 1./sqrt(sum(Rd.^2,1));
Ed(:,1,:) = bsxfun(@times,Ed(:,1,:),reshape(tmp,[1 1 num_dip]));

Ed(:,3,:) = bsxfun(@times,Rd,tmp); % 法線ベクトル

% ダイポール位置がz軸上（特異点）の場合
tmp = sum(lv);
if tmp>0
    Ed(:,1,lv) = repmat([1;0;0],[1 1 tmp]);
    Ed(:,2,lv) = repmat([0;1;0],[1 1 tmp]);
    Ed(:,3,lv) = repmat([0;0;1],[1 1 tmp]);
end
end