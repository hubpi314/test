function L = computeMEGLeadField(Rq,Rc,c)
% ダイポール磁場のリードフィールド行列を計算する．
% 入力引数
% Rq: ダイポールの位置ベクトルを各列に持つ行列．サイズは[3 num_dip]
% Rc: 観測点の位置ベクトルで，サイズは[3, num_ch]
% 出力引数
% L: リードフィールド行列で，サイズは[3,num_ch,3,num_dip]．
% squeeze(L(:,ch,:,q))*Jで，モーメントJのダイポールqがチャネルchに作り出す
% 3次元の磁場ベクトルを計算することができる．

num_dip = size(Rq,2);
num_ch = size(Rc,2);

Diff = bsxfun(@minus,Rc,reshape(Rq,[3,1,1,num_dip]));
Diff = bsxfun(@times,Diff,sum(Diff.*Diff,1).^-1.5);
L = 1e-3*cross(repmat(reshape(eye(3),[3,1,3]),[1,num_ch,1,num_dip]),...
    repmat(Diff,[1,1,3,1]));

end