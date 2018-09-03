function Ed = dipoleAxesOnShpere(X,c)
% ���S�̈ʒu�x�N�g����c�̋��̕\�ʂɐڂ�����W���W�n�̒P�ʃx�N�g����Ԃ��D
% X: �e��Ɋe�_��3�������W�����ׂ�ꂽ�s��
% c: ���̒��S
num_dip = size(X,2);
Rd = bsxfun(@minus,X,c);

Ed = zeros(3,3,num_dip);
Ed(1:2,1,:) = bsxfun(@times,Rd(1:2,:),Rd(3,:));
Ed(1,2,:) = -Rd(2,:);
Ed(2,2,:) = Rd(1,:);
tmp = sqrt(sum(Rd(1:2,:).^2,1)); % xy���ʐ����̑傫��
lv = tmp<1e-6;
Ed(3,1,:) = -tmp;
Ed(1:2,1,:) = bsxfun(@rdivide,Ed(1:2,1,:),reshape(tmp,[1 1 num_dip]));
Ed(1:2,2,:) = bsxfun(@rdivide,Ed(1:2,2,:),reshape(tmp,[1 1 num_dip]));
tmp = 1./sqrt(sum(Rd.^2,1));
Ed(:,1,:) = bsxfun(@times,Ed(:,1,:),reshape(tmp,[1 1 num_dip]));

Ed(:,3,:) = bsxfun(@times,Rd,tmp); % �@���x�N�g��

% �_�C�|�[���ʒu��z����i���ٓ_�j�̏ꍇ
tmp = sum(lv);
if tmp>0
    Ed(:,1,lv) = repmat([1;0;0],[1 1 tmp]);
    Ed(:,2,lv) = repmat([0;1;0],[1 1 tmp]);
    Ed(:,3,lv) = repmat([0;0;1],[1 1 tmp]);
end
end