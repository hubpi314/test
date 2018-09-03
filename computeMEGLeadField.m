function L = computeMEGLeadField(Rq,Rc,c)
% �_�C�|�[������̃��[�h�t�B�[���h�s����v�Z����D
% ���͈���
% Rq: �_�C�|�[���̈ʒu�x�N�g�����e��Ɏ��s��D�T�C�Y��[3 num_dip]
% Rc: �ϑ��_�̈ʒu�x�N�g���ŁC�T�C�Y��[3, num_ch]
% �o�͈���
% L: ���[�h�t�B�[���h�s��ŁC�T�C�Y��[3,num_ch,3,num_dip]�D
% squeeze(L(:,ch,:,q))*J�ŁC���[�����gJ�̃_�C�|�[��q���`���l��ch�ɍ��o��
% 3�����̎���x�N�g�����v�Z���邱�Ƃ��ł���D

num_dip = size(Rq,2);
num_ch = size(Rc,2);

Diff = bsxfun(@minus,Rc,reshape(Rq,[3,1,1,num_dip]));
Diff = bsxfun(@times,Diff,sum(Diff.*Diff,1).^-1.5);
L = 1e-3*cross(repmat(reshape(eye(3),[3,1,3]),[1,num_ch,1,num_dip]),...
    repmat(Diff,[1,1,3,1]));

end