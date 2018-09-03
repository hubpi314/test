function [M,Rq,L,snr,ev] = mneForSphereModel(Rc,N,w,center,r,C,lamda)
% Minimum Norm Estimates (MNE)�@�ɂ��MEG�t���̃\���o�[
% ���͈���
% Rc: �Z���T�[�ʒu�̍��W�D�T�C�Y��[3,num_pt,num_ch]
% N: �e��ɃZ���T�ʒu�̖@���x�N�g�������s��
% w: �Z���T�o�͂��V�~�����[�g���邽�߂̏d��
% center: ���̋����f���̒��S���W
% r : ���̋����f���̔��a
% C: �m�C�Y�̕��U�����U�s��
% lamda: �������p�����[�^

num_ch = size(C,1);
num_pt = length(w);

% ���̋��̕\�ʂ���ѓ�����̃O���b�h�����悻1cm�Ԋu�ō쐬
Rq = megSphereModel(center,r,1e-2);
num_dip = size(Rq,2);

% fprintf('Computation of lead field matrix\n');

% ���[�h�t�B�[���h�s��̌v�Z
L = reshape(...
    computeMEGLeadField(Rq,reshape(Rc,[3,num_pt*num_ch]),center),...
    [3,num_pt,num_ch,3,num_dip]);
% �Z���T�ʂ̖@�������̎��ꐬ�������o���悤�ɂ���D
L = squeeze(sum(bsxfun(@times,L,reshape(N,[3,1,num_ch])),1));
% �Z���T�Ōv�������M�������o���悤�ɂ���D
L = squeeze(sum(bsxfun(@times,L,reshape(w,[2,1])),1));
% �����܂ł�L�̃T�C�Y��[num_ch,3,num_dip]�ɂȂ�D
% �����p���ēd����3�����x�N�g������v�����ꂽ�M���ւ̏���肪������D

% L�𐮌`
L = reshape(L,[num_ch,3*num_dip]);

% �e�_�C�|�[���ʒu�ɂ����锼�a�����ɐ����ȕ��ʂ̊��𓾂�D
Ed = dipoleAxesOnShpere(Rq,center);
% �_�C�|�[�����Ƃ�2�������W�n���猳��3�������W�n�ւ̕ϊ��s��
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

