% Minimum Norm Estimates (MNE)�̃e�X�g�X�N���v�g

close all;
clear variables;

addpath('../meglib','../fifffile');
filename = '../yano_aircon_sound/shukunami_atsuto/141224/cool02081632.fif';
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


% Fiff�t�@�C���̓ǂݍ���
[S, status] = read_fif(filename);
Rch = S.meas_info.ch_info.r0(:,1:num_ch); % �Z���T�̈ʒu�x�N�g��
Ex = S.meas_info.ch_info.ex(:,1:num_ch); % 
Ey = S.meas_info.ch_info.ey(:,1:num_ch);
N = S.meas_info.ch_info.ez(:,1:num_ch); % �e�Z���T�̖@���x�N�g��

% HPI�R�C���̈ʒu�x�N�g��
Rhpi = zeros(3,3);
Rhpi(:,1)=S.meas_info.hpi_result.dig_point(1).r;
Rhpi(:,2)=S.meas_info.hpi_result.dig_point(2).r;
Rhpi(:,3)=S.meas_info.hpi_result.dig_point(3).r;

% Device���W�n����Head���W�n�ɕϊ������]�s��ƕ��s�ړ�
Rot = S.meas_info.hpi_result.coord_trans.rot;
mv = S.meas_info.hpi_result.coord_trans.move;


% ���Z�g�`�̓ǂݍ���
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


% �Z���T�ʒu�Ƃ��̍��W�n���擾
ind_odd = 1:2:num_ch;% ��`���l���̃C���f�b�N�X
Rs = Rch(:,ind_odd);
Ex_s = Ex(:,ind_odd);
Ey_s = Ey(:,ind_odd);
N_s = N(:,ind_odd);


% �`���l���̏�Ԃ������_���x�N�g��
sta_ch = false(num_ch,1);
sta_ch(ind_ch) = true;
sta_ch(badch) = false;
num_ch = nnz(sta_ch);

% �Z���T�������p�b�`
% [Xs, Ys, Zs] = MEGSensorPatch(Rch,Ex);

% �g�p����`���l���̌���
Rch = Rch(:,sta_ch);
Ex = Ex(:,sta_ch);
N = N(:,sta_ch);
X = X(sta_ch,:);


% �t�B���^�����O
Y = neuromagfilter2(X,wc,wt,'bandpass',fftsize,2);
% �x�[�X���C������
Y = bsxfun(@minus,Y,mean(Y(:,1:ind_base),2));

% �m�C�Y�̕��U�����U�s��̐���
C = Y(:,1:ind_base)*Y(:,1:ind_base)'/(ind_base-1);
% ���U�����U�s��̐�����
ep = 0.05;
C = C + ep*trace(C)*eye(num_ch)/num_ch;

% Neuromag122�̌Œ�p�����[�^
dx = 8.1e-3; % �Z���T�̒��S����R�C���̒��S�܂ł̋�����8.1mm
w = 0.5./[dx -dx]; % �Z���T�ɂ������Ԕ������ߎ����邽�߂̏d��


% �R�C���̒��S�̈ʒu�x�N�g���C�T�C�Y��[3, 2, num_ch]
Rc = reshape(bsxfun(@plus,Rch(:),Ex(:)*[dx -dx]),[3,num_ch*2]);


% Head���W�n�ւ̕ϊ�
Rc = permute(reshape(bsxfun(@plus,Rot*Rc,mv),...
    [3,num_ch,2]),[1,3,2]);
N = Rot*N;
Rhpi = bsxfun(@plus,Rot*Rhpi,mv);

% ���̋����f���̃p�����[�^
center = [0;0;0.04];% Head���W�n�ɂ����钆�S
% ���̋����f���̒��S����HPI�R�C���܂ł̋���
dist = sqrt(sum(bsxfun(@minus,Rhpi,center).^2,1));


% MNE�̐������p�����[�^
lamda = 7;


% Head���W�n��MNE�̉����v�Z�D
[M,Rq,L,snr,ev] = mneForSphereModel(Rc,N,w,center,0.9*mean(dist),C,lamda);
q = M*Y(:,ind_t);

num_dip = size(Rq,2);
fprintf('Number of grids: %d\n',num_dip);

% �덷�̌v�Z
y_p = L*q;
err = Y(:,ind_t) - y_p;
var_err = err'*C^-1*err;
fprintf('Power SNR: %f\n',snr);
fprintf('Loss: %f\n',var_err);
figure;
plot(ev);

% Device���W�n�ɕϊ�
Q = Rot'*reshape(q,[3,size(Rq,2)]);


% Device���W�n�ɍ��W�ϊ�
Rq = Rot'*bsxfun(@minus,Rq,mv);
center = Rot'*(center-mv);


% �����E���}���v�Z
num_div = 64; % 2�����O���b�h�̕�����
normalMagneticContour(Q,Rq,center,Rs,Ex_s,Ey_s,N_s,num_div);


% �d�����z�̃v���b�g
figure;
quiver3(Rq(1,:),Rq(2,:),Rq(3,:),Q(1,:),Q(2,:),Q(3,:));
line('XData',Rs(1,:),'YData',Rs(2,:),'ZData',Rs(3,:),...
    'LineStyle','none','Marker','o');
% patch(Xs,Ys,Zs,'green','FaceAlpha',0);
line('XData',center(1),'YData',center(2,:),'ZData',center(3,:),...
    'LineStyle','none','Marker','o','Color','r');
