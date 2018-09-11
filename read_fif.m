function [D,status] = read_fif( filen )
%READ_FIF read FIFF file(.fif). 
%   read FIFF file(.fif).

% file open
fid = fopen(filen,'r');
D = struct;
status = 0;

%% read FILE_START
t_id = zeros(1,2);
for i=1:2
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	if feof(fid)
		status = -2;
		break;
	end
	fseek(fid,tag(3),'cof');
	
	t_id(i) = tag(1);
	if tag(4)==-1
		status = -1;
		break;
	end
end

if t_id(1)~=100 || t_id(2)~=101 || tag(4)<0
	% error message
	status = -1;
end


%% read measurement data
while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	if feof(fid)
		status = -2;
		break;
	end
	
	if tag(1)==104
		b_id = fread(fid,1,'int32',0,'b');
		if b_id==100
			[D,tag,status] = read_meas(fid, tag);
		else
			[tag,status] = skip_block(fid, tag, b_id);
		end
	else
		fseek(fid,tag(3),'cof');
	end
	
	if tag(4)<0 && tag(1)==108
		break;
	elseif tag(4)<0
		status = -1;
	end
end

fclose(fid);
end



function [tag,status] = skip_block(fid,tag,b_id)
%% 
status=0;
while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	if feof(fid)
		status = -2;
		break;
	end
	
	if tag(1)==104
		b_id2 = fread(fid,1,'int32',0,'b');
		[tag,status] = skip_block(fid, tag, b_id2);
	elseif tag(1)==105
		b_id2 = fread(fid,1,'int32',0,'b');
		if b_id2~=b_id
			status = -3;
		end
		break;
	else
		fseek(fid,tag(3),'cof');
	end
end
end

function [D,tag,status] = read_meas(fid, tag)
%% 
status=0;
read_f = false;
D = struct;
while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	if tag(1)==104
		b_id = fread(fid,1,'int32',0,'b');
		if b_id==101
			[D.meas_info,tag,status] = read_meas_info(fid,tag);
 		elseif b_id==102 && ~read_f
			[D.raw_data,tag,status] = read_raw_data(fid,tag,D.meas_info.nchan);
			read_f = true;
		elseif b_id==103 && ~read_f
			[D.processed_data,tag,status] = read_processed_data(fid,tag,D.meas_info.nchan);
			read_f = true;
		else
			[tag, status] = skip_block(fid,tag,b_id);
		end
	elseif tag(1)==105
		b_id = fread(fid,1,'int32',0,'b');
		if b_id~=100
			status = -3;
		end
		break;
	else
		fseek(fid,tag(3),'cof');
	end	
end
end


function [D,tag,status] = read_meas_info(fid,tag)
%% 
D = struct;
status=0;
while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	switch tag(1)
		case 104
			b_id = fread(fid,1,'int32',0,'b');
			
			if b_id==106 % subject
				[D.subject,tag,status] = read_subject(fid,tag);
			elseif b_id==109 % hpi_result
				[D.hpi_result,tag,status] = read_hpi_result(fid,tag);
            % elseif b_id==108 % hpi_meas
            %     [D.hpi_result,tag,status] = read_hpi_meas(fid,tag);
            elseif b_id==107 % isotrak
                [D.isotrak,tag,status] = read_isotrak(fid,tag);
            elseif b_id==117
                [D.dacq_pars,tag,status] = read_dacq_pars(fid,tag);
			else
				[tag,status] = skip_block(fid, tag, b_id);
			end
		case 105
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=101
				status = -3;
			end
			break;
		case 200 % nchan
			D.nchan = fread(fid,1,'int32',0,'b');
		case 201 % sfreq
			D.sfreq = fread(fid,1,'float',0,'b');
		case 202 % data_pack
			D.data_pack = fread(fid,1,'int32',0,'b');
		case 203 % ch_info
			[D.ch_info,tag,status] = read_ch_info(fid,tag,D.nchan);
		case 204 % meas_date
			len = tag(3)/4;
			D.meas_date = fread(fid,len,'int32',0,'b')';
		case 206 % experimenter
			D.experimenter = fread(fid,tag(3),'uchar=>char',0,'b')';
		case 212 % description
			D.description = fread(fid,tag(3),'uchar=>char',0,'b')';
		case 219 % lowpass
			D.lowpass = fread(fid,1,'float',0,'b');
		case 220 % bad_chs
			len = tag(3)/4;
			D.bad_chs = fread(fid,len,'int32',0,'b')';
		case 223 % highpass
			D.highpass = fread(fid,1,'float',0,'b');
		case 500 % proj_id
			D.proj_id = fread(fid,1,'int32',0,'b');
		case 501 % proj_name
			D.proj_name = fread(fid,tag(3),'uchar=>char',0,'b')';
		otherwise
			fseek(fid,tag(3),'cof');
	end
end
end

function [ch_info,tag,status] = read_ch_info(fid,tag,nchan)
%% 
status = 0;

ch_info.scanNo = zeros(1,nchan);
ch_info.kind = zeros(1,nchan);
ch_info.range = zeros(1,nchan);
ch_info.cal = zeros(1,nchan);
ch_info.coil_type = zeros(1,nchan);
ch_info.r0 = zeros(3,nchan);
ch_info.ex = zeros(3,nchan);
ch_info.ey = zeros(3,nchan);
ch_info.ez = zeros(3,nchan);
ch_info.unit = zeros(1,nchan);
ch_info.unit_mul = zeros(1,nchan);
ch_info.val13 = zeros(1,nchan);
ch_info.val14 = zeros(1,nchan);
ch_info.val15 = zeros(1,nchan);
ch_info.val16 = zeros(1,nchan);

ch_info.scanNo(1) = fread(fid,1,'int32',4,'b');
ch_info.kind(1) = fread(fid,1,'int32',0,'b');
ch_info.range(1) = fread(fid,1,'float',0,'b');
ch_info.cal(1) = fread(fid,1,'float',0,'b');
ch_info.coil_type(1) = fread(fid,1,'int32',0,'b');
ch_info.r0(:,1) = fread(fid,3,'float',0,'b');
ch_info.ex(:,1) = fread(fid,3,'float',0,'b');
ch_info.ey(:,1) = fread(fid,3,'float',0,'b');
ch_info.ez(:,1) = fread(fid,3,'float',0,'b');
ch_info.unit(1) = fread(fid,1,'int32',0,'b');
ch_info.unit_mul(1) = fread(fid,1,'int32',0,'b');
ch_info.val13(1) = fread(fid,1,'float',0,'b');
ch_info.val14(1) = fread(fid,1,'float',0,'b');
ch_info.val15(1) = fread(fid,1,'int32',0,'b');
ch_info.val16(1) = fread(fid,1,'int32',0,'b');

for i=2:nchan
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	if tag(1)~=203
		status = -4;
		break;
	end
	
	ch_info.scanNo(i) = fread(fid,1,'int32',4,'b');
	ch_info.kind(i) = fread(fid,1,'int32',0,'b');
	ch_info.range(i) = fread(fid,1,'float',0,'b');
	ch_info.cal(i) = fread(fid,1,'float',0,'b');
	ch_info.coil_type(i) = fread(fid,1,'int32',0,'b');
	ch_info.r0(:,i) = fread(fid,3,'float',0,'b');
	ch_info.ex(:,i) = fread(fid,3,'float',0,'b');
	ch_info.ey(:,i) = fread(fid,3,'float',0,'b');
	ch_info.ez(:,i) = fread(fid,3,'float',0,'b');
	ch_info.unit(i) = fread(fid,1,'int32',0,'b');
	ch_info.unit_mul(i) = fread(fid,1,'*int32',0,'b');
	ch_info.val13(i) = fread(fid,1,'float',0,'b');
	ch_info.val14(i) = fread(fid,1,'float',0,'b');
	ch_info.val15(i) = fread(fid,1,'int32',0,'b');
	ch_info.val16(i) = fread(fid,1,'int32',0,'b');
end
end

function [S,tag,status] = read_subject(fid, tag)
%% 
S = struct;
status = 0;
% lv = false(1,3);

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
		
	switch tag(1)
		case 104
			b_id = fread(fid,1,'int32',0,'b');
			[tag,status] = skip_block(fid, tag, b_id);
		case 105
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=106
				status = -3;
			end
			break;
		case 400
			S.id = fread(fid,1,'int32',0,'b');
			% lv(1) = true;
		case 401
			S.first_name = fread(fid,tag(3),'uchar=>char',0,'b')';
			% lv(2) = true;
		case 402
			S.middle_name = fread(fid,tag(3),'uchar=>char',0,'b')';
		case 403
			S.last_name = fread(fid,tag(3),'uchar=>char',0,'b')';
			% lv(3) = true;
		case 404
			S.birth_day = fread(fid,1,'int32',0,'b');
		case 405
			S.sex = fread(fid,1,'int32',0,'b');
		case 406
			S.hand = fread(fid,1,'int32',0,'b');
		case 407
			S.weight = fread(fid,1,'float',0,'b');
		case 408
			S.height = fread(fid,1,'float',0,'b');
		case 409
			S.comment = fread(fid,tag(3),'uchar=>char',0,'b')';
		case 410
			S.his_id = fread(fid,tag(3),'uchar=>char',0,'b')';
		otherwise
			fseek(fid,tag(3),'cof');
	end
end
end

function [S,tag,status] = read_hpi_result(fid, tag)
S = struct;
status = 0;
i = 1;

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	next = tag(3) + ftell(fid);
	
	switch tag(1)
		case 104 % block_start
			b_id = fread(fid,1,'int32',0,'b');
			[tag,status] = skip_block(fid, tag, b_id);
		case 105 % block_end
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=109
				status = -3;
			end
			break
		case 213 % dig_point
			S.dig_point(i).kind = fread(fid,1,'int32',0,'b');
			S.dig_point(i).ident = fread(fid,1,'int32',0,'b');
			S.dig_point(i).r = fread(fid,3,'float',0,'b')';
			i = i + 1;
		case 222 % coord_trans
			S.coord_trans.from = fread(fid,1,'int32',0,'b');
			S.coord_trans.to = fread(fid,1,'int32',0,'b');
			S.coord_trans.rot = fread(fid,[3 3],'float',0,'b')';
			S.coord_trans.move = fread(fid,3,'float',0,'b');
			S.coord_trans.invrot = fread(fid,[3 3],'float',0,'b')';
			S.coord_trans.invmove = fread(fid,3,'float',0,'b');
		case 240 % hpi_coil_moments
			pos = ftell(fid);
			fseek(fid,tag(3)-4,'cof');
			ndim = fread(fid,1,'int32',0,'b');
			fseek(fid,-4*(ndim+1),'cof');
			sz = fread(fid,ndim,'int32',0,'b')';
			fseek(fid,pos,'bof');
			S.hpi_coil_moments = fread(fid,sz,'float',0,'b')';
			fseek(fid,4*(ndim+1),'cof');
		case 241 % hpi_fit_goodness
			len = tag(3)/4;
			S.fit_goodness = fread(fid,len,'float',0,'b')';
		case 242 % hpi_fit_accept
			S.fit_accept = fread(fid,1,'int32',0,'b');
		case 243 % hpi_fit_good_limit
			S.fit_good_limit = fread(fid,1,'float',0,'b');
		case 244 % hpi_fit_dist_limit
			S.fit_dist_limit = fread(fid,1,'float',0,'b');
		case 246 % hpi_coils_used
			len = tag(3)/4;
			S.coils_used = fread(fid,len,'int32',0,'b')';
		case 247 % hpi_digitization_order
			len = tag(3)/4;
			S.digitization_order = fread(fid,len,'int32',0,'b')';
		otherwise
			fseek(fid,tag(3),'cof');
	end
	
	if ~(next==ftell(fid) || tag(1)==104)
		fprintf('reading of data is incorrect(%d)D\n',tag(1));
		status = -1;
	end
end
end


function [S,tag,status] = read_isotrak(fid, tag)
S = struct;
status = 0;
i = 1;

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	next = tag(3) + ftell(fid);
	
	switch tag(1)
		case 105 % block_end
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=107
				status = -3;
			end
			break;
		case 213 % dig_point
			S.dig_point(i).kind = fread(fid,1,'int32',0,'b');
			S.dig_point(i).ident = fread(fid,1,'int32',0,'b');
			S.dig_point(i).r = fread(fid,3,'float',0,'b')';
			i = i + 1;
        otherwise
            fseek(fid,tag(3),'cof');
			status = -1;
            break;
	end
	
	if ~(next==ftell(fid) || tag(1)==104)
		fprintf('reading of data is incorrect(%d)D\n',tag(1));
		status = -1;
	end
end
end


function [S,tag,status] = read_dacq_pars(fid, tag)
S = struct;
status = 0;

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	next = tag(3) + ftell(fid);
	
	switch tag(1)
		case 105 % block_end
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=117
				status = -3;
			end
			break;
		case 150 % dig_point
            S.dacq_pars = fread(fid,tag(3),'uchar=>char',0,'b')';
        otherwise
            fseek(fid,tag(3),'cof');
	end
	
	if ~(next==ftell(fid) || tag(1)==104)
		fprintf('reading of data is incorrect(%d)D\n',tag(1));
		status = -1;
	end
end
end


function [raw_data,tag,status] = read_raw_data(fid, tag, nchan)
%% 
status = 0;
raw_data = int16(zeros(nchan,0));
len = 0;
samp_sz = 2*nchan;
pos = ftell(fid);

while status==0
	if tag(4)~=0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	switch tag(1)
		case 105
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=102
				status = -3;
			else
				fseek(fid,pos,'bof');
				raw_data = int16(zeros(nchan,len));
			end
			break;
		case 300 % data_buffer
			len = len + tag(3)/samp_sz;
			fseek(fid,tag(3),'cof');
% 		case 301 % data_skip
% 			
% 		case 303 % data_skip_samp
% 			
		otherwise
			fseek(fid,tag(3),'cof');
	end	
end

c_st = 1;
c_end = 0;
while status==0
	tag = fread(fid,4,'int32',0,'b');
	
	switch tag(1)
		case 105
			fseek(fid,tag(3),'cof');
			break;
		case 300 % data_buffer
			nsamp = tag(3)/samp_sz;
			c_end = c_end + nsamp;
			raw_data(:,c_st:c_end) = fread(fid,[nchan nsamp],'*int16',0,'b');
			c_st = c_end + 1;
% 		case 301 % data_skip
% 			
% 		case 303 % data_skip_samp
% 			
		otherwise
			fseek(fid,tag(3),'cof');
	end	
end

end

function [D,tag,status] = read_processed_data(fid,tag,nchan)
status = 0;
D = struct;
nevoked = 1;

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	if tag(1)==104
		b_id = fread(fid,1,'int32',0,'b');
		if b_id==104
			[D.evoked(nevoked),tag,status] = read_evoked(fid,tag,nchan);
			nevoked = nevoked + 1;
 		elseif b_id==112
			[D.continuous_data,tag,status] = read_raw_data(fid,tag,nchan);
		else
			[tag, status] = skip_block(fid,tag,b_id);
		end
	elseif tag(1)==105
		b_id = fread(fid,1,'int32',0,'b');
		if b_id~=103
			status = -3;
		end
		break;
	else
		fseek(fid,tag(3),'cof');
	end	
end
end

function [D,tag,status] = read_evoked(fid,tag,nchan)
status = 0;
D = struct;
naspect = 1;

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	switch tag(1)
		case 104
			b_id = fread(fid,1,'int32',0,'b');
			if b_id==105
				[D.aspect(naspect),tag,status] = read_aspect(fid,tag,nchan);
				naspect = naspect + 1;
			else
				[tag, status] = skip_block(fid,tag,b_id);
			end
		case 105
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=104
				status = -3;
			end
			break;
		case 206
			D.description = fread(fid,tag(3),'uchar=>char',0,'b')';
		case 208
			D.first_samp = fread(fid,1,'int32',0,'b');
		case 209
			D.last_samp = fread(fid,1,'int32',0,'b');
		case 211
			D.ref_event = fread(fid,1,'int32',0,'b');
		case 227
			D.event_comment = fread(fid,tag(3),'uchar=>char',0,'b')';
		otherwise
			fseek(fid,tag(3),'cof');
	end	
end
end

function [D,tag,status] = read_aspect(fid,tag,nchan)
status = 0;
D = struct;
flag = false;
D.epoch = zeros(nchan,0);
i = 1;

while status==0
	if tag(4)>0
		fseek(fid,tag(4),'bof');
	elseif tag(4)<0
		status = -1;
		break;
	end
	
	if feof(fid)
		status = -2;
		break;
	end
	
	tag = fread(fid,4,'int32',0,'b');
	
	if feof(fid)
		status = -2;
		break;
	end
	
	switch tag(1)
		case 105
			b_id = fread(fid,1,'int32',0,'b');
			if b_id~=105
				status = -3;
			end
			break;
		case 210
			D.aspect_kind = fread(fid,1,'int32',0,'b');
		case 207
			D.nave = fread(fid,1,'int32',0,'b');
		case 302
			if ~flag
				nsamp = (tag(3)-8)/2;
				D.epoch = zeros(nchan,nsamp);
				% buffer = int16(zeros(nsamp,1));
				flag = true;
			end
			temp = fread(fid,2,'float=>double',0,'b');
			buffer = fread(fid,nsamp,'int16=>double',0,'b');
			D.epoch(i,:) = temp(2)*buffer' + temp(1);
			i = i + 1;
		otherwise
			fseek(fid,tag(3),'cof');
	end	
end
end
