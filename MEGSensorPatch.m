function [X,Y,Z] = MEGSensorPatch(Rc,Ex)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a = 14e-3;
E = reshape(Ex,[3 2 61]);
R = zeros(3,61,4);
R(:,:,1) = Rc(:,1:2:122) - a*squeeze(E(:,1,:)+E(:,2,:));
R(:,:,2) = Rc(:,1:2:122) + a*squeeze(E(:,1,:)-E(:,2,:));
R(:,:,3) = Rc(:,1:2:122) + a*squeeze(E(:,1,:)+E(:,2,:));
R(:,:,4) = Rc(:,1:2:122) + a*squeeze(-E(:,1,:)+E(:,2,:));
R = permute(R,[3 2 1]);
X = R(:,:,1);
Y = R(:,:,2);
Z = R(:,:,3);
end

