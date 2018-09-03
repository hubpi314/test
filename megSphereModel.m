function Rg = megSphereModel(c,r,d)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

th_max = 90*pi/180;

r_sub = r:-d:d/2;
div_th = ceil(th_max*r_sub/d);

num_point = 0;
for i=1:length(r_sub);
    theta = (1:div_th(i))*th_max/div_th(i);
    div_ph = ceil(r_sub(i)*sin(theta)*2*pi/d);
    num_point = num_point + sum(div_ph) + 1;
end

Rg = zeros(3,sum(num_point));

st = 1;
for i=1:length(r_sub);
    theta = (1:div_th(i))*th_max/div_th(i);
    div_ph = ceil(r_sub(i)*sin(theta)*2*pi/d);
    
    for j=1:length(div_ph)
        phi = (0:div_ph(j)-1)*2*pi/div_ph(j);
        
        en = st + div_ph(j) - 1;
        Rg(1,st:en) = r_sub(i)*sin(theta(j))*cos(phi);
        Rg(2,st:en) = r_sub(i)*sin(theta(j))*sin(phi);
        Rg(3,st:en) = r_sub(i)*cos(theta(j));
        st = en + 1;
    end
    
    Rg(3,st) = r_sub(i);
    st = st + 1;
end

Rg = bsxfun(@plus,Rg,c);

end

