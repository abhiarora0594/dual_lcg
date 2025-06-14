clear all;

ngp = 25; % no of gps

% scal1 = sqrt((3/7) - (2/7)*sqrt(6/5));
% scal2 = sqrt((3/7) + (2/7)*sqrt(6/5));
% 
% scal3 = (18 + sqrt(30))/36;
% scal4 = (18 - sqrt(30))/36;

scal1 = 0;
scal2 = (1/3)*sqrt(5 - 2*sqrt(10/7));
scal3 = 1/3*sqrt(5 + 2*sqrt(10/7));

scal4 = 128/225;
scal5 = (322+13*sqrt(70))/900;
scal6 = (322-13*sqrt(70))/900;

coords = zeros(4,3); % x, y & z coordinates for an element

% coords(1,:) = [0.9692931180000000, -0.0446239933000000, 0.0000000000000000]; 
% coords(2,:) = [0.9720588330000000, 0.0000000000000000, 0.0000000000000000];
% coords(3,:) = [0.9441176650000001, 0.0000000000000000, 0.0000000000000000];
% coords(4,:) = [0.9364338520000000, -0.0524140708000000, 0.0000000000000000];
% 
% coords(1,:) = [0.9991992710000001 -0.0400096066000000 0.0000000000000000]; 
% coords(2,:) = [1.0000000000000000 0.0000000000000000 0.0000000000000000];
% coords(3,:) = [0.9720588330000000, 0.0000000000000000, 0.0000000000000000];
% coords(4,:) = [0.9692931180000000, -0.0446239933000000, 0.0000000000000000];

% coords(1,:) = [0.9720588330000000, 0.0000000000000000, 0.0000000000000000]; 
% coords(2,:) = [1.0000000000000000 0.0000000000000000 0.0000000000000000];
% coords(3,:) = [0.9991992710000001 0.0400096066000000 0.0000000000000000];
% coords(4,:) = [0.9715160130000000 0.0393345691000000 0.0000000000000000];

% coords(1,:) = [0.9442617890000000 0.0397517942000000 0.0000000000000000]; 
% coords(2,:) = [0.9441176650000001, 0.0000000000000000, 0.0000000000000000];
% coords(3,:) = [0.9720588330000000, 0.0000000000000000, 0.0000000000000000];
% coords(4,:) = [0.9715160130000000 0.0393345691000000 0.0000000000000000];


coords(1,:) = [0.9997997880000000 -0.0200088099000000 0.0000000000000000]; 
coords(2,:) = [1.0000000000000000 0.0000000000000000 0.0000000000000000];
coords(3,:) = [0.9858208890000000 0.0000000000000000 0.0000000000000000];
coords(4,:) = [0.9857366680000000 -0.0201696847000000 0.0000000000000000];

coords(1,:) = [0.9858208890000000 0.0000000000000000 0.0000000000000000]; 
coords(2,:) = [0.9716417790000000 0.0000000000000000 0.0000000000000000];
coords(3,:) = [0.9717174170000000 -0.0207224451000000 0.0000000000000000];
coords(4,:) = [0.9857366680000000 -0.0201696847000000 0.0000000000000000];
% % 
coords(1,:) = [0.9845005270000000 0.0222788677000000 0.0000000000000000]; 
coords(2,:) = [0.9696943760000000 0.0240827166000000 0.0000000000000000];
coords(3,:) = [0.9716417790000000 0.0000000000000000 0.0000000000000000];
coords(4,:) = [0.9858208890000000 0.0000000000000000 0.0000000000000000];
% 
coords(1,:) = [0.9858208890000000 0.0000000000000000 0.0000000000000000]; 
coords(2,:) = [1.0000000000000000 0.0000000000000000 0.0000000000000000];
coords(3,:) = [0.9997997880000000 0.0200088099000000 0.0000000000000000];
coords(4,:) = [0.9845005270000000 0.0222788677000000 0.0000000000000000];


xi = zeros(ngp,2); % gps coordinates in master element 
wt = zeros(ngp,1); % wt of gps

if (ngp == 4)
    xi(1,1) = -0.577350269189626; xi(1,2) = -0.577350269189626;
    xi(2,1) = 0.577350269189626; xi(2,2) = -0.577350269189626;
    xi(3,1) = 0.577350269189626; xi(3,2) = 0.577350269189626;
    xi(4,1) = -0.577350269189626; xi(4,2) = 0.577350269189626;

    wt(1,1) = 1.0; wt(2,1) = 1.0; wt(3,1) = 1.0; wt(4,1) = 1.0;

elseif (ngp == 16)
    
    xi(1,1) = -scal2; xi(1,2) = -scal2;
    xi(2,1) = -scal1; xi(2,2) = -scal2;
    xi(3,1) = scal1; xi(3,2) = -scal2;
    xi(4,1) = scal2; xi(4,2) = -scal2;

    xi(5,1) = -scal2; xi(5,2) = -scal1;
    xi(6,1) = -scal1; xi(6,2) = -scal1;
    xi(7,1) = scal1; xi(7,2) = -scal1;
    xi(8,1) = scal2; xi(8,2) = -scal1;

    xi(9,1) = -scal2; xi(9,2) = scal1;
    xi(10,1) = -scal1; xi(10,2) = scal1;
    xi(11,1) = scal1; xi(11,2) = scal1;
    xi(12,1) = scal2; xi(12,2) = scal1;

    xi(13,1) = -scal2; xi(13,2) = scal2;
    xi(14,1) = -scal1; xi(14,2) = scal2;
    xi(15,1) = scal1; xi(15,2) = scal2;
    xi(16,1) = scal2; xi(16,2) = scal2;
    
    wt(1,1) = scal4*scal4;
    wt(2,1) = scal3*scal4;
    wt(3,1) = scal3*scal4;
    wt(4,1) = scal4*scal4;

    wt(5,1) = scal4*scal3;
    wt(6,1) = scal3*scal3;
    wt(7,1) = scal3*scal3;
    wt(8,1) = scal4*scal3;

    wt(9,1) = scal4*scal3;
    wt(10,1) = scal3*scal3;
    wt(11,1) = scal3*scal3;
    wt(12,1) = scal4*scal3;

    wt(13,1) = scal4*scal4;
    wt(14,1) = scal3*scal4;
    wt(15,1) = scal3*scal4;
    wt(16,1) = scal4*scal4;

elseif (ngp == 25)
    
    xi(1,1) = -scal3; xi(1,2) = -scal3; wt(1,1) = scal6*scal6;
    xi(2,1) = -scal2; xi(2,2) = -scal3; wt(2,1) = scal5*scal6;
    xi(3,1) = scal1; xi(3,2) = -scal3;  wt(3,1) = scal4*scal6;
    xi(4,1) = scal2; xi(4,2) = -scal3;  wt(4,1) = scal5*scal6;
    xi(5,1) = scal3; xi(5,2) = -scal3;  wt(5,1) = scal6*scal6;

    xi(6,1) = -scal3; xi(6,2) = -scal2; wt(6,1) = scal6*scal5;
    xi(7,1) = -scal2; xi(7,2) = -scal2; wt(7,1) = scal5*scal5;
    xi(8,1) = scal1; xi(8,2) = -scal2;  wt(8,1) = scal4*scal5;
    xi(9,1) = scal2; xi(9,2) = -scal2; wt(9,1) = scal5*scal5;
    xi(10,1) = scal3; xi(10,2) = -scal2; wt(10,1) = scal6*scal5;

    xi(11,1) = -scal3; xi(11,2) = scal1; wt(11,1) = scal6*scal4;
    xi(12,1) = -scal2; xi(12,2) = scal1; wt(12,1) = scal5*scal4;
    xi(13,1) = scal1; xi(13,2) = scal1; wt(13,1) = scal4*scal4;
    xi(14,1) = scal2; xi(14,2) = scal1; wt(14,1) = scal5*scal4;
    xi(15,1) = scal3; xi(15,2) = scal1; wt(15,1) = scal6*scal4;

    xi(16,1) = -scal3; xi(16,2) = scal2; wt(16,1) = scal6*scal5;
    xi(17,1) = -scal2; xi(17,2) = scal2; wt(17,1) = scal5*scal5;
    xi(18,1) = scal1; xi(18,2) = scal2; wt(18,1) = scal4*scal5;
    xi(19,1) = scal2; xi(19,2) = scal2; wt(19,1) = scal5*scal5;
    xi(20,1) = scal3; xi(20,2) = scal2; wt(20,1) = scal6*scal5; 

    xi(21,1) = -scal3; xi(21,2) = scal3; wt(21,1) = scal6*scal6;
    xi(22,1) = -scal2; xi(22,2) = scal3; wt(22,1) = scal5*scal6;
    xi(23,1) = scal1; xi(23,2) = scal3; wt(23,1) = scal4*scal6;
    xi(24,1) = scal2; xi(24,2) = scal3; wt(24,1) = scal5*scal6;
    xi(25,1) = scal3; xi(25,2) = scal3; wt(25,1) = scal6*scal6;

end




psi = zeros(ngp,4);
dpsi = zeros(ngp,4,2);
  
% shape function values at different gps
for i=1:ngp
    psi(i,1) =  0.25*(1.0-xi(i,1))*(1.0-xi(i,2));
    psi(i,2) =  0.25*(1.0+xi(i,1))*(1.0-xi(i,2));
    psi(i,3) =  0.25*(1.0+xi(i,1))*(1.0+xi(i,2));
    psi(i,4) =  0.25*(1.0-xi(i,1))*(1.0+xi(i,2));
end

% shape function derivative values at different gps
for i=1:ngp

    dpsi(i,1,1) =  -0.25*(1.0-xi(i,2));
    dpsi(i,2,1) =  0.25*(1.0-xi(i,2));
    dpsi(i,3,1) =  0.25*(1.0+xi(i,2));
    dpsi(i,4,1) =  -0.25*(1.0+xi(i,2));

    dpsi(i,1,2) =  -0.25*(1.0-xi(i,1));
    dpsi(i,2,2) =  -0.25*(1.0+xi(i,1));
    dpsi(i,3,2) =  0.25*(1.0+xi(i,1));
    dpsi(i,4,2) =  0.25*(1.0-xi(i,1));

end

% real basis vectors
E_real= zeros(ngp,6);
for i=1:ngp    
    for alpha = 1:2
        for node_id =1:4
            E_real(i,3*(alpha-1)+1) = E_real(i,3*(alpha-1)+1) + coords(node_id,1)*dpsi(i,node_id,alpha); 
            E_real(i,3*(alpha-1)+2) = E_real(i,3*(alpha-1)+2) + coords(node_id,2)*dpsi(i,node_id,alpha);
            E_real(i,3*(alpha-1)+3) = E_real(i,3*(alpha-1)+3) + coords(node_id,3)*dpsi(i,node_id,alpha);
        end
    end
end

% Jacobian
det_E = zeros(ngp,1);
for i=1:ngp
    
    tmp1 = E_real(i,1:3);
    tmp2 = E_real(i,4:6);

    det_E(i,1) = norm(cross(tmp1,tmp2),2);

end

% x_gp and F_gp at gps
x_gp = zeros(ngp,3);
F_gp = zeros(ngp,3,2);

for i=1:ngp
    for p = 1:node_id
        x_gp(i,1) = x_gp(i,1) + coords(p,1)*psi(i,p);
        x_gp(i,2) = x_gp(i,2) + coords(p,2)*psi(i,p);
        x_gp(i,3) = x_gp(i,3) + coords(p,3)*psi(i,p);

        F_gp(i,1,1) = F_gp(i,1,1) + coords(p,1)*dpsi(i,p,1);
        F_gp(i,1,2) = F_gp(i,1,2) + coords(p,1)*dpsi(i,p,2);

        F_gp(i,2,1) = F_gp(i,2,1) + coords(p,2)*dpsi(i,p,1);
        F_gp(i,2,2) = F_gp(i,2,2) + coords(p,2)*dpsi(i,p,2);

        F_gp(i,3,1) = F_gp(i,3,1) + coords(p,3)*dpsi(i,p,1);
        F_gp(i,3,2) = F_gp(i,3,2) + coords(p,3)*dpsi(i,p,2);

    end 
end

x_gp_av = zeros(3,1);

for i=1:ngp
    x_gp_av(1,1) = x_gp_av(1,1) + x_gp(i,1);
    x_gp_av(2,1) = x_gp_av(2,1) + x_gp(i,2);
    x_gp_av(3,1) = x_gp_av(3,1) + x_gp(i,3);
end

x_gp_av(1,1) = x_gp_av(1,1)/ngp;
x_gp_av(2,1) = x_gp_av(2,1)/ngp;
x_gp_av(3,1) = x_gp_av(3,1)/ngp;

%residual vector
res = zeros(4,1);
fe = zeros(4,1);

for i=1:ngp
    
    for p=1:4
        
        fe(p,1) = psi(i,p)*F_gp(i,1,1) + dpsi(i,p,1)*x_gp_av(1,1);
        
    end

    for p=1:4
        res(p,1) = res(p,1) - fe(p,1)*det_E(i)*wt(i,1);
    end

end





















