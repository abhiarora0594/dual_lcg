clear all;

ngp = 4; % no of gps
coords = zeros(9,3); % x, y & z coordinates for an element

hx = 0.2;
hy = 0.15;

for i=1:3
    for j=1:3

        coords(i+(j-1)*3,1) = (i-1)*hx;
        coords(i+(j-1)*3,2) = (j-1)*hy;
        coords(i+(j-1)*3,3) = 0.0;
    
    end
end

% coords(5,1) = coords(5,1) - 0.2*hx;
% coords(5,2) = coords(6,1) + 0.3*hy;

% connectivity matrix
nel = 4;
el_conn = zeros(nel,4);
el_conn(1,:) = [1 2 5 4];
el_conn(2,:) = [2 3 6 5];
el_conn(3,:) = [4 5 8 7];
el_conn(4,:) = [6 9 8 5];

xi = zeros(ngp,2); % gps coordinates in master element 
wt = zeros(ngp,1); % wt of gps

if (ngp == 4)
    xi(1,1) = -0.577350269189626; xi(1,2) = -0.577350269189626;
    xi(2,1) = 0.577350269189626; xi(2,2) = -0.577350269189626;
    xi(3,1) = 0.577350269189626; xi(3,2) = 0.577350269189626;
    xi(4,1) = -0.577350269189626; xi(4,2) = 0.577350269189626;

    wt(1,1) = 1.0; wt(2,1) = 1.0; wt(3,1) = 1.0; wt(4,1) = 1.0;

elseif (ngp == 16)
    
    scal1 = sqrt((3/7) - (2/7)*sqrt(6/5));
    scal2 = sqrt((3/7) + (2/7)*sqrt(6/5));
    scal3 = (18 + sqrt(30))/36;
    scal4 = (18 - sqrt(30))/36;

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

    scal1 = 0;
    scal2 = (1/3)*sqrt(5 - 2*sqrt(10/7));
    scal3 = 1/3*sqrt(5 + 2*sqrt(10/7));
    
    scal4 = 128/225;
    scal5 = (322+13*sqrt(70))/900;
    scal6 = (322-13*sqrt(70))/900;
    
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

psi_g = zeros(ngp,4);
dpsi_g = zeros(ngp,4,2);
  
% shape function values at different gps
for i=1:ngp
    psi(i,1) =  0.25*(1.0-xi(i,1))*(1.0-xi(i,2));
    psi(i,2) =  0.25*(1.0+xi(i,1))*(1.0-xi(i,2));
    psi(i,3) =  0.25*(1.0+xi(i,1))*(1.0+xi(i,2));
    psi(i,4) =  0.25*(1.0-xi(i,1))*(1.0+xi(i,2));
    
    psi_g(i,1) =  0.25*(xi(i,1).^2 - xi(i,1))*(xi(i,2).^2 - xi(i,2));
    psi_g(i,2) =  0.25*(xi(i,1).^2 + xi(i,1))*(xi(i,2).^2 - xi(i,2));
    psi_g(i,3) =  0.25*(xi(i,1).^2 + xi(i,1))*(xi(i,2).^2 + xi(i,2));
    psi_g(i,4) =  0.25*(xi(i,1).^2 - xi(i,1))*(xi(i,2).^2 + xi(i,2));

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

    dpsi_g(i,1,1) =  0.25*(2.0*xi(i,1) - 1.0)*(xi(i,2).^2 - xi(i,2));
    dpsi_g(i,2,1) =  0.25*(2.0*xi(i,1) + 1.0)*(xi(i,2).^2 - xi(i,2));
    dpsi_g(i,3,1) =  0.25*(2.0*xi(i,1) + 1.0)*(xi(i,2).^2 + xi(i,2));
    dpsi_g(i,4,1) =  0.25*(2.0*xi(i,1) - 1.0)*(xi(i,2).^2 + xi(i,2));

    dpsi_g(i,1,2) =  0.25*(xi(i,1).^2 - xi(i,1))*(2.0*xi(i,2) - 1.0);
    dpsi_g(i,2,2) =  0.25*(xi(i,1).^2 + xi(i,1))*(2.0*xi(i,2) - 1.0);
    dpsi_g(i,3,2) =  0.25*(xi(i,1).^2 + xi(i,1))*(2.0*xi(i,2) + 1.0);
    dpsi_g(i,4,2) =  0.25*(xi(i,1).^2 - xi(i,1))*(2.0*xi(i,2) + 1.0);

end

% residual vector global
res_g = zeros(9,1);

% loop over the elements
for ie=1:nel


    %residual vector
    res_e = zeros(4,1);
    fe = zeros(4,1);

    % real basis vectors
    E_real= zeros(ngp,6);
    for i=1:ngp    
        for alpha = 1:2
            for node_id=1:4
                E_real(i,3*(alpha-1)+1) = E_real(i,3*(alpha-1)+1) + coords(el_conn(ie,node_id),1)*dpsi(i,node_id,alpha); 
                E_real(i,3*(alpha-1)+2) = E_real(i,3*(alpha-1)+2) + coords(el_conn(ie,node_id),2)*dpsi(i,node_id,alpha);
                E_real(i,3*(alpha-1)+3) = E_real(i,3*(alpha-1)+3) + coords(el_conn(ie,node_id),3)*dpsi(i,node_id,alpha);
            end
        end
    end
    
    % Jacobian
    det_E = zeros(ngp,1);
    for i=1:ngp
        
        tmp1 = E_real(i,1:3);
        tmp2 = E_real(i,4:6);
        
        % cross(tmp1,tmp2)
        det_E(i,1) = norm(cross(tmp1,tmp2));
    
    end
    
    % x_gp and F_gp at gps
    x_gp = zeros(ngp,3);
    F_gp = zeros(ngp,3,2);
    
    for i=1:ngp
        for p=1:4
            
            x_gp(i,1) = x_gp(i,1) + coords(el_conn(ie,p),1)*psi(i,p);
            x_gp(i,2) = x_gp(i,2) + coords(el_conn(ie,p),2)*psi(i,p);
            x_gp(i,3) = x_gp(i,3) + coords(el_conn(ie,p),3)*psi(i,p);
    
            F_gp(i,1,1) = F_gp(i,1,1) + coords(el_conn(ie,p),1)*dpsi(i,p,1);
            F_gp(i,1,2) = F_gp(i,1,2) + coords(el_conn(ie,p),1)*dpsi(i,p,2);
    
            F_gp(i,2,1) = F_gp(i,2,1) + coords(el_conn(ie,p),2)*dpsi(i,p,1);
            F_gp(i,2,2) = F_gp(i,2,2) + coords(el_conn(ie,p),2)*dpsi(i,p,2);
    
            F_gp(i,3,1) = F_gp(i,3,1) + coords(el_conn(ie,p),3)*dpsi(i,p,1);
            F_gp(i,3,2) = F_gp(i,3,2) + coords(el_conn(ie,p),3)*dpsi(i,p,2);
        
        end
    end
    
    % sum over the gps
    for i=1:ngp
        
        for p=1:4
            
            fe(p,1) = psi(i,p)*F_gp(i,1,1) + dpsi(i,p,1)*x_gp(i,1);
            
        end
    
        for p=1:4
            res_e(p,1) = res_e(p,1) - fe(p,1)*det_E(i)*wt(i,1);
        end
    
    end

    % res_e

    % assembly of residual vector
    % el_conn(ie,:)
    res_g(el_conn(ie,:),1) = res_g(el_conn(ie,:),1) + res_e;

end
