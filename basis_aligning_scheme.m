clear all;

ngp = 4; % no of gps

coords = load('./input/coordinates.in');
el_conn = load('./input/element_connectivity.in');

[nel,mm] = size(el_conn);

for ie=1:nel
    
    el_conn(ie,1) = el_conn(ie,1)+1;
    el_conn(ie,2) = el_conn(ie,2)+1;
    el_conn(ie,3) = el_conn(ie,3)+1;
    el_conn(ie,4) = el_conn(ie,4)+1;

end

% coords = zeros(16,3); % x, y & z coordinates for an element
% 
% hx = 0.2;
% hy = 0.15;
% 
% for i=1:4
%     for j=1:4
% 
%         coords(i+(j-1)*3,1) = (i-1)*hx;
%         coords(i+(j-1)*3,2) = (j-1)*hy;
%         coords(i+(j-1)*3,3) = 0.0;
% 
%     end
% end
% 
% coords(6,1) = coords(6,1) - 0.2*hx;
% coords(6,2) = coords(6,2) + 0.3*hy;
% 
% coords(8,1) = coords(8,1) + 0.1*hx;
% coords(8,2) = coords(8,2) - 0.1*hy;
% 
% % connectivity matrix
% nel = 9;
% el_conn = zeros(nel,4);
% el_conn(1,:) = [1 2 5 6];
% el_conn(2,:) = [2 3 7 6];
% el_conn(3,:) = [3 4 8 7];
% el_conn(4,:) = [5 6 10 9];
% el_conn(5,:) = [10 6 7 11];
% el_conn(6,:) = [7 8 12 11];
% el_conn(7,:) = [9 10 14 13];
% el_conn(8,:) = [10 11 15 14];
% el_conn(9,:) = [11 12 16 15];

% edges list for each element
edges_list = zeros(4*nel,2);

for ie=1:nel
    edges_list(4*(ie-1)+1,:) = [el_conn(ie,1:2)];
    edges_list(4*(ie-1)+2,:) = [el_conn(ie,2:3)];
    edges_list(4*(ie-1)+3,:) = [el_conn(ie,3:4)];
    edges_list(4*(ie-1)+4,:) = [el_conn(ie,4),el_conn(ie,1)];
end

% signed edges list for each element
signs_edges_list = zeros(nel,4);

% boolean list for elements
logical_array = false(nel,1);

xi = zeros(ngp,2); % gps coordinates in master element 
wt = zeros(ngp,1); % wt of gps

if (ngp == 4)
    xi(1,1) = -0.577350269189626; xi(1,2) = -0.577350269189626;
    xi(2,1) = 0.577350269189626; xi(2,2) = -0.577350269189626;
    xi(3,1) = 0.577350269189626; xi(3,2) = 0.577350269189626;
    xi(4,1) = -0.577350269189626; xi(4,2) = 0.577350269189626;

    wt(1,1) = 1.0; wt(2,1) = 1.0; wt(3,1) = 1.0; wt(4,1) = 1.0;

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


% assigning signs to edges
ie = 1; % first element
sum = 1;

% getting the basis vector at 1st gp
E_real= zeros(6,1); 
for alpha = 1:2
    for node_id=1:4
        E_real(3*(alpha-1)+1) = E_real(3*(alpha-1)+1,1) + coords(el_conn(ie,node_id),1)*dpsi(1,node_id,alpha); 
        E_real(3*(alpha-1)+2) = E_real(3*(alpha-1)+2,1) + coords(el_conn(ie,node_id),2)*dpsi(1,node_id,alpha);
        E_real(3*(alpha-1)+3) = E_real(3*(alpha-1)+3,1) + coords(el_conn(ie,node_id),3)*dpsi(1,node_id,alpha);
    end
end

v1 = zeros(3,1); v2 = zeros(3,1); 
v3 = zeros(3,1); v4 = zeros(3,1);

for i=1:3
   v1(i,1) = coords(edges_list(4*(ie-1)+1,2),i) - coords(edges_list(4*(ie-1)+1,1),i);
   v2(i,1) = coords(edges_list(4*(ie-1)+2,2),i) - coords(edges_list(4*(ie-1)+2,1),i);
   v3(i,1) = coords(edges_list(4*(ie-1)+3,2),i) - coords(edges_list(4*(ie-1)+3,1),i);
   v4(i,1) = coords(edges_list(4*(ie-1)+4,2),i) - coords(edges_list(4*(ie-1)+4,1),i);
end

% angles b/w v and E_basis vectors
theta = zeros(4,1);
theta(1,1)  = acosd(dot(v1,E_real(1:3,1))/(norm(v1)*norm(E_real(1:3,1))));
theta(2,1)  = acosd(dot(v2,E_real(1:3,1))/(norm(v2)*norm(E_real(1:3,1))));
theta(3,1)  = acosd(dot(v3,E_real(1:3,1))/(norm(v3)*norm(E_real(1:3,1))));
theta(4,1)  = acosd(dot(v4,E_real(1:3,1))/(norm(v4)*norm(E_real(1:3,1))));

[M1,I1] = min(theta);
theta2 = 180*ones(4,1) - theta; 

[M2, I2] = min(theta2);

if (M1 <= M2)

    if (I1 == 1)
        signs_edges_list(ie,1) = +1;
        signs_edges_list(ie,2) = +2;
        signs_edges_list(ie,3) = -1;
        signs_edges_list(ie,4) = -2;
    elseif (I1 == 2)
        signs_edges_list(ie,1) = -2;
        signs_edges_list(ie,2) = +1;
        signs_edges_list(ie,3) = +2;
        signs_edges_list(ie,4) = -1;
    elseif (I1 == 3)
        signs_edges_list(ie,1) = -1;
        signs_edges_list(ie,2) = -2;
        signs_edges_list(ie,3) = +1;
        signs_edges_list(ie,4) = +2;
    elseif (I1 == 4)
        signs_edges_list(ie,1) = +2;
        signs_edges_list(ie,2) = -1;
        signs_edges_list(ie,3) = -2;
        signs_edges_list(ie,4) = +1;
    end

elseif (M1 > M2)

    if (I2 == 1)
        signs_edges_list(ie,1) = -1;
        signs_edges_list(ie,2) = -2;
        signs_edges_list(ie,3) = +1;
        signs_edges_list(ie,4) = +2;
    elseif (I2 == 2)
        signs_edges_list(ie,1) = +2;
        signs_edges_list(ie,2) = -1;
        signs_edges_list(ie,3) = -2;
        signs_edges_list(ie,4) = +1;
    elseif (I2 == 3)
        signs_edges_list(ie,1) = +1;
        signs_edges_list(ie,2) = +2;
        signs_edges_list(ie,3) = -1;
        signs_edges_list(ie,4) = -2;
    elseif (I2 == 4)
        signs_edges_list(ie,1) = -2;
        signs_edges_list(ie,2) = +1;
        signs_edges_list(ie,3) = +2;
        signs_edges_list(ie,4) = -1;
    end

end

logical_array(ie,1) = true;
iter = 1;
sum_p = sum;

total_iter = 1;

while sum < nel
    
    % pick an arbitrary edge
    r = randi([1 4],1,1);
    
    % nodes of that egde
    n1 = edges_list(4*(ie-1)+r,1);
    n2 = edges_list(4*(ie-1)+r,2);
    
    % elements associated with the edge
    [ie1,c1] = find(el_conn == n1);
    [ie2,c2] = find(el_conn == n2);
    
    count = 0;
    for i=1:length(ie1)
        % ind in the ie2 array
        ie_ind = find(ie2 == ie1(i));

        if (ie2(ie_ind) - ie ~=0)
            new_ie = ie2(ie_ind);
            count = 1;
            break;
        end
    end
    
    % for boundary element
    if (count == 0)
        new_ie = ie;
    end

    
    % if false then list of edges is unsigned
    if (logical_array(new_ie,1) == false)
        
        % summing to count the loop
        sum = sum + 1;

        % assign directions for the new element
        list_of_nodes = el_conn(new_ie,:);

        ind1 = find(list_of_nodes == n1);
        ind2 = find(list_of_nodes == n2);

        if ((ind1 == 1 && ind2 == 2) || (ind1 ==2 && ind2 ==1))
            new_r = 1;
        elseif ((ind1 == 2 && ind2 == 3) || (ind1 ==3 && ind2 ==2))
            new_r = 2;
        elseif ((ind1 == 3 && ind2 == 4) || (ind1 ==4 && ind2 ==3))
            new_r = 3;
        elseif ((ind1 == 1 && ind2 == 4) || (ind1 ==4 && ind2 ==1))
            new_r = 4;
        end
        
        new_sign_edge = -1.0*signs_edges_list(ie,r);
        signs_edges_list(new_ie,new_r) = new_sign_edge;
        
        if (new_r ==1)
            new_list_ind = [2 3 4];
        elseif (new_r == 2)
            new_list_ind = [3 4 1];
        elseif (new_r == 3)
            new_list_ind = [4 1 2];
        elseif (new_r == 4)
            new_list_ind = [1 2 3];
        end


        if (new_sign_edge == +1)
            signs_edges_list(new_ie,new_list_ind(1)) = +2;
            signs_edges_list(new_ie,new_list_ind(2)) = -1;
            signs_edges_list(new_ie,new_list_ind(3)) = -2;
        elseif (new_sign_edge == +2)
            signs_edges_list(new_ie,new_list_ind(1)) = -1;
            signs_edges_list(new_ie,new_list_ind(2)) = -2;
            signs_edges_list(new_ie,new_list_ind(3)) = +1;
        elseif (new_sign_edge == -1)
            signs_edges_list(new_ie,new_list_ind(1)) = -2;
            signs_edges_list(new_ie,new_list_ind(2)) = +1;
            signs_edges_list(new_ie,new_list_ind(3)) = +2;
        elseif (new_sign_edge == -2)
            signs_edges_list(new_ie,new_list_ind(1)) = +1;
            signs_edges_list(new_ie,new_list_ind(2)) = +2;
            signs_edges_list(new_ie,new_list_ind(3)) = -1;
        end
        
        % assigning next element as the current
        ie = new_ie;
        logical_array(ie,1) = true;
    end

    % if all adjoining elements already have signed edges
    if (sum - sum_p == 0)
        iter = iter + 1;
    else
        iter = 1;
        sum_p = sum;
    end

    if (iter > 4)
        
       ie_list = find(logical_array == true);
       size_list = length(ie_list);
       ie_ind  = randi([1 size_list],1,1);
       ie = ie_list(ie_ind);

    end

    total_iter = total_iter + 1;
    

end

% correcting the connectivity matrix

new_el_conn = zeros(nel,4);

for ie=1:nel

    if (signs_edges_list(ie,1) == 1)        
        new_el_conn(ie,1) = el_conn(ie,1);
        new_el_conn(ie,2) = el_conn(ie,2);
        new_el_conn(ie,3) = el_conn(ie,3);
        new_el_conn(ie,4) = el_conn(ie,4);
    elseif (signs_edges_list(ie,1) == 2)        
        new_el_conn(ie,1) = el_conn(ie,4);
        new_el_conn(ie,2) = el_conn(ie,1);
        new_el_conn(ie,3) = el_conn(ie,2);
        new_el_conn(ie,4) = el_conn(ie,3);
    elseif (signs_edges_list(ie,1) == -1)        
        new_el_conn(ie,1) = el_conn(ie,3);
        new_el_conn(ie,2) = el_conn(ie,4);
        new_el_conn(ie,3) = el_conn(ie,1);
        new_el_conn(ie,4) = el_conn(ie,2);
    elseif (signs_edges_list(ie,1) == -2)        
        new_el_conn(ie,1) = el_conn(ie,2);
        new_el_conn(ie,2) = el_conn(ie,3);
        new_el_conn(ie,3) = el_conn(ie,4);
        new_el_conn(ie,4) = el_conn(ie,1);
    end

end


% writing a new file
fileID1= fopen('element_connectivity.in','w');

for i=1:nel
    fprintf(fileID1,'%12d %12d %12d %12d\r\n',new_el_conn(i,1)-1,new_el_conn(i,2)-1,new_el_conn(i,3)-1,new_el_conn(i,4)-1); 
end

fclose(fileID1);



% % residual vector global
% res_g = zeros(9,1);
% 
% % loop over the elements
% for ie=1:nel
% 
% 
%     %residual vector
%     res_e = zeros(4,1);
%     fe = zeros(4,1);
% 
%     % real basis vectors
%     E_real= zeros(ngp,6);
%     for i=1:ngp    
%         for alpha = 1:2
%             for node_id=1:4
%                 E_real(i,3*(alpha-1)+1) = E_real(i,3*(alpha-1)+1) + coords(el_conn(ie,node_id),1)*dpsi(i,node_id,alpha); 
%                 E_real(i,3*(alpha-1)+2) = E_real(i,3*(alpha-1)+2) + coords(el_conn(ie,node_id),2)*dpsi(i,node_id,alpha);
%                 E_real(i,3*(alpha-1)+3) = E_real(i,3*(alpha-1)+3) + coords(el_conn(ie,node_id),3)*dpsi(i,node_id,alpha);
%             end
%         end
%     end
% 
%     % Jacobian
%     det_E = zeros(ngp,1);
%     for i=1:ngp
% 
%         tmp1 = E_real(i,1:3);
%         tmp2 = E_real(i,4:6);
% 
%         % cross(tmp1,tmp2)
%         det_E(i,1) = norm(cross(tmp1,tmp2));
% 
%     end
% 
%     % x_gp and F_gp at gps
%     x_gp = zeros(ngp,3);
%     F_gp = zeros(ngp,3,2);
% 
%     for i=1:ngp
%         for p=1:4
% 
%             x_gp(i,1) = x_gp(i,1) + coords(el_conn(ie,p),1)*psi(i,p);
%             x_gp(i,2) = x_gp(i,2) + coords(el_conn(ie,p),2)*psi(i,p);
%             x_gp(i,3) = x_gp(i,3) + coords(el_conn(ie,p),3)*psi(i,p);
% 
%             F_gp(i,1,1) = F_gp(i,1,1) + coords(el_conn(ie,p),1)*dpsi(i,p,1);
%             F_gp(i,1,2) = F_gp(i,1,2) + coords(el_conn(ie,p),1)*dpsi(i,p,2);
% 
%             F_gp(i,2,1) = F_gp(i,2,1) + coords(el_conn(ie,p),2)*dpsi(i,p,1);
%             F_gp(i,2,2) = F_gp(i,2,2) + coords(el_conn(ie,p),2)*dpsi(i,p,2);
% 
%             F_gp(i,3,1) = F_gp(i,3,1) + coords(el_conn(ie,p),3)*dpsi(i,p,1);
%             F_gp(i,3,2) = F_gp(i,3,2) + coords(el_conn(ie,p),3)*dpsi(i,p,2);
% 
%         end
%     end
% 
%     % sum over the gps
%     for i=1:ngp
% 
%         for p=1:4
% 
%             fe(p,1) = psi(i,p)*F_gp(i,1,1) + dpsi(i,p,1)*x_gp(i,1);
% 
%         end
% 
%         for p=1:4
%             res_e(p,1) = res_e(p,1) - fe(p,1)*det_E(i)*wt(i,1);
%         end
% 
%     end
% 
%     % res_e
% 
%     % assembly of residual vector
%     % el_conn(ie,:)
%     res_g(el_conn(ie,:),1) = res_g(el_conn(ie,:),1) + res_e;
% 
% end
