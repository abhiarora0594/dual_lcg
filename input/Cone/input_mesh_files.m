clc
% close all
clear all

aa = load('nodes.txt');
b = load('connectivity.txt');
Boundary_nodes=load('boundary.txt');
Boundary_nodes2 = load('boundary2.txt');
% initial_guess = load('./hat_shape/initial_guess.txt');

[m,n] = size(Boundary_nodes);
Boundary_nodes = reshape(Boundary_nodes,m*n,1);
ind = find(Boundary_nodes==0);
Boundary_nodes(ind) = [];

[m,n] = size(Boundary_nodes2);
Boundary_nodes2 = reshape(Boundary_nodes2,m*n,1);
ind = find(Boundary_nodes2==0);
Boundary_nodes2(ind) = [];

[m,n] = size(aa);
a = zeros(m,n);
a(:,1:2) = aa(:,1:2); % 1 node no 2 x-coord
a(:,3) = aa(:,4); % 3 is y-coord but in abaqus it was z
a(:,4) = aa(:,3); % 4 is z-coord but in abaqus it was y

maxx = max(a(:,1));

a_full = zeros(maxx,1);

[m,n] = size(a);
%ndof = m*3;

for i=1:m
    a_full(a(i,1),1) = i; 
end

for i=1:m
    a(i,1) = a_full(a(i,1),1);
end

[m,n] = size(b);

for i=1:m
    
    b(i,1)=i;
    b(i,2)=a_full(b(i,2),1); 
    b(i,3)=a_full(b(i,3),1); 
    b(i,4)=a_full(b(i,4),1); 
    b(i,5)=a_full(b(i,5),1); 
end


nb=size(Boundary_nodes,1);

for i=1:nb
    Boundary_nodes(i,1) = a_full(Boundary_nodes(i,1),1);
end

for i=1:nb
    Boundary_nodes(i,2)=a(Boundary_nodes(i,1),2);
    Boundary_nodes(i,3)=a(Boundary_nodes(i,1),3); 
    % Boundary nodes contains X1 and X2 co-ordinates of all nodes at boundary in its 2nd and 3rd column
    % to arrange all boundary nodes in a sequential manner, one nodes
    % next to another nodes by sorting it based on angles.
    if (Boundary_nodes(i,2)<0&&Boundary_nodes(i,3)>=0) % X1, X2 in second quadrant
        Boundary_nodes(i,4)=180+atand(Boundary_nodes(i,3)/Boundary_nodes(i,2)); % 4th column forms angle/slope tan
    else
         if (Boundary_nodes(i,2)<0&&Boundary_nodes(i,3)<0) % X1 and X2 in 3rd quadrant
            Boundary_nodes(i,4)=-180+atand(Boundary_nodes(i,3)/Boundary_nodes(i,2));
         else
            Boundary_nodes(i,4)=atand(Boundary_nodes(i,3)/Boundary_nodes(i,2));
         end
    end
end


Boundary_nodes=sortrows(Boundary_nodes,4);

Boundary_elems = zeros(nb,2);

for i=1:nb-1
    Boundary_elems(i,1) = Boundary_nodes(i,1);
    Boundary_elems(i,2) = Boundary_nodes(i+1,1);
end

Boundary_elems(nb,1) = Boundary_nodes(nb,1);
Boundary_elems(nb,2) = Boundary_nodes(1,1);

% fileID1= fopen('coordinates.in','w');
% fileID2= fopen('element_connectivity.in','w');
% fileID3= fopen('boundary_elems.in','w');
% fileID4= fopen('mesh_info.in','w');
fileID5=fopen('initial_guess.in','w');

% for i=1:size(a,1)
%     fprintf(fileID1,'%12.16f %12.16f %12.16f\r\n',a(i,2),a(i,3),0.0); 
% end
% 
% for i=1:size(b,1)
%     fprintf(fileID2,'%12d %12d %12d %12d\r\n',b(i,2)-1,b(i,3)-1,b(i,4)-1,b(i,5)-1); 
% end
% 
% for i=1:nb
%     fprintf(fileID3,'%12d %12d \r\n',Boundary_elems(i,1)-1,Boundary_elems(i,2)-1); 
% end
% 
% fprintf(fileID4,'%12d %12d %12d\r\n',size(a,1),size(b,1),size(Boundary_elems,1));

for i=1:size(a,1)
    fprintf(fileID5,'%12.16f %12.16f %12.16f\r\n',a(i,2),a(i,3),a(i,4)); 
end
    
% figure;
% for i=1:100%size(b,1)
%     
%     xx = [a(b(i,2),2),a(b(i,3),2),a(b(i,4),2),a(b(i,5),2),a(b(i,2),2)];
%     yy = [a(b(i,2),3),a(b(i,3),3),a(b(i,4),3),a(b(i,5),3),a(b(i,2),3)];
%     plot(xx,yy,'-r');
%     hold on;
%     
% end

% fclose(fileID1);
% fclose(fileID2);
% fclose(fileID3);
% fclose(fileID4);
fclose(fileID5);


