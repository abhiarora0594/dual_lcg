clc; 
clear all;
close all;

coord = load('coordinates.txt');
el_conn = load('connectivity.txt');
bound_tmp = load('bound_nodes.txt');

numofzeros =8;

bound_tmp2= reshape(bound_tmp,[],1);

bound_nodes=zeros(size(bound_tmp2,1)-numofzeros,1);
count = 0;

for i=1:size(bound_tmp2,1)
   if (bound_tmp2(i,1)>0) 
        count= count +1 ;
        bound_nodes(count,1) = bound_tmp2(i,1);  
   end
end

Bound_nodes = zeros(size(bound_nodes,1),2);

theta = 0;

% correcting the orientation of nodes to be in anticlockwise sense
for i=1:size(bound_nodes,1)
    
    Bound_nodes(i,1) = bound_nodes(i,1);
    
    Bound_nodes(i,2) = atan2(coord(bound_nodes(i,1),2),coord(bound_nodes(i,1),3));
    
end

Bound_nodes_sorted = sortrows(Bound_nodes,2);

for i=1:size(bound_nodes,1)
   
    bound_nodes(i,1) = Bound_nodes_sorted(i,1);
end


ndof = size(coord,1);
nel = size(el_conn,1);
nbel = size(bound_nodes,1);

coord = sortrows(coord);
el_conn = sortrows(el_conn);

fileID1= fopen('coordinates.in','w');
fileID2= fopen('element_connectivity.in','w');
fileID3= fopen('bound_nodes.in','w');
fileID4= fopen('mesh_info.in','w');

fprintf(fileID4,'%12d %12d %12d\r\n',ndof,nel,nbel);

for i=1:ndof
 fprintf(fileID1,'%12.12f %12.12f %12.12f\r\n',coord(i,2),coord(i,3),0.0);
end

for i=1:nel
   fprintf(fileID2,'%12d %12d %12d %12d\r\n',el_conn(i,2)-1,el_conn(i,3)-1,el_conn(i,4)-1,el_conn(i,5)-1); 
end

for i=1:nbel
   fprintf(fileID3,'%12d\r\n',bound_nodes(i,1)-1); 
end

fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
fclose(fileID4);

