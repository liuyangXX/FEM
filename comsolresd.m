clear all
close all

fp = fopen('F:\FEM\model_normal.mphtxt','r');
for i = 1:18
    tline = fgets(fp);
end
num_nodes = fscanf(fp,'%d # number of mesh points\n',[1,1]);
pts_ind = fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);
xy = fscanf(fp,'%lf %lf \n',[2,num_nodes]);
xy = xy';
x = xy(:,1);
y = xy(:,2);

for i = 1:7
    fgets(fp);
end
num_vtx_ns = fscanf(fp,'%d # number of nodes per element\n',[1,1]);
num_vtx_ele = fscanf(fp,'%d # number of elements\n',[1.1]);
fgets(fp);
vtx = fscanf(fp,'%d \n',[1,num_vtx_ele]);
vtx = vtx'+1;

num_vtx_ele2 = fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
vtx2 = fscanf(fp,'%d \n',[1,num_vtx_ele2]);
vtx2 = vtx2'+1;

for i = 1:5
    fgets(fp);
end
num_bdr_ns = fscanf(fp,'%d # number of nodes per element\n',[1,1]);
num_bdr_ele = fscanf(fp,'%d # number of elements\n',[1,1]);
fgets(fp);
p = fscanf(fp,'%d %d\n',[2,num_bdr_ele]);
p = p'+1;
NDPP = p(:,1);

num_entity = fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
entity = fscanf(fp,'%d \n',[1,num_entity]);
entity = entity'+1;

for i = 1:5
    fgets(fp);
end
ns_per_element = fscanf(fp,'%d # number of nodes per element\n',[1,1]);
NE = fscanf(fp,'%d # number of elements\n',[1,1]);
fgets(fp);
NL = fscanf(fp,'%d %d %d\n',[3,NE]);
NL = NL'+1;                         %List of nodes for each element

num_domain = fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
Domain = fscanf(fp,'%d\n',[1,num_domain]);
Domain = Domain';
fclose(fp);

NL(:,4) = Domain;
NP = length(NDPP);              %边界上节点的个数 存入NP
NF = num_nodes-NP;              %自由节点的个数
%%%设定存储的系数矩阵
Ce = zeros(3,3);
P = zeros(3,1);
Q = zeros(3,1);
C = zeros(num_nodes,num_nodes);
XL = zeros(1,3);
YL = zeros(1,3);
for i=1:NE                               
    for j=1:3
        k = NL(i,j);
        XL(j) = x(k);
        YL(j) = y(k);
    end

    P(1)=YL(2)-YL(3);
    P(2)=YL(3)-YL(1);
    P(3)=YL(1)-YL(2);

    Q(1)=XL(3)-XL(2);
    Q(2)=XL(1)-XL(3);
    Q(3)=XL(2)-XL(1);
    A=1/2*(P(2)*Q(3)-P(3)*Q(2)); %area of elements
    
    for m=1:3
        for n=1:3
           Ce(m,n) = 1/(4*A)*(P(m)*P(n)+Q(m)*Q(n));     %element coefficient matrix
        end
    end 
    for j=1:3
        IR = NL(i,j);
        for t=1:3
            IC=NL(i,t);
            C(IR,IC) = C(IR,IC) + Ce(j,t);           %global coefficient matrix
        end
    end 
end

























