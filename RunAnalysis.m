function RunAnalysis()
%{
Supporting Information to: "Harnessing Geometric Frustration to Form Band
Gaps in Acoustic Networks"

By: Pai Wang, Yue Zheng, Matheus C. Fernandes,
Yushen Sun, Kai Xu, Sijie Sun, Sung Hoon Kang,
Vincent Tournat, Katia Bertoldi.

doi:10.1103/PhysRevLett.118.084302

Harvard University, February 2017

Notes:
Run the function as is and use ABAQUS input files for geometry
definitions. Please refer to the .inp for the format requested. This code
will only work with 1D geometries with specified Block-type boundary
conditions. For the numerical deptiction of this code, please refer to the
supplemental section of the article.

Corresponding code author: m@fer.me
%}

clc;close all;clear vars
answer = inputdlg({'Enter input file path (.inp)','Enter file path for k-vector','Number of Modes Desired'},'answer',1,{'test.inp','test.txt','20'});
inputfile=answer{1};
KVect=importdata(answer{2});
mode=str2double(answer{3});

[NodeCoor,SymNodes,KGlobal,MGlobal]=FEMGener(inputfile,1,1);
for i=1:length(KVect(:,1))
    w(i,:)= Program(KVect(i,:),mode,KGlobal,MGlobal,NodeCoor,SymNodes);
end

%% Plotting
figure
hold all
for i=1:mode
    plot(1:length(KVect(:,1)),w(:,i),'k','linewidth',2)
end
xlabel('Reduced Wave Vector, k','fontsize',14,'fontweight','bold')
ylabel('Normalized Frequency, \Omega','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold','linewidth',1.5)
axis([1 length(KVect(:,1)) min(min(w)) max(max(w))])
box on
axis square

end

function [NodeCoor,SymNodes,KGlobal,MGlobal]=FEMGener(NameOfInputFile,Density,cSquared)
%% Reading Input files and defining Global Stiffness and Mass Matrices
[NodeCoor,ElemConnectivity,NumElem,SymNodes]=ReadInputAbaqus(NameOfInputFile);
fclose all;

%Number of Nodes
NumNode=size(NodeCoor);
NumNode=NumNode(1);

%Creating a Zero square matrix with 2 times to account for the X&Y
%component of each node.
KGlobal=zeros(NumNode);
MGlobal=zeros(NumNode);

for i=1:NumElem
    %Node 1 Global Number
    Node1GlobNum=ElemConnectivity(i,2);
    %Node 2 Global Number
    Node2GlobNum=ElemConnectivity(i,3);
    
    %Node1 [x,y] coordinates
    NodeC1=NodeCoor(Node1GlobNum,[2,3]);
    %Node2 [x,y] coordinates
    NodeC2=NodeCoor(Node2GlobNum,[2,3]);
    
    [KLocal,MLocal]=Truss2D(NodeC1,NodeC2,cSquared,Density);
    
    index=[Node1GlobNum Node2GlobNum];
    
    KGlobal(index,index)=KGlobal(index,index)+KLocal;
    MGlobal(index,index)=MGlobal(index,index)+MLocal;
    
end

end

function [NodeCoor,ElemConnectivity,numelements,SymNodes]=ReadInputAbaqus(NameTextFile)
ind=0;
fid=fopen(NameTextFile,'r');

while ~feof(fid)
    
    lineread= lower(fgetl(fid));
    
    if numel(lineread)>=5
        switch lineread
            case '*node'
                NodeCoor=fscanf(fid,'%f, %f, %f',[3 inf])';
        end
    end
    if numel(lineread)>=8
        switch lineread(1:8)
            case '*element'
                ElemConnectivity=fscanf(fid,'%f, %f, %f',[3 inf])';
        end
    end
    if numel(lineread)>=5
        switch lineread(1:5)
            
            case '*nset'
                ind=ind+1;
                SymNodes(ind,:)=fscanf(fid,'%f, %f',[2 inf]);
        end
    end
end

fclose(fid);
NodeCoor=sortrows(NodeCoor);

numelements=size(ElemConnectivity);
numelements=numelements(1);

for i=1:numelements
    elementcoor(i,1)=ElemConnectivity(i,1);
    elementcoor(i,2)=NodeCoor(ElemConnectivity(i,2),2);
    elementcoor(i,3)=NodeCoor(ElemConnectivity(i,2),3);
    elementcoor(i,4)=NodeCoor(ElemConnectivity(i,3),2);
    elementcoor(i,5)=NodeCoor(ElemConnectivity(i,3),3);
end

end

function [w]= Program(KVect,mode,KGlobal,MGlobal,NodeCoor,SymNodes)
%% Applying Periodic Boundary Coditions
%Number of symetry nodes
numSymNodes=size(SymNodes);
numSymNodes=numSymNodes(1);

for MA=1:numSymNodes
    clear SIndex MIndex
    MIndex=find(NodeCoor(:,1)==SymNodes(MA,1));
    SIndex=find(NodeCoor(:,1)==SymNodes(MA,2));
    
    MasterXCoor=NodeCoor(MIndex,2);
    MasterYCoor=NodeCoor(MIndex,3);
    SlaveXCoor=NodeCoor(SIndex,2);
    SlaveYCoor=NodeCoor(SIndex,3);
    
    R=[MasterXCoor-SlaveXCoor; MasterYCoor-SlaveYCoor];
    
    BCVect=KVect*R;
    
    KGlobal(MIndex,:)=KGlobal(MIndex,:)+exp(1i*BCVect)*KGlobal(SIndex,:);
    KGlobal(:,MIndex)=KGlobal(:,MIndex)+exp(-1i*BCVect)*KGlobal(:,SIndex);
    
    MGlobal(MIndex,:)=MGlobal(MIndex,:)+exp(1i*BCVect)*MGlobal(SIndex,:);
    MGlobal(:,MIndex)=MGlobal(:,MIndex)+exp(-1i*BCVect)*MGlobal(:,SIndex);
    
    KGlobal(SIndex,:)=[];
    KGlobal(:,SIndex)=[];
    
    MGlobal(SIndex,:)=[];
    MGlobal(:,SIndex)=[];
    
    for AS=MA:length(SymNodes)
        for AF=1:2
            if SymNodes(AS,AF)==SymNodes(MA,2)
                SymNodes(AS,AF)=SymNodes(MA,1);
            end
        end
    end
    NodeCoor(SIndex,:)=[];
end
%% Solving Eigenvalue problem
w=eig(KGlobal,MGlobal);
w=real(sqrt(w))/(2*pi);
w=sort(w);
w=w(1:mode);
end

function [Klocal,Mlocal]=Truss2D(Coor1,Coor2,cSquared,Density)

%Subtratcting coordinates
XDiff=Coor2(1)-Coor1(1); %x-coordinate
YDiff=Coor2(2)-Coor1(2); %y-coordinate

%Finding length of elemet using Pythagorean
l=sqrt(XDiff^2+YDiff^2);

%local element Stiffness matrix
Klocal=(cSquared/l)*[1 -1;-1 1];
%local element mass matrix
Mlocal=(Density*l)/6*[2 1; 1 2];
end