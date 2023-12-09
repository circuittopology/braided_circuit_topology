function [br, slicedbraidcoordinates, idxlist] = data_to_braid(braidlist, braidcoordinates, nr_beads)
import braidlab.*

bl = length(braidlist);

%%%%%%%% End-to-end vectors and planar orthonormal vectors %%%%%%%%%%%%%%%%
R = zeros(bl,3);
endpoints = zeros(2,bl,3);

for i=1:bl
    endpoints(1,i,:) = braidcoordinates(1,:,i);
    endpoints(2,i,:) = braidcoordinates(end,:,i);
    R(i,:) = braidcoordinates(end,:,i) - braidcoordinates(1,:,i);
end
meanstart = reshape(mean(endpoints(1,:,:)),[1 3]);
meanend = reshape(mean(endpoints(2,:,:)),[1 3]);

Rm = meanend-meanstart; % mean end-to end vector; becomes normal vector for scanning plane
Rnormal = Rm/norm(Rm); % Normalise R
Q = (reshape(endpoints(1,1,:),[1 3])-meanstart)/norm(reshape(endpoints(1,1,:),[1 3])-meanstart); % Normalised vector from meanstart to first terminal
portho1 = (Q - dot(Q,Rnormal)*Rnormal)/norm(Q - dot(Q,Rnormal)*Rnormal);
portho2 = cross(portho1,Rnormal)/norm(cross(portho1,Rnormal));

A = transpose([portho1; portho2; Rnormal]);
transformedbraidcoordinates = zeros(nr_beads,3,bl);

for i=1:bl
    for j=1:nr_beads
    transformedbraidcoordinates(j,:,i) = transpose(A*transpose(braidcoordinates(j,:,i))) ;
    end
end
meanstarttransform = transpose(A*transpose(meanstart));
meanendtransform = transpose(A*transpose(meanend));

%%%%%%%%%%%% Slice constant z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transformedendpoints = [reshape(transformedbraidcoordinates(1,3,:),[1 bl]); reshape(transformedbraidcoordinates(end,3,:),[1 bl])];
Zend = min(transformedendpoints,[],'all');
Zstart = max(transformedendpoints,[],'all');

dz = abs(Zstart-Zend)/(2*nr_beads-1);
Lslice = length(Zstart:-dz:Zend);
count =1;

slicedbraidcoordinates = [];
for slice=Zstart:-dz:Zend
    for j = 1:bl
        Xarray = transformedbraidcoordinates(:,1,j);
        Yarray = transformedbraidcoordinates(:,2,j);
        Zarray = transformedbraidcoordinates(:,3,j);
        diffValues = Zarray-slice;
        diffValues(diffValues >0) =-inf;
        [~, indexOfMax] = max(diffValues);
        if(indexOfMax==1)
            break;
        end
        yminus = Yarray(indexOfMax);
        yplus = Yarray(indexOfMax-1);
        xminus = Xarray(indexOfMax);
        xplus = Xarray(indexOfMax-1);
        zminus = Zarray(indexOfMax);
        zplus = Zarray(indexOfMax-1);
        x0 = xminus + ((slice-zminus)/(zplus-zminus))*(xplus-xminus);
        y0 = yminus + ((slice-zminus)/(zplus-zminus))*(yplus-yminus);
        slicedbraidcoordinates(count,1,j) = x0;
        slicedbraidcoordinates(count,2,j) = y0;
        slicedbraidcoordinates(count,3,j) = slice;
    end
    count = count+1;
end

lastline = length(slicedbraidcoordinates);
removeFirst = [];
removeLast = [];
for i=1:bl
    tmptotal = sum(slicedbraidcoordinates(:,:,i),2);
    removeFirst = [removeFirst find(tmptotal~=0, 1, 'first')];
    removeLast = [removeLast find(tmptotal~=0, 1, 'last')];
end

if max(removeFirst)~=1 % If there are leading zero rows, remove them
    slicedbraidcoordinates(1:(max(removeFirst)-1),:,:) =[];
end
if min(removeLast)<lastline % If there are trailing zero rows, remove them but take care of shifted index from initial removal of leading rows
    slicedbraidcoordinates = slicedbraidcoordinates(1:(end-abs(lastline-min(removeLast))),:,:);
end

idxlist=[];
for i=1:bl
    idxlist = [idxlist find(all(slicedbraidcoordinates(:,:,i) == 0,2))];
end
if isempty(idxlist)==0
    slicedbraidcoordinates(max(idxlist),:,:) = [];
end


br = braid(slicedbraidcoordinates(:,1:2,:));
end