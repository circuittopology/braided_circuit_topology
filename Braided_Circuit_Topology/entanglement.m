clear all; clc; close all;
import braidlab.*
prop('BraidPlotDir','tb')

fpar = fopen('main_batch.topo','r');
for i=1:6
    tline = fgets(fpar);
end
A = fgets(fpar);
N = str2double(regexp(A,'\d*','Match'));

fk = fopen('ktmp.txt');
kappa = str2double(fgets(fk));
str = ['Bending strength:' num2str(kappa,'%3.1f')];
display(str)

strGR = ['RG_' num2str(kappa,'%3.1f') '.txt'];
fGR = fopen(strGR,'r');

RG = zeros(N,1);
for i=1:(3*N)
    RGtmp = str2double(fgets(fGR));
    if mod(i,3)==0
        RG(floor(i/3),1)=RGtmp;
    end
end
fullmeasure = zeros(N,6);
fullmeasure(:,4) = RG;

errors = 0;
f = waitbar(0, 'Starting...');
for j=1:N
    try
    filename = ['dump.ppa_' num2str(j) '.lammpstrj'];
    fid = fopen(filename,'r');
    for k=1:5 % Skip unnecessary lines
        tline = fgets(fid);
    end
    %------------------------ Get the box sizes -------------------------------
    
    xsize = str2num(fgets(fid));
    ysize = str2num(fgets(fid));
    zsize = str2num(fgets(fid));
    tline = fgets(fid);
    
    %------------------------ Get the coordinates -----------------------------
    coordinates = [];
    
    while ischar(tline)
        coordinates = [coordinates; str2num(tline)];
        tline = fgets(fid);
    end
    fclose(fid);
    
    nr_polymers = coordinates(end,2);
    nr_beads = coordinates(end,1)/nr_polymers;
    pol = coordinates(1:nr_beads,3:5) + [coordinates(1:nr_beads,6)*abs(xsize(1)-xsize(2)) coordinates(1:nr_beads,7)*abs(ysize(1)-ysize(2)) coordinates(1:nr_beads,8)*abs(zsize(1)-zsize(2))];
    
    for i=2:nr_polymers
        pol(:,:,i)= coordinates((1+nr_beads*(i-1)):(nr_beads*i),3:5) + [coordinates((1+nr_beads*(i-1)):(nr_beads*i),6)*abs(xsize(1)-xsize(2)) coordinates((1+nr_beads*(i-1)):(nr_beads*i),7)*abs(ysize(1)-ysize(2)) coordinates((1+nr_beads*(i-1)):(nr_beads*i),8)*abs(zsize(1)-zsize(2))];
    end
    cset = nchoosek(1:nr_polymers,2);
    cr = 0;
        for i=1:(length(nchoosek(1:nr_polymers,2)))
            pol1 = cset(i,1);
            pol2 = cset(i,2);
            v1 = [pol(end,1,pol1)-pol(1,1,pol1) pol(end,2,pol1)-pol(1,2,pol1) pol(end,3,pol1)-pol(1,3,pol1)];
            v1 = v1/norm(v1);
            v2 = [pol(end,1,pol2)-pol(1,1,pol2) pol(end,2,pol2)-pol(1,2,pol2) pol(end,3,pol2)-pol(1,3,pol2)];
            v2 = v2/norm(v2);
            cr = cr+(dot(v1,v2)^2);
        end
    cr = (2/(nr_polymers*(nr_polymers-1)))*cr;
    fullmeasure(j,5)=cr;
    
    fullbraid = compact(data_to_braid(1:nr_polymers,pol,nr_beads));
    plot(fullbraid);
    fullmeasure(j,1) = writhe(fullbraid);
    fullmeasure(j,2) = complexity(fullbraid);
    fullmeasure(j,3) = fullbraid.length;

    tnt = train(fullbraid).tntype;
    
    switch tnt
        case 'pseudo-Anosov'
            fullmeasure(j,6) = 1;
        case 'finite-order'
            fullmeasure(j,6) = 2;
        case 'reducible'
            fullmeasure(j,6) = 3;
    end
    catch
        errors = errors+1;
    end
    waitbar(j/N, f, sprintf('Run: %d \t Progress: %d %% \t Errors: %d', j,floor(j/N*100),errors));
end
close(f);

% fig1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile;
%     histogram(fullmeasure(:,1),'Normalization','probability')
%     pbaspect([2 1 1])
%     xlabel('Writhe');
%     ylabel('Probability')
% nexttile;
%     histogram(fullmeasure(:,2),'Normalization','probability')
%     pbaspect([2 1 1])
%     xlabel('Complexity');
%     ylabel('Probability')
% nexttile;
%     histogram(fullmeasure(:,3),'Normalization','probability')
%     pbaspect([2 1 1])
%     xlabel('Length');
%     ylabel('Probability')
% exportgraphics(fig1,'4strandhistograms.pdf','BackgroundColor','none');

str2 = ['braid_measures_K' num2str(kappa,'%3.1f') '.txt'];
writematrix(fullmeasure,str2,'Delimiter',"\t");
close all;