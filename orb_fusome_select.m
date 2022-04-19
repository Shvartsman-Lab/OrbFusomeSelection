%orb_fusome_test3

clear; clc;
% Initialize fusome/ring system

x_voxel = 0.0696;
z_voxel = 0.2098;
voxelSize = [x_voxel x_voxel z_voxel];

%all samples should be 16 cells - get file for analysis
sample = '9_2';
smooth = 5; %smoothing factor for fusome/ring splitting (can adjust)
uiopen(strcat('C:\Users\rockyd\Desktop\Orb_Images\Fusome_Ring_Prob\Shot_YFP_Pav_RFP_',sample,'_Probabilities.h5'),1);
orb_file = strcat('C:\Users\rockyd\Desktop\Orb_Images\Orb Images\Orb_FISH_',sample,'.tif'); %also input file with orb data

ImageInfo = imfinfo(orb_file);
x_length = ImageInfo(1).Width;
y_length = ImageInfo(1).Height;
z_length = length(ImageInfo);

% preparing probabilities:
image_germ = double(I.img);
image_germ = permute(image_germ, [1 2 4 3]);
image_germ = image_germ/max(image_germ(:));
ring = image_germ(:,:,:,2);
fusome = image_germ(:,:,:,1);

% Clean up rings and segment fusome
bw_fusome = fusome > 0.5 | ring > 0.5;
frame = true(size(bw_fusome));
frame(2:end-1, 2:end-1, :) = false;
labels = bwlabeln(bw_fusome);
to_remove = integer_unique(labels(frame));
labels(ismember(labels, to_remove)) = 0;
labels(~isolate_lcc(labels)) = 0;
newfus = isolate_lcc(labels > 0 & fusome > 0.5);
newring = bwlabeln(labels > 0 & ring > 0.5);

%remove small objects that are not rings
s = bwskel(newring > 0);
change_occur = true;
while change_occur
    new_s = s;
    new_s(bwmorph3(s,'endpoints')) = false;
    change_occur = ~isequal(new_s, s);
    s = new_s;
end

% finish removing things that are not rings...
s = bwmorph3(s, 'clean');
rings_to_keep = integer_unique(nonzeros(newring(s)));
newring(~ismember(newring, rings_to_keep)) = 0;
newring = bwlabeln(newring > 0);
num_rings = length(unique(newring))-1;
fus_parts = length(unique(newfus))-1;
[x,y,z] = ind2sub(size(s), find(s));
x = x * voxelSize(1);
y = y * voxelSize(2);
z = z * voxelSize(3);

%expand out size to make isotropic directions
in_x = (0:size(newfus,1)-1)*voxelSize(1);
in_y = (0:size(newfus,2)-1)*voxelSize(2);
in_z = (0:size(newfus,3)-1)*voxelSize(3);

out_x = (0:size(newfus,1)-1)*voxelSize(1);
out_y = (0:size(newfus,2)-1)*voxelSize(1);
out_z = (0:size(newfus,3)-1)*voxelSize(1);

interp_newfus = medfilt3(interp3(single(newfus), 1:size(newfus,2), (1:size(newfus,1))', 1 : voxelSize(1)/voxelSize(3) : size(newfus,3), 'linear') >= 0.5);
interp_newring = medfilt3(bwlabeln(interp3(single(newring > 0), 1:size(newring,2), (1:size(newring,1))', 1 : voxelSize(1)/voxelSize(3) : size(newfus,3), 'linear') >= 0.5));
%% visualization before segmenting

for i = 1:num_rings
    p(i) = isosurface(interp_newring == i);
    p(i).vertices(:,1:2) = p(i).vertices(:,1:2)*voxelSize(1);
    p(i).vertices(:,3) = p(i).vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
end

figure;
for j = 1:length(p)
    patch(p(j),'FaceColor','k','EdgeColor','none'); axis image; alpha(0.6)
end

p1 = isosurface(interp_newfus == 1);
p1.vertices(:,1:2) = p1.vertices(:,1:2)*voxelSize(1);
p1.vertices(:,3) = p1.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p1,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)
box on; grid on;
view([100 15])
%% split fusome and rings into component parts
clc;
newfus_labeled = split_fusome(0, interp_newfus, interp_newring, smooth);
%% Find adjacencies between objects
bw = imdilate(interp_newring, create_ball(1));
newimage = newfus_labeled.*(bw > 0);

if exist('list','var')
    clear list
end

for i = 1:(num_rings)
    check_obj = bwconncomp(bw == i,26);
    props = regionprops(check_obj,'Area','BoundingBox','Centroid','PixelIdxList','PixelList');
    list{i} = props.PixelIdxList;
end

matchup = zeros(length(list),2);

for j = 1:length(list)
    output = newimage(list{j});

    number = zeros(1,num_rings);
    for k = 1:num_rings+1
        number(k) = sum(output == k);
    end

    if nnz(number) < 2
        matchup(j,1) = 0;
        matchup(j,2) = 0;
    elseif nnz(find(number == max(number))) == 2
        [indices] = find(number == max(number));
        matchup(j,1) = indices(1);
        matchup(j,2) = indices(2);
    elseif max(number) == max(number(number < max(number)))
        matchup(j,1) = find(number == max(number));
        matchup(j,2) = find(number == max(number));
    else
        matchup(j,1) = find(number == max(number));
        matchup(j,2) = find(number == max(number(number<max(number))));
    end
end

% Remove zero rows and rows with one element
rule = matchup(:,2) == 0;
matchup(rule,:) = [];

% Remove zero columns
matchup(:,all(~matchup,1)) = [];

matchupswitch = [matchup(:,2) matchup(:,1)];
matchup = unique([matchup; matchupswitch], 'rows');

%Remove duplicate pairs
pair = matchup(:,1) < matchup(:,2);
matchup = matchup.*pair;
rule = matchup(:,2) == 0;
matchup(rule,:) = [];

total_seg = max(matchup(:));

%get fusome volumes for each portion
fus_vol_out = zeros(total_seg,1);
for i = 1:total_seg
    fus_vol_out(i) = nnz(newfus_labeled(:) == i).*voxelSize(1).*voxelSize(2).*voxelSize(3).*voxelSize(1)./voxelSize(3);
end
oocyte = find(fus_vol_out == max(fus_vol_out));

matchupmatrix = zeros(total_seg);
for i = 1:size(matchup,1)
    matchupmatrix(matchup(i,1),matchup(i,2)) = 1;
    matchupmatrix(matchup(i,2),matchup(i,1)) = 1;
end

%compare adjacency matrix with known 16-cell adjacency matrix
truematchup = [1 2;1 3;1 5;1 9;2 4; 2 6;2 10;3 7;3 11;4 12;4 8;5 13;6 14;7 15;8 16];
truematchupmatrix = zeros(16);

for i = 1:size(truematchup,1)
    truematchupmatrix(truematchup(i,1),truematchup(i,2)) = 1;
    truematchupmatrix(truematchup(i,2),truematchup(i,1)) = 1;
end

matchupmatrix = zeros(total_seg);
for i = 1:size(matchup,1)
    matchupmatrix(matchup(i,1),matchup(i,2)) = 1;
    matchupmatrix(matchup(i,2),matchup(i,1)) = 1;
end

matchupmatrix(oocyte,:) = 2*matchupmatrix(oocyte,:);
matchupmatrix(:,oocyte) = 2*matchupmatrix(:,oocyte);

truematchupmatrix(1,:) = 2*truematchupmatrix(1,:);
truematchupmatrix(:,1) = 2*truematchupmatrix(:,1);

[P1, ~] = eig(truematchupmatrix);
[P2, ~] = eig(matchupmatrix);

P1 = abs(P1);
P2 = abs(P2);

outputmatrix = P2*P1';
[permutation, ~] = munkres(outputmatrix);
permutation = ((total_seg+1) - permutation)';

% Relabel fusome pieces
output_fus = zeros(size(newfus_labeled));
fusome_vol = zeros(size(fus_vol_out));
for m = 1:(num_rings+1)
    output_fus(newfus_labeled == m) = permutation(m);
    fusome_vol(permutation(m)) = fus_vol_out(m);
end
%%
% Create 3D reconstruction

for i = 1:num_rings
    p(i) = isosurface(interp_newring == i);
    p(i).vertices(:,1:2) = p(i).vertices(:,1:2)*voxelSize(1);
    p(i).vertices(:,3) = p(i).vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
end

figure;
for j = 1:length(p)
    patch(p(j),'FaceColor','k','EdgeColor','none'); axis image; alpha(0.6)
end

p1 = isosurface(output_fus == 1);
p1.vertices(:,1:2) = p1.vertices(:,1:2)*voxelSize(1);
p1.vertices(:,3) = p1.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p1,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none'); axis image; alpha(0.6)

p2 = isosurface(output_fus == 2);
p2.vertices(:,1:2) = p2.vertices(:,1:2)*voxelSize(1);
p2.vertices(:,3) = p2.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p2,'FaceColor','b','EdgeColor','none'); axis image; alpha(0.6)

p3 = isosurface(output_fus == 3);
p3.vertices(:,1:2) = p3.vertices(:,1:2)*voxelSize(1);
p3.vertices(:,3) = p3.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p3,'FaceColor','b','EdgeColor','none'); axis image; alpha(0.6)

p4 = isosurface(output_fus == 4);
p4.vertices(:,1:2) = p4.vertices(:,1:2)*voxelSize(1);
p4.vertices(:,3) = p4.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p4,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)

p5 = isosurface(output_fus == 5);
p5.vertices(:,1:2) = p5.vertices(:,1:2)*voxelSize(1);
p5.vertices(:,3) = p5.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p5,'FaceColor','b','EdgeColor','none'); axis image; alpha(0.6)

p6 = isosurface(output_fus == 6);
p6.vertices(:,1:2) = p6.vertices(:,1:2)*voxelSize(1);
p6.vertices(:,3) = p6.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p6,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)

p7 = isosurface(output_fus == 7);
p7.vertices(:,1:2) = p7.vertices(:,1:2)*voxelSize(1);
p7.vertices(:,3) = p7.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p7,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)

p8 = isosurface(output_fus == 8);
p8.vertices(:,1:2) = p8.vertices(:,1:2)*voxelSize(1);
p8.vertices(:,3) = p8.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p8,'FaceColor','g','EdgeColor','none'); axis image; alpha(0.6)

p9 = isosurface(output_fus == 9);
p9.vertices(:,1:2) = p9.vertices(:,1:2)*voxelSize(1);
p9.vertices(:,3) = p9.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p9,'FaceColor','b','EdgeColor','none'); axis image; alpha(0.6)

p10 = isosurface(output_fus == 10);
p10.vertices(:,1:2) = p10.vertices(:,1:2)*voxelSize(1);
p10.vertices(:,3) = p10.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p10,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)

p11 = isosurface(output_fus == 11);
p11.vertices(:,1:2) = p11.vertices(:,1:2)*voxelSize(1);
p11.vertices(:,3) = p11.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p11,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)

p12 = isosurface(output_fus == 12);
p12.vertices(:,1:2) = p12.vertices(:,1:2)*voxelSize(1);
p12.vertices(:,3) = p12.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p12,'FaceColor','g','EdgeColor','none'); axis image; alpha(0.6)

p13 = isosurface(output_fus == 13);
p13.vertices(:,1:2) = p13.vertices(:,1:2)*voxelSize(1);
p13.vertices(:,3) = p13.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p13,'FaceColor','r','EdgeColor','none'); axis image; alpha(0.6)

p14 = isosurface(output_fus == 14);
p14.vertices(:,1:2) = p14.vertices(:,1:2)*voxelSize(1);
p14.vertices(:,3) = p14.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p14,'FaceColor','g','EdgeColor','none'); axis image; alpha(0.6)

p15 = isosurface(output_fus == 15);
p15.vertices(:,1:2) = p15.vertices(:,1:2)*voxelSize(1);
p15.vertices(:,3) = p15.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p15,'FaceColor','g','EdgeColor','none'); axis image; alpha(0.6)

p16 = isosurface(output_fus == 16);
p16.vertices(:,1:2) = p16.vertices(:,1:2)*voxelSize(1);
p16.vertices(:,3) = p16.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p16,'FaceColor','y','EdgeColor','none'); axis image; alpha(0.6)

box on; grid on; axis off; axis image
view([100 15])
%%
% quantifying orb
orbThresh = 0.2;
orb_image=zeros(y_length,x_length,z_length,'double');
for i=1:z_length
    orb_image(:,:,i)=imread(orb_file,'Index',i);
end
%swap axes of image compared with MATLAB output
orb_image = permute(orb_image,[2 1 3]);
orb_image = orb_image./max(orb_image(:));
orb_image = orb_image > orbThresh; %find where orb is above threshold level

p = isosurface(orb_image);
p.vertices(:,1:2) = p.vertices(:,1:2)*voxelSize(1);
p.vertices(:,3) = p.vertices(:,3)*voxelSize(3);
patch(p,'FaceColor','c','EdgeColor','none'); axis image; alpha(0.6)

orb_image_interp = medfilt3(interp3(single(orb_image), 1:size(orb_image,2), (1:size(orb_image,1))', 1 : voxelSize(1)/voxelSize(3) : size(orb_image,3), 'nearest') >= orbThresh);
p = isosurface(orb_image_interp);
p.vertices = p.vertices*voxelSize(1);
patch(p,'FaceColor','c','EdgeColor','none'); axis image; alpha(0.6)

box on; grid on;
view([100 15])

%quantify where orb and fusome portion in two central cells meet
vol1 = nnz(orb_image_interp.*(output_fus == 1))*voxelSize(1).*voxelSize(1).*voxelSize(1);
vol2 = nnz(orb_image_interp.*(output_fus == 2))*voxelSize(1).*voxelSize(1).*voxelSize(1);

%get volume and orb ratio of two central cells
vol12_ratio_orb = vol1./vol2;
vol12_ratio_fus = fusome_vol(1)./fusome_vol(2);
