function [fusome_labeled, last_used_index] = split_fusome(last_used_index, fusome, all_rings, smooth) %split fusome into connected components within cells

% if there's no ring - update ind to ind+1 and label all voxels with ind, return the labeled fusome and ind
% else: pick a ring, split the fusome into two parts (this is the recursion)
ring_idcs = nonzeros(integer_unique(all_rings));
if isempty(ring_idcs) % STOP CONDITION MET!!
    last_used_index = last_used_index + 1;
    fusome_labeled = fusome .* last_used_index;
else
    [fus1, fus2] = split_fusome_at_ring(fusome, all_rings == ring_idcs(1), smooth);
    [ring1, ring2] = split_rings_at_ring(fus1, fus2, all_rings-ring_idcs(1).*(all_rings == ring_idcs(1)), smooth);
    [fus1_labeled, last_used_index] = split_fusome(last_used_index, fus1, ring1, smooth);
    [fus2_labeled, last_used_index] = split_fusome(last_used_index, fus2, ring2, smooth);
    fusome_labeled = fus1_labeled + fus2_labeled;
end
end

function [fus1, fus2] = split_fusome_at_ring(fusome, splitting_ring, smooth) %use the ring to split fusome into 2 segments

%make convex hull and use as way to remove these pixels from fusome
out = regionprops3(splitting_ring, "all");
xyz = out.ConvexImage{1,1};
[x,y,z] = ind2sub(size(xyz),find(xyz));
box = out.BoundingBox(1,:);

ring_hull = zeros(size(splitting_ring));

for i = 1:length(x)
    ring_hull(x(i)+box(2)-0.5,y(i)+box(1)-0.5,z(i)+box(3)-0.5) = 1;
end

%close and smooth the convex hull of the ring for further divisions
ring_hull = imclose(ring_hull,ones([30 30 30]));
ring_hull = imdilate(ring_hull,ones([smooth smooth smooth]));

[rxyz(:,1),rxyz(:,2),rxyz(:,3)] = ind2sub(size(ring_hull),find(ring_hull));

fusome_copy = fusome;
for j = 1:length(rxyz)
    fusome_copy(rxyz(j,1),rxyz(j,2),rxyz(j,3)) = 0;
end

%find leftover pieces of fusome that were removed
[fxyz(:,1),fxyz(:,2),fxyz(:,3)] = ind2sub(size(fusome-fusome_copy),find(fusome-fusome_copy));

%split fusome into 2 parts and see if you only end with 2 pieces
obj = bwconncomp(fusome_copy, 26);
[~, order] = sort(cellfun(@length, obj.PixelIdxList), 'descend');
obj.PixelIdxList = obj.PixelIdxList(order);
split_fus = regionprops3(obj,"VoxelIdxList","VoxelList", "Centroid", "BoundingBox");
fus1_pts = split_fus.VoxelList{1};
fus2_pts = split_fus.VoxelList{2};

%test one way to match the pieces together
fus1 = zeros(size(fusome));
for i = 1:length(fus1_pts)
    fus1(fus1_pts(i,2),fus1_pts(i,1),fus1_pts(i,3)) = 1;
end

fus2 = zeros(size(fusome));
for i = 1:length(fus2_pts)
    fus2(fus2_pts(i,2),fus2_pts(i,1),fus2_pts(i,3)) = 1;
end


dist1 = (fusome-fusome_copy).*bwdist(fus1,'euclidean');
dist2 = (fusome-fusome_copy).*bwdist(fus2,'euclidean');

for i = 1:length(fxyz)
    if dist1(fxyz(i,1),fxyz(i,2),fxyz(i,3)) <= dist2(fxyz(i,1),fxyz(i,2),fxyz(i,3))
        fus1(fxyz(i,1),fxyz(i,2),fxyz(i,3)) = 1;
    elseif dist1(fxyz(i,1),fxyz(i,2),fxyz(i,3)) > dist2(fxyz(i,1),fxyz(i,2),fxyz(i,3))
        fus2(fxyz(i,1),fxyz(i,2),fxyz(i,3)) = 1;
    end
end


%output is two binary fusome portions
fus1 = logical(fus1);
fus2 = logical(fus2);
end

function [ring1, ring2] = split_rings_at_ring(fus1, fus2, rings_left, smooth) %divide up the remaining rings for use in the recursion

fus1_increase = imdilate(fus1,ones(2,2,2));
fus1_rings = integer_unique(nonzeros(fus1_increase.*rings_left));


fus2_increase = imdilate(fus2,ones(2,2,2));
fus2_rings = integer_unique(nonzeros(fus2_increase.*rings_left));


%scalar mult retains the rings along with their identity for further segmentation
ring1 = zeros(size(rings_left));
for i = 1:length(fus1_rings)
ring1(rings_left == fus1_rings(i)) = fus1_rings(i);
end

ring2 = zeros(size(rings_left));
for i = 1:length(fus2_rings)
ring2(rings_left == fus2_rings(i)) = fus2_rings(i);
end
end








