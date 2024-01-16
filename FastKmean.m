function [centers,mincenter,mindist,q2,quality] = FastKmean(data,initcenters,method)
% This program compute the fast k-means as described in the paper
% "Using the Triangle Inequality to Accelerate k-Means"
% published in Proceedings of the Twentieth International Conference on Machine Learning (ICML'03).


if nargin < 3 method = 2; end
[n,dim] = size(data);

if max(size(initcenters)) == 1
    k = initcenters;
    [centers, mincenter, mindist, lower, computed] = anchors(mean(data),k,data);
    total = computed;
    skipestep = 1;
else
    centers = initcenters;
    mincenter = zeros(n,1);
    total = 0;
    skipestep = 0;
    [k,dim2] = size(centers);
    if dim ~= dim2 error('dim(data) ~= dim(centers)'); end;
end

nchanged = n;
iteration = 0;
oldmincenter = zeros(n,1);

while nchanged > 0
    % do one E step, then one M step
    computed = 0;
    
    if method == 0 & ~skipestep
        for i = 1:n
            for j = 1:k
                distmat(i,j) = calcdist(data(i,:),centers(j,:));
            end
        end
        [mindist,mincenter] = min(distmat,[],2);
        computed = k*n;
        
    elseif (method == 1 | (method == 2 & iteration == 0)) & ~skipestep
        mindist = Inf*ones(n,1);
        lower = zeros(n,k);
        for j = 1:k
            jdist = calcdist(data,centers(j,:));
            lower(:,j) = jdist;
            track = find(jdist < mindist);
            mindist(track) = jdist(track);
            mincenter(track) = j;
        end
        computed = k*n;
        
    elseif method == 2 & ~skipestep
        computed = 0;
        
        
        nndist = min(centdist,[],2);
        mobile = find(mindist > nndist(mincenter));
        
        mdm = mindist(mobile);
        mcm = mincenter(mobile);
        
        for j = 1:k
            track = find(mdm > centdist(mcm,j));
            if isempty(track) continue; end
            alt = find(mdm(track) > lower(mobile(track),j));
            if isempty(alt) continue; end
            track1 = mobile(track(alt));
            
            
            redo = find(~recalculated(track1));
            redo = track1(redo);
            c = mincenter(redo);
            computed = computed + size(redo,1);
            for jj = unique(c)'
                rp = redo(find(c == jj));
                udist = calcdist(data(rp,:),centers(jj,:));
                lower(rp,jj) = udist;
                mindist(rp) = udist;
            end
            recalculated(redo) = 1;
            
            track2 = find(mindist(track1) > centdist(mincenter(track1),j));
            track1 = track1(track2);
            if isempty(track1) continue; end
            
            % calculate exact distances to center j
            track4 = find(lower(track1,j) < mindist(track1));
            if isempty(track4) continue; end
            track5 = track1(track4);
            jdist = calcdist(data(track5,:),centers(j,:));
            computed = computed + size(track5,1);
            lower(track5,j) = jdist;
            
            % find which points really are assigned to center j
            track2 = find(jdist < mindist(track5));
            track3 = track5(track2);
            mindist(track3) = jdist(track2);
            mincenter(track3) = j;
        end % for j=1:k
    end % if method
    
    oldcenters = centers;
    
    
    
    diff = find(mincenter ~= oldmincenter);
    diffj = unique([mincenter(diff);oldmincenter(diff)])';
    diffj = diffj(find(diffj > 0));
    
    if size(diff,1) < n/3 & iteration > 0
        for j = diffj
            plus = find(mincenter(diff) == j);
            minus = find(oldmincenter(diff) == j);
            oldpop = pop(j);
            pop(j) = pop(j) + size(plus,1) - size(minus,1);
            if pop(j) == 0 continue; end
            centers(j,:) = (centers(j,:)*oldpop + sum(data(diff(plus),:),1) - sum(data(diff(minus),:),1))/pop(j);
        end
    else
        for j = diffj
            track = find(mincenter == j);
            pop(j) = size(track,1);
            if pop(j) == 0 continue; end
            % it's correct to have mean(data(track,:),1) but this can make answer worse!
            centers(j,:) = mean(data(track,:),1);
        end
    end
    
    if method == 2
        for j = diffj
            offset = calcdist(centers(j,:),oldcenters(j,:));
            computed = computed + 1;
            if offset == 0 continue; end
            track = find(mincenter == j);
            mindist(track) = mindist(track) + offset;
            lower(:,j) = max(lower(:,j) - offset,0);
        end
        
        % compute distance between each pair of centers
        % modify centdist to make "find" using it faster
        recalculated = zeros(n,1);
        realdist = alldist(centers);
        centdist = 0.5*realdist + diag(Inf*ones(k,1));
        computed = computed + k + k*(k-1)/2;
    end
    
    nchanged = size(diff,1) + skipestep;
    iteration = iteration+1;
    skipestep = 0;
    oldmincenter = mincenter;
    
    
    total = total + computed;
end % while nchanged > 0

udist = calcdist(data,centers(mincenter,:));
quality = mean(udist);
q2 = mean(udist.^2);
% [iteration toc quality q2 total]

