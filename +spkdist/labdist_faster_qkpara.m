function labdist_faster_qk_mat = labdist_faster_qkpara(sa,la,sb,lb,q,k)
% LABDIST_FASTER_QKPARA(SA,LA,SB,LB,Q,K), parallel in q and k
% Calculates the multi-unit metric distance between two spike trains
% Uses a fast version of the algorithm 
% SA, SB - spike times on the two spike trains
% LA, LB - spike labels (positive integers) 
% Q - vector of timing precision parameters 
% K - vector of label reassigning parameters 
%
% Thomas Kreuz, 10/19/08; based on code by Dmitriy Aronov, 6/20/01;

%Assign labels in the form 1,2,...,L and count spikes of each label 
lbs = unique([la lb]); 
L = size(lbs,2);
for c = 1:L 
    j = find(la==lbs(c));
    la(j) = c; 
    numa(c) = size(j,2); 
    j = find(lb==lbs(c)); 
    lb(j) = c; 
    numb(c) = size(j,2); 
end
%Choose the spike train to separate to subtrains 
if prod(numb+1)*(sum(numa)+1) > prod(numa+1)*(sum(numb)+1)
    t = la; 
    la = lb; 
    lb = t; 
    t = sa; 
    sa = sb; 
    sb = t; 
    t = numa; 
    % numa=numb;
    numb = t;
end
tb=zeros(L,max(numb));
for c = 1:L 
    tb(c,1:numb(c))=sb(logical(lb==c));
end
%Set up an indexing system 
ind = []; 
for c = 1:L 
    j = repmat(0:numb(c),prod(numb(c+1:end)+1),1);
    j = repmat(reshape(j,numel(j),1),prod(numb(1:c-1)+1),1);
    ind = [ind j];
end
ind = sortrows([sum(ind,2) ind]); 
ind = ind(:,2:end);
%Initialize the array
m = zeros(size(ind,1),size(sa,2)+1); 
m(1,:) = 0:size(sa,2); 
m(:,1) = sum(ind,2); 
m = repmat(shiftdim(m,-2),[length(k),length(q),1,1]);
%Perform the calculation
for v = 2:size(m,3)
    fa2=find(shiftdim(m(1,1,:,1),2)==m(1,1,v,1)-1);
    fa=fa2(logical(sum(ind(fa2,:)-repmat(ind(v,:),length(fa2),1)==0,2)==L-1));
    fth=find(ind(v,:)>0)';
    bsv=diag(tb(fth,ind(v,fth)));
    for w = 2:size(m,4)       
        m(:,:,v,w)=min(cat(3,m(:,:,v,w-1)+1,m(:,:,fa,w)+1,m(:,:,fa,w-1)+ ...
            repmat(shiftdim(q'*abs(sa(w-1)-bsv'),-1),[length(k),1,1])+ ...
            permute(repmat(shiftdim(k'*not(la(w-1)==fth)',-1),[length(q),1,1]),[2 1 3])),[],3);
    end
end
labdist_faster_qk_mat = m(:,:,end,end);