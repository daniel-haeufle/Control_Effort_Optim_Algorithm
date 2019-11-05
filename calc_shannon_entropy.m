function H = calc_shannon_entropy(x)

if size(x,1)<size(x,2)
    x=x';
end

x = double(x);

[x_unique,~,ic] = unique(x,'rows');
for i=1:size(x_unique,1)
    p_x(i,:)=[ x_unique(i,:),length(find(ic==i))/size(x,1) ]; 
end

H = 0;
for i = 1:size(p_x,1)
        H = H - p_x(i,2) * log2(p_x(i,2));
end