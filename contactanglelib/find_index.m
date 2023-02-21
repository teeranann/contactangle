function index=find_index(trace,p1,p2)
% Small function that finds the coordinate in a list (trace) closest to the line
% between points p1 and p2. This is done by calculating the distance between
% all points in the list and the line
% All coordinates are in x,y

a=trace;
b=bsxfun(@plus,zeros(size(a)),p1);
c=bsxfun(@plus,zeros(size(a)),p2);
dist=((a(:,2)-b(:,2)).*(c(:,1)-b(:,1))-(a(:,1)-b(:,1)).*(c(:,2)-b(:,2)))/norm(p2-p1);
[~,index]=min(abs(dist));