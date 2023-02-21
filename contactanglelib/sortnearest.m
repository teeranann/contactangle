function [locsorted]=sortnearest(trace)

nearestnum=2;

disall=sqrt((trace(:,1)-trace(:,1)').^2+(trace(:,2)-trace(:,2)').^2);
locsorted=zeros(size(trace,1),nearestnum^nearestnum);

j=1;
pathdist=zeros(1,size(trace,1));
for kk=1:size(trace,1)
    disallloop=disall;
    ii=1;
    for i=1:size(trace,1)
        if i==1
            a=0;
            b=kk;%size(trace,1); %loc starting point
        else
%             k=1;
%             bs=sort(disallloop(ii,:));
%             a=bs(k);
            [a,b]=min(disallloop(ii,:));
        end
        pathdist(kk)=pathdist(kk)+a;
        locsorted(i,j)=b;
        disallloop(:,b)=NaN;
        ii=b;
    end
end
[~,startloc]=min(pathdist);

ii=1;
looppp=unique(nchoosek(repmat(1:nearestnum, 1,4), nearestnum), 'rows');
pathdist=zeros(1,nearestnum^nearestnum);
for j=1:nearestnum^nearestnum
    disallloop=disall;
    for i=1:size(trace,1)
        if i==1
            a=0;
            b=startloc;%size(trace,1); %loc starting point
        else
            k=looppp(j,mod(i-2,nearestnum)+1);
            bs=sort(disallloop(ii,:));
            a=bs(k);
            [~,b]=(min(disallloop(ii,:)-bs(k)));
        end
        pathdist(j)=pathdist(j)+a;
        locsorted(i,j)=b;
        disallloop(:,b)=NaN;
        ii=b;
    end
end
[~,d]=min(pathdist);
locsorted=locsorted(:,d)';