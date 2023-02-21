close all
clearvars -except pathnamm

% 21/3/23 
% T. Nongnual, K. Saelee, T. Pakminakhom
% Burapha University

% The script can be divided into several steps
% 1. Loading all data (video, tilt values, callibration functions)
% 2. Edge detection and baseline determination in each frame
% 3. Shift baseline coordinates from each frame into coordinate system of
%       first frame. Then all baseline coordinates are fitted into a single
%       baseline used for all frames
% 4. Calculate contact angles for all frames using the determined baseline.
%       In this part of the script there is the possibility to plot each
%       frame with the corresponding contact angles and fit.
% 5. Plot the fitted contact angles and tripple line displacements as a
%       function of tilt.

% Additional feature (21 Feb 2023)
% 1. A checkbox used to choose between the method to find the base line and
%       the option to crop out only the drop.
% 2. The checkbox also come with the option to trim only the 1st frame to calculate.
% 3. The method to find the baseline using reflection.
% 4. The method to find the baseline using our device.
% 5. The calculated contact angle also come with the built in delete
%       outlier.

polycal=0;
ellipcal=1;
zpixelcutoff=10 ;
caadsorpcutoff=10 ;
contactcalsec=1;

myOutput=checkboxtest;
Qrotstage=myOutput(1);
Qcrop=myOutput(3);
Qmanual=myOutput(4);
Qreflection=myOutput(5);
Qneedle=myOutput(6);

warning off MATLAB:subscripting:noSubscriptsSpecified
cc = lines(25);

%%
% %---------------------------------------------------------------------------
% %                1.+2.       Loading of data + Edge detection
% %---------------------------------------------------------------------------

edgedet

%%
%---------------------------------------------------------------------------
% %        3.       Average baselines found in all frames.
% %-------------------------------------------------------------------------

%---------------------------------------------------------------------------
% %        3.2       Crop
% %-------------------------------------------------------------------------

if Qcrop
    figure(1)
    imagesc(mov{1})
    axis ij
    axis equal
    hold on
    colormap gray
    xlim([1 camwidth]);
    ylim([1 camheight]);
    xlim([0 camwidth]);
    ylim([0 camheight]);
    set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.02,0.02])
    set(gcf,'windowstyle','normal')
    xlabel('x (pixel)')
    ylabel('y (pixel)')

    disp('Please draw a rectangle around the drop form top left to buttom right then right click and choose crop.')
    [xrefout,yrefout]=imcrop;
    xi=yrefout(1,1);
    xf=yrefout(1,1)+yrefout(1,3);
    yi=yrefout(1,2);
    yf=yrefout(1,2)+yrefout(1,4);

    sumcrop=(sum(mov{1,1}(round(yi):round(yf),round(xi):round(xf)),2))/255;
    if mean(sumcrop(end-2:end))<50
        [~,yfcutoffdark]=findpeaks(abs(gradient(sumcrop(round(size(sumcrop,1)*1/2):end))),'SortStr', 'descend',"NPeaks",1);
        yf=yf-(round(size(sumcrop,1)*1/2)-yfcutoffdark-0);
        disp('Dark sample surface in crop area is detected. Removed.')
        disp(' ')
    end

    if logical(Qcrop*Qmanual)
        figure(10)
        imagesc(mov{1})
        axis ij
        axis equal
        hold on
        colormap gray
        xlim([xi xf]);
        ylim([yi yf+30]);
        % xlim([0 camwidth]);
        % ylim([0 camheight]);
        set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.02,0.02])
        set(gcf,'windowstyle','normal')
        xlabel('x (pixel)')
        ylabel('y (pixel)')
    else
        figure(1);plot([xi xi xf xf xi], [yf yi yi yf yf],'r') % plot crop box of the drop
    end


    EdgeDrop=cell(camframetot,2);
    EdgeDropTotal=cell(camframetot,2);
    for i=1:camframetot
        EdgeDrop{i,1}.x=[];EdgeDrop{i,1}.y=[];EdgeDrop{i,1}.position=[];
        for j=1:size(EdgeCellSALR{i}.x,1)
            if EdgeCellSALR{i}.y(j)>yi && EdgeCellSALR{i}.y(j)<yf && EdgeCellSALR{i}.x(j)>xi && EdgeCellSALR{i}.x(j)<xf
                EdgeDrop{i,1}.x(end+1,1)=EdgeCellSALR{i}.x(j);
                EdgeDrop{i,1}.y(end+1,1)=EdgeCellSALR{i}.y(j);
                EdgeDrop{i,1}.position(end+1,1)=EdgeCellSALR{i}.position(j);
            end
        end

        EdgeDrop{i,1}=findlongestedge(EdgeDrop{i,1},[camheight camwidth],marginedgediff);
        if i==1
            eLRDataD=EllipticFitVolnotilt(EdgeDrop{i,1},EdgeDrop{i,1},[xi,yf],[xf,yf],zpixelcutoff+50);
            crix=eLRDataD.ellipse{1,1}.X0_in;
            criy=eLRDataD.ellipse{1,1}.Y0_in-eLRDataD.ellipse{1,1}.b;
        end
        locdropL=find(EdgeDrop{i,1}.x<=crix);
        locdropR=find(EdgeDrop{i,1}.x>crix);
        EdgeDroptemp=EdgeDrop;
        EdgeDrop{i,1}.x=EdgeDroptemp{i,1}.x(locdropL,1);    %left
        EdgeDrop{i,2}.x=EdgeDroptemp{i,1}.x(locdropR,1);   %right
        EdgeDrop{i,1}.y=EdgeDroptemp{i,1}.y(locdropL,1);    %left
        EdgeDrop{i,2}.y=EdgeDroptemp{i,1}.y(locdropR,1);   %right
        EdgeDrop{i,1}.position=EdgeDroptemp{i,1}.position(locdropL,1);   %left
        EdgeDrop{i,2}.position=EdgeDroptemp{i,1}.position(locdropR,1);   %right

        locsortL=sortnearest([EdgeDrop{i,1}.x,EdgeDrop{i,1}.y]);
        EdgeDropLR{i,1}.x=EdgeDrop{i,1}.x(locsortL);
        EdgeDropLR{i,1}.y=EdgeDrop{i,1}.y(locsortL);
        EdgeDropLR{i,1}.position=EdgeDrop{i,1}.position(locsortL);

        locsortR=sortnearest([EdgeDrop{i,2}.x,EdgeDrop{i,2}.y]);
        EdgeDropLR{i,2}.x=EdgeDrop{i,2}.x(locsortR);
        EdgeDropLR{i,2}.y=EdgeDrop{i,2}.y(locsortR);
        EdgeDropLR{i,2}.position=EdgeDrop{i,2}.position(locsortR);

        EdgeDropTotal{i,1}.x=[flip(EdgeDropLR{i,1}.x);EdgeDropLR{i,2}.x];
        EdgeDropTotal{i,1}.y=[flip(EdgeDropLR{i,1}.y);EdgeDropLR{i,2}.y];
        EdgeDropTotal{i,1}.position=[flip(EdgeDropLR{i,1}.position);EdgeDropLR{i,2}.position];
    end

    EdgeCellLR=EdgeDropLR;
    EdgeCellTotal=EdgeDropTotal;
end

if Qneedle
    % EDGE NEEDLE
    EdgeNeedleL=cell(camframetot,2);
    EdgeNeedleLL=cell(camframetot,2);
    EdgeNeedleLR=cell(camframetot,2);
    EdgeNeedleR=cell(camframetot,2);
    EdgeNeedleRL=cell(camframetot,2);
    EdgeNeedleRR=cell(camframetot,2);
    % plot crop box of the needle (L/R)
    needlecutoff=zpixelcutoff+35;
    figure(1);plot([0 0 xi xi 0], [yf-needlecutoff yi yi yf-needlecutoff yf-needlecutoff],'b');
    figure(1);plot([xf xf camwidth camwidth xf], [yf-needlecutoff yi yi yf-needlecutoff yf-needlecutoff],'b');

    for i=1:camframetot
        EdgeNeedleL{i,1}.x=[];EdgeNeedleL{i,1}.y=[];EdgeNeedleL{i,1}.position=[];
        for j=1:size(EdgeCellSALR{i}.x,1)
            if EdgeCellSALR{i}.y(j)>yi && EdgeCellSALR{i}.y(j)<yf-needlecutoff && EdgeCellSALR{i}.x(j)<xi
                EdgeNeedleL{i,1}.x(end+1,1)=EdgeCellSALR{i}.x(j);
                EdgeNeedleL{i,1}.y(end+1,1)=EdgeCellSALR{i}.y(j);
                EdgeNeedleL{i,1}.position(end+1,1)=EdgeCellSALR{i}.position(j);
            end
        end

        EdgeNeedleR{i,1}.x=[];EdgeNeedleR{i,1}.y=[];EdgeNeedleR{i,1}.position=[];
        for j=1:size(EdgeCellSALR{i}.x,1)
            if EdgeCellSALR{i}.y(j)>yi && EdgeCellSALR{i}.y(j)<yf-needlecutoff && EdgeCellSALR{i}.x(j)>xf
                EdgeNeedleR{i,1}.x(end+1,1)=EdgeCellSALR{i}.x(j);
                EdgeNeedleR{i,1}.y(end+1,1)=EdgeCellSALR{i}.y(j);
                EdgeNeedleR{i,1}.position(end+1,1)=EdgeCellSALR{i}.position(j);
            end
        end
    end


    % Determine a middle point of a needle and create the linear equation to splitthe needle into 4 linear equation
    SumNeedleLlow=sum(mov{1,1}(round(yf-needlecutoff):round(yf-needlecutoff+2),1:round(xi)),1);
    [pSumNeedleLlow,locsLl]=findpeaks(abs(gradient(SumNeedleLlow)),"MinPeakDistance",5,"SortStr","descend","NPeaks",2);
    coLl=[(locsLl(1,1)+locsLl(1,2))/2,yf-needlecutoff];

    SumNeedleLhigh=sum(mov{1,1}(round(yi):round(yi+2),1:round(xi)),1);
    [pSumNeedleLhigh,locsLh]=findpeaks(abs(gradient(SumNeedleLhigh)),"MinPeakDistance",5,"SortStr","descend","NPeaks",2);
    coLh=[(locsLh(1,1)+locsLh(1,2))/2,yi];

    slopeL=(coLl(1,2)-coLh(1,2))/(coLl(1,1)-coLh(1,1));

    SumNeedleRlow=sum(mov{1,1}(round(yf-needlecutoff):round(yf-needlecutoff+2),round(xf):camwidth),1);
    [pSumNeedleRlow,locsRl]=findpeaks(abs(gradient(SumNeedleRlow)),"MinPeakDistance",5,"SortStr","descend","NPeaks",2);
    coRl=[((locsRl(1,1)+xf)+(locsRl(1,2)+xf))/2,yf-needlecutoff];

    SumNeedleRhigh=sum(mov{1,1}(round(yi):round(yi+2),round(xf):camwidth),1);
    [pSumNeedleRhigh,locsRh]=findpeaks(abs(gradient(SumNeedleRhigh)),"MinPeakDistance",5,"SortStr","descend","NPeaks",2);
    coRh=[((locsRh(1,1)+xf)+(locsRh(1,2)+xf))/2,yi];

    slopeR=(coRl(1,2)-coRh(1,2))/(coRl(1,1)-coRh(1,1));

    for i=1:camframetot
        EdgeNeedleLL{i,1}.x=[];EdgeNeedleLL{i,1}.y=[];EdgeNeedleLL{i,1}.position=[];
        EdgeNeedleLR{i,1}.x=[];EdgeNeedleLR{i,1}.y=[];EdgeNeedleLR{i,1}.position=[];

        EdgeNeedleLL{i,1}.x=EdgeNeedleL{i}.x(EdgeNeedleL{i}.x<((EdgeNeedleL{i,1}.y-coLh(1,2))/slopeL)+coLh(1,1));
        EdgeNeedleLL{i,1}.y=EdgeNeedleL{i}.y(EdgeNeedleL{i}.x<((EdgeNeedleL{i,1}.y-coLh(1,2))/slopeL)+coLh(1,1));
        EdgeNeedleLL{i,1}.position=EdgeNeedleL{i}.position(EdgeNeedleL{i}.x<((EdgeNeedleL{i,1}.y-coLh(1,2))/slopeL)+coLh(1,1));
        EdgeNeedleLL{i,1}=findlongestedge(EdgeNeedleLL{i,1},[camheight camwidth],3);

        EdgeNeedleLR{i,1}.x=EdgeNeedleL{i}.x(EdgeNeedleL{i}.x>((EdgeNeedleL{i,1}.y-coLh(1,2))/slopeL)+coLh(1,1));
        EdgeNeedleLR{i,1}.y=EdgeNeedleL{i}.y(EdgeNeedleL{i}.x>((EdgeNeedleL{i,1}.y-coLh(1,2))/slopeL)+coLh(1,1));
        EdgeNeedleLR{i,1}.position=EdgeNeedleL{i}.position(EdgeNeedleL{i}.x>((EdgeNeedleL{i,1}.y-coLh(1,2))/slopeL)+coLh(1,1));
        EdgeNeedleLR{i,1}=findlongestedge(EdgeNeedleLR{i,1},[camheight camwidth],3);


        EdgeNeedleRL{i,1}.x=[];EdgeNeedleRL{i,1}.y=[];EdgeNeedleRL{i,1}.position=[];
        EdgeNeedleRR{i,1}.x=[];EdgeNeedleRR{i,1}.y=[];EdgeNeedleRR{i,1}.position=[];

        EdgeNeedleRL{i,1}.x=EdgeNeedleR{i}.x(EdgeNeedleR{i}.x<((EdgeNeedleR{i,1}.y-coRh(1,2))/slopeR)+coRh(1,1));
        EdgeNeedleRL{i,1}.y=EdgeNeedleR{i}.y(EdgeNeedleR{i}.x<((EdgeNeedleR{i,1}.y-coRh(1,2))/slopeR)+coRh(1,1));
        EdgeNeedleRL{i,1}.position=EdgeNeedleR{i}.position(EdgeNeedleR{i}.x<((EdgeNeedleR{i,1}.y-coRh(1,2))/slopeR)+coRh(1,1));
        EdgeNeedleRL{i,1}=findlongestedge(EdgeNeedleRL{i,1},[camheight camwidth],3);

        EdgeNeedleRR{i,1}.x=EdgeNeedleR{i}.x(EdgeNeedleR{i}.x>((EdgeNeedleR{i,1}.y-coRh(1,2))/slopeR)+coRh(1,1));
        EdgeNeedleRR{i,1}.y=EdgeNeedleR{i}.y(EdgeNeedleR{i}.x>((EdgeNeedleR{i,1}.y-coRh(1,2))/slopeR)+coRh(1,1));
        EdgeNeedleRR{i,1}.position=EdgeNeedleR{i}.position(EdgeNeedleR{i}.x>((EdgeNeedleR{i,1}.y-coRh(1,2))/slopeR)+coRh(1,1));
        EdgeNeedleRR{i,1}=findlongestedge(EdgeNeedleRR{i,1},[camheight camwidth],3);
    end

    % Find the coordinate of the baseline
    for i=1:camframetot
        [mcLL,~]=polyfit(EdgeNeedleLL{i,1}.x,EdgeNeedleLL{i,1}.y,1);
        [mcLR,~]=polyfit(EdgeNeedleLR{i,1}.x,EdgeNeedleLR{i,1}.y,1);
        Lx(i)=((mcLR(1,2)-mcLL(1,2))/(mcLL(1,1)-mcLR(1,1)));
        Ly(i)=((mcLR(1,1)*Lx(i))+mcLR(1,2));
        Lxy(i,1:2)=[Lx(i),Ly(i)];

        [mcRL,~]=polyfit(EdgeNeedleRL{i,1}.x,EdgeNeedleRL{i,1}.y,1);
        [mcRR,~]=polyfit(EdgeNeedleRR{i,1}.x,EdgeNeedleRR{i,1}.y,1);
        Rx(i)=((mcRL(1,2)-mcRR(1,2))/(mcRR(1,1)-mcRL(1,1)));
        Ry(i)=((mcRL(1,1)*Rx(i))+mcRL(1,2));
        Rxy(i,1:2)=[Rx(i),Ry(i)];
    end
end

%%
%---------------------------------------------------------------------------
% %        3.3       Manually set baseline
% %---------------------------------------------------------------------------
if Qrotstage
    camframetot=1;
end
if logical(Qcrop*Qmanual)
    figure(10)
    plotint=round(camframetot/10);
    for n_frame=1:camframetot
        if mod(n_frame,plotint)==0 || n_frame==1 || n_frame==camframetot
            scatter(EdgeCellLR{n_frame,1}.x,EdgeCellLR{n_frame,1}.y,1,'m')
            scatter(EdgeCellLR{n_frame,2}.x,EdgeCellLR{n_frame,2}.y,1,'c')
        end
    end
else
    figure(2); imagesc(mov{1})
    axis ij; axis equal
    hold on
    colormap gray
    xlim([1 camwidth]);
    ylim([1 camheight]);
    % plot([xi xi xf xf xi], [yf yi yi yf yf],'r') % plot crop box of the drop
    set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.02,0.02])
    set(gcf,'windowstyle','normal')
    xlabel('x (pixel)')
    ylabel('y (pixel)')
    plotint=round(camframetot/10);
    for n_frame=1:camframetot
        if mod(n_frame,plotint)==0 || n_frame==1 || n_frame==camframetot
            scatter(EdgeCellLR{n_frame,1}.x,EdgeCellLR{n_frame,1}.y,1,'m')
            scatter(EdgeCellLR{n_frame,2}.x,EdgeCellLR{n_frame,2}.y,1,'c')
        end
    end
end

AverageBaseVecY=zeros(camframetot,4);
if Qmanual
    disp('Please manually pick two baseline points (two left clicks, then one right click)')
    [x0L,y0L,~]=ginput(1);
    [x0R,y0R,~]=ginput(1);
    AverageBaseVecY(:,1)=x0L;
    AverageBaseVecY(:,2)=y0L;
    AverageBaseVecY(:,3)=x0R;
    AverageBaseVecY(:,4)=y0R;
    plot(AverageBaseVecY(1,[1 3]),AverageBaseVecY(1,[2 4]),'b','LineWidth',1)
    w = waitforbuttonpress;
    fprintf('Baseline Estimation Method : Manual')
end

if Qreflection
    for i=1:camframetot
        [x01L,y01L,indexL]=findreflection([EdgeDropLR{i,1}.x,EdgeDropLR{i,1}.y],5,round(eLRDataD.ellipse{1,1}.Y0_in)+50);
        [x01R,y01R,indexR]=findreflection([EdgeDropLR{i,2}.x,EdgeDropLR{i,2}.y],5,round(eLRDataD.ellipse{1,1}.Y0_in)+50);
        AverageBaseVecY(i,1:4)=[x01L y01L x01R y01R];
    end
    fprintf('Baseline Estimation Method : Reflection')
end

if Qneedle
    fprintf('Needle Cutoff = %1.0i pixels\n',needlecutoff)
    AverageBaseVecY(:,1:2)=Lxy;
    AverageBaseVecY(:,3:4)=Rxy;
    plot(AverageBaseVecY(1,[1 3]),AverageBaseVecY(1,[2 4]),'b','LineWidth',1)
    fprintf('Baseline Estimation Method : Guiding Needle')
end
fprintf('\n')

%%
%---------------------------------------------------------------------------
%        4.       Fit contact angles using averaged baseline
%---------------------------------------------------------------------------
% predefine variables prior to loop
TLvec=zeros(camframetot,4);
PolyCAS=zeros(camframetot,2);
eLRCA=zeros(camframetot,2);
eCA=zeros(camframetot,2);
EllipseDIAM=zeros(camframetot,4);
edropvolume=zeros(camframetot,1);
edropvolumecap=zeros(camframetot,1);
edropvolumewholeellipse=zeros(camframetot,1);
edropvolumesphereref=zeros(camframetot,1);
edropvolumespherepoly=zeros(camframetot,1);
ellipsurf=zeros(camframetot,1);
Sa=zeros(camframetot,1);
ellipsurfsphere=zeros(camframetot,1);

% fit contact angles using the edges stored in EdgeCell as in example 1
disp('Contact angle calculations')
fprintf('Z Cutoff = %1.0i pixels\n',zpixelcutoff)
ii=1;
tababcdefg.a=[];tababcdefg.b=[];tababcdefg.c=[];tababcdefg.d=[];
tababcdefg.e=[];tababcdefg.f=[];tababcdefg.g=[];tababcdefg.dv=[];
tababcdefg.xL=[];tababcdefg.xR=[];tababcdefg.zL=[];tababcdefg.zR=[];

%% Ellipse
if ellipcal
    disp('Ellipse/Polynomial fitting')
    adsorptime=camduration;
    eTotalDataBaselineL=zeros(camframetot,2);
    eTotalDataBaselineR=zeros(camframetot,2);

    for n_frame=1:camframetot
        x0L=AverageBaseVecY(n_frame,1);y0L=AverageBaseVecY(n_frame,2);
        x0R=AverageBaseVecY(n_frame,3);y0R=AverageBaseVecY(n_frame,4);
        if min(EdgeCellLR{n_frame,1}.y)>min([y0L,y0R])-zpixelcutoff-caadsorpcutoff
            disp('')
            disp('adsorption detected.')
            fprintf('Adsorption frame = %1.0i (t=%1.3fs)\n', n_frame-1,(n_frame-1)/camfrate)
            adsorptime=(n_frame-1)/camfrate;
            break
        end
        %%
        if mod(n_frame,10^(max([1 floor(log10(camframetot)-2)])))==0
            fprintf('%1.0i ', n_frame)
            if mod(n_frame,10^(max([2 floor(log10(camframetot)-2)])))==0 || n_frame==camframetot
                fprintf('\n')
            end
        end
        %%
        % VOLUME + CA + Surface area + Surface contact area
        % fit edge of whole "L+R" to ellipsoid
        eTotalData=EllipticFitVolnotilt(EdgeCellLR{n_frame,1},EdgeCellLR{n_frame,1},[x0L,y0L],[x0R,y0R],zpixelcutoff);

        eTotalOutput.phi(n_frame,1)=eTotalData.ellipse{1,1}.phi;
        eTotalOutput.a(n_frame,1)=eTotalData.ellipse{1,1}.a; eTotalOutput.a(n_frame,2)=eTotalData.ellipse{1,2}.a;    %%%%aL %%%aR
        eTotalOutput.c(n_frame,1)=eTotalData.ellipse{1,1}.b; eTotalOutput.c(n_frame,2)=eTotalData.ellipse{1,2}.b;  %%%bL %%%bR
        eTotalOutput.X0_in(n_frame,1)=eTotalData.ellipse{1,1}.X0_in; eTotalOutput.Y0_in(n_frame,1)=eTotalData.ellipse{1,1}.Y0_in;

        eTotalDataBaselineL(n_frame,1:2)=[eTotalData.TLL ];
        eTotalDataBaselineR(n_frame,1:2)=[eTotalData.TLR ];

        %Recenter the line to x0,y0
        xy0LN1 = [x0L;y0L]-[eTotalOutput.X0_in(1,1);eTotalOutput.Y0_in(1,1)];
        xy0RN1 = [x0R;y0R]-[eTotalOutput.X0_in(1,1);eTotalOutput.Y0_in(1,1)];

        %Rotate
        xy0LN2 = [cos(eTotalOutput.phi(n_frame,1)) -sin(eTotalOutput.phi(n_frame,1)); sin(eTotalOutput.phi(n_frame,1)) cos(eTotalOutput.phi(n_frame,1))]...
            *[xy0LN1(1);xy0LN1(2)];
        xy0RN2 = [cos(eTotalOutput.phi(n_frame,1)) -sin(eTotalOutput.phi(n_frame,1)); sin(eTotalOutput.phi(n_frame,1)) cos(eTotalOutput.phi(n_frame,1))]...
            *[xy0RN1(1);xy0RN1(2)];

        %Recenter the line to xo(1)-x0(n_frame)
        xy0LN = xy0LN2-[eTotalOutput.X0_in(n_frame,1)-eTotalOutput.X0_in(1,1);eTotalOutput.Y0_in(n_frame,1)-eTotalOutput.Y0_in(1,1)];
        xy0RN = xy0RN2-[eTotalOutput.X0_in(n_frame,1)-eTotalOutput.X0_in(1,1);eTotalOutput.Y0_in(n_frame,1)-eTotalOutput.Y0_in(1,1)];

        eTotalOutput.lineX0LN(n_frame,1)=xy0LN(1); eTotalOutput.lineY0LN(n_frame,1)=xy0LN(2);
        eTotalOutput.lineX0RN(n_frame,1)=xy0RN(1); eTotalOutput.lineY0RN(n_frame,1)=xy0RN(2);

        %Parameters for calculating volume/surface area
        aV = eTotalOutput.a(n_frame,1); bV = aV; cV = eTotalOutput.c(n_frame,1); cVR = eTotalOutput.c(n_frame,2);
        x1=xy0LN(1);y1=0;z1=xy0LN(2); x2=xy0RN(1);y2=0;z2=xy0RN(2); x3=xy0RN(1);y3=1;z3=xy0RN(2);
        dV = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
        eV = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
        fV = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
        gV = dV*x1 + eV*y1 + fV*z1;

        %baseline
        m1=(z2-z1)/(x2-x1);
        b1=z1-m1*x1;
        hV=cV-(-b1);

        DV=  -abs(gV) / sqrt( (dV*aV)^2+(eV*bV)^2+(fV*cV)^2 );

        %volume, volumewhole
        edropvolume(n_frame,1) = aV*bV*cV*pi/3*( DV^3 - 3*DV +2 );
        edropvolumewholeellipse(n_frame,1) = 4/3*pi*aV*bV*cV;
        edropvolumecap(n_frame,1) = pi*aV*bV*hV^2/3/cV^2*(3*cV-hV);
        if b1<0
            edropvolume(n_frame,1)=edropvolumewholeellipse(n_frame,1)-edropvolume(n_frame,1);
        end

        %intersection of baseline and ellipse
        aITC=cV^2+aV^2*m1^2;
        bITC=2*aV^2*m1*b1;
        cITC=aV^2*(b1^2-cV^2);
        xL=(-bITC-sqrt(bITC^2-4*aITC*cITC))/(2*aITC);
        zL=xL*m1+b1;
        xR=(-bITC+sqrt(bITC^2-4*aITC*cITC))/(2*aITC);
        zR=xR*m1+b1;

        rSa=sqrt(((xR-xL)^2)+((zR-zL)^2))/2;

        %Surface area Sa
        Sa(n_frame,1) =pi*rSa^2;

        %CA
        m2L=-(xL*cV^2)/(aV^2*zL);
        m2R=-(xR*cV^2)/(aV^2*zR);
        if m2L>0
            eLRCA(n_frame,1)= (pi-abs(atan((m2L-m1)/(1+(m1*m2L)))))/pi*180;
        else
            eLRCA(n_frame,1)= (abs(atan((m2L-m1)/(1+(m1*m2L)))))/pi*180;
        end
        if m2R<0
            eLRCA(n_frame,2)= (pi-abs(atan((m2R-m1)/(1+(m1*m2R)))))/pi*180;
        else
            eLRCA(n_frame,2)= (abs(atan((m2R-m1)/(1+(m1*m2R)))))/pi*180;
        end
        CAx=mean(eLRCA(n_frame,1:2))/180*pi;

        %volumesphereref
        edropvolumesphereref(n_frame,1)=pi*rSa^3*(2-3*cos(CAx)+(cos(CAx))^3)/(3*(sin(CAx))^3); % edropvolumesphere(n_frame,1)=pi*rb^3*(2-3*cos(CAx)+(cos(CAx))^3)/(3*(sin(CAx))^3);

        %save these out
        tababcdefg.a(end+1,1)=aV; tababcdefg.b(end+1,1)=bV; tababcdefg.c(end+1 ,1)=cV; tababcdefg.d(end+1,1)=dV;
        tababcdefg.e(end+1,1)=eV; tababcdefg.f(end+1,1)=fV; tababcdefg.g(end+1,1)=gV; tababcdefg.dv(end+1,1)=DV;
        tababcdefg.xL(end+1,1)=xL; tababcdefg.xR(end+1,1)=xR; tababcdefg.zL(end+1,1)=zL; tababcdefg.zR(end+1,1)=zR;

        %surface area Sv
        zi = z1 - m1*x1;
        zi = -zi;
        ebsi = sqrt((cV^2-aV^2)/cV^2);
        ellipsurf(n_frame,1) = pi*aV*cV/ebsi*(  asin(ebsi) + ebsi*sqrt(1-ebsi^2) - ( asin(ebsi*zi/aV)+ebsi*zi/aV*sqrt(1-(ebsi*zi/aV)^2)   ));
        ellipsurfsphere(n_frame,1) = 4*pi*((aV^3.2+2*(aV*cV)^1.6)/3)^(1/1.6);

        %%
        % VOLUME + CA + Surface area + Surface contact area
        % fit edge of "separate L / R" to ellipsoid
        warning('off')
        eLRData=EllipticFitVol(EdgeCellLR{n_frame,1},EdgeCellLR{n_frame,2},[x0L,y0L],[x0R,y0R],zpixelcutoff);
        warning('on')

        eLROutput.phi(n_frame,1)=eLRData.ellipse{1,1}.phi;
        eLROutput.a(n_frame,1)=eLRData.ellipse{1,1}.a; eLROutput.a(n_frame,2)=eLRData.ellipse{1,2}.a;   %%%%aL %%%aR
        eLROutput.c(n_frame,1)=eLRData.ellipse{1,1}.b; eLROutput.c(n_frame,2)=eLRData.ellipse{1,2}.b;  %%%bL %%%bR
        eLROutput.X0_in(n_frame,1)=eLRData.ellipse{1,1}.X0_in; eLROutput.Y0_in(n_frame,1)=eLRData.ellipse{1,1}.Y0_in;

        aVL = eLROutput.a(n_frame,1);
        aVR = eLROutput.a(n_frame,2);
        bVL = aVL;
        bVR = aVR;
        cVL = eLROutput.c(n_frame,1);
        cVR = eLROutput.c(n_frame,2);

        % CA
        m11=((xy0RN(2,:)-xy0LN(2,:))/(xy0RN(1,:)-xy0LN(1,:)));

        b11=xy0RN(2,:)-m11*xy0LN(2,:);

        aITCL=cVL^2+aVL^2*m11^2;
        bITCL=2*aVL^2*m11*b11;
        cITCL=aVL^2*(b11^2-cVL^2);
        xL=(-bITCL-sqrt(bITCL^2-4*aITCL*cITCL))/(2*aITCL);
        zL=xL*m11+b11;

        aITCR=cVR^2+aVR^2*m11^2;
        bITCR=2*aVR^2*m11*b11;
        cITCR=aVR^2*(b11^2-cVR^2);
        xR=(-bITCR+sqrt(bITCR^2-4*aITCR*cITCR))/(2*aITCR);
        zR=xR*m11+b11;

        m22L=-(xL*cVL^2)/(aVL^2*zL);
        m22R=-(xR*cVR^2)/(aVR^2*zR);

        baseL=[xL,zL];
        baseR=[xR,zR];

        if m22L>0
            eCA(n_frame,1)= 180-(abs(atan((m11-m22L)/(1+(m11*m22L)))) *180 /pi);
        else
            eCA(n_frame,1)= (abs(atan((m11-m22L)/(1+(m11*m22L)))) *180 /pi);
        end
        if m22R<0
            eCA(n_frame,2)= 180-(abs(atan((m11-m22R)/(1+(m11*m22R)))) *180 /pi);
        else
            eCA(n_frame,2)= (abs(atan((m11-m22R)/(1+(m11*m22R)))) *180 /pi);
        end
        %%
        if Qrotstage
            figure(3);
            hold on
            t=linspace(-3,3);
            plot((xy0RN(1)-xy0LN(1))*t+xy0LN(1),(xy0RN(2)-xy0LN(2))*t+xy0LN(2),'color',cc(ii,:),'LineWidth',1); hold on
            radius=60;
            imrotation=atand((y0R-y0L)/(x0R-x0L));
            %             ellipserimlr=@(e,t) cos(t)*e.a+sin(t)*e.b;
            ellipserimlr=@(e,t) ones(size(t))*[0,0]+cos(t)*[1,0]*e.a(n_frame)+sin(t)*[0,1]*e.c(n_frame);
            PlotEllipseDataLR=ellipserimlr(eTotalOutput,linspace(0,2*pi,1000)');
            plot(PlotEllipseDataLR(:,1),PlotEllipseDataLR(:,2),'color',cc(ii,:),'LineWidth',1)
            scatter(tababcdefg.xL(n_frame,1),tababcdefg.zL(n_frame,1))
            scatter(tababcdefg.xR(n_frame,1),tababcdefg.zR(n_frame,1))
            axis ij
            axis equal
            xlim([min(PlotEllipseDataLR(:,1))-50 max(PlotEllipseDataLR(:,1))+50])
            ylim([min(PlotEllipseDataLR(:,2))-50 max(PlotEllipseDataLR(:,2))+50])
            ylabel('z (pixel)');
            xlabel('x (pixel)');
            axis square
            box on
        else
            if  n_frame==1 % plot ellipse data
                figure(3); subplot(2,2,1)
                imagesc(mov{n_frame})
                colormap gray
                axis ij
                axis equal
                xlim([1 camwidth]);
                ylim([1 camheight]);
                hold on
                title('elliptical L/R')
                t=linspace(-3,3);
                plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r','LineWidth',1)
                radius=60;
                imrotation=atand((y0R-y0L)/(x0R-x0L));
                ellipserim=@(e,t) ones(size(t))*[e.X0_in,e.Y0_in]+cos(t)*[cos(-e.phi),sin(-e.phi)]*e.a+sin(t)*[-sin(-e.phi),cos(-e.phi)]*e.b;
                lw=2;
                plot([eLRData.TLR(1),eLRData.TLR(1)-radius*cosd(eLRData.CAR+imrotation)],[eLRData.TLR(2),eLRData.TLR(2)-radius*sind(eLRData.CAR+imrotation)],'r','LineWidth', 1)
                plot([eLRData.TLL(1),eLRData.TLL(1)+radius*cosd(eLRData.CAL-imrotation)],[eLRData.TLL(2),eLRData.TLL(2)-radius*sind(eLRData.CAL-imrotation)],'b','LineWidth', 1)
                PlotEllipseDataL=ellipserim(eLRData.ellipse{1},linspace(0,2*pi,1000)');
                PlotEllipseDataR=ellipserim(eLRData.ellipse{2},linspace(0,2*pi,1000)');
                plot(PlotEllipseDataL(:,1),PlotEllipseDataL(:,2),'c','LineWidth',1)
                plot(PlotEllipseDataR(:,1),PlotEllipseDataR(:,2),'m','LineWidth',1)
            end

            if n_frame == 1
                figure(3); subplot(2,2,2)
                imagesc(mov{n_frame})
                axis ij
                axis equal
                hold on
                colormap gray
                title('elliptical total')
                xlim([1 camwidth]);
                ylim([1 camheight]);
            end
            if n_frame ==1 || (mod(n_frame,plotint)==0 && isfield(eLRData,'TLR')==1)
                figure(3); subplot(2,2,1)
                t=linspace(-3,3);
                plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r','LineWidth',1)
                radius=60;
                imrotation=atand((y0R-y0L)/(x0R-x0L));
                ellipserim=@(e,t) ones(size(t))*[e.X0_in,e.Y0_in]+cos(t)*[cos(-e.phi),sin(-e.phi)]*e.a+sin(t)*[-sin(-e.phi),cos(-e.phi)]*e.b;
                lw=2;
                plot([eLRData.TLR(1),eLRData.TLR(1)-radius*cosd(eLRData.CAR+imrotation)],[eLRData.TLR(2),eLRData.TLR(2)-radius*sind(eLRData.CAR+imrotation)],'r','LineWidth', 1)
                plot([eLRData.TLL(1),eLRData.TLL(1)+radius*cosd(eLRData.CAL-imrotation)],[eLRData.TLL(2),eLRData.TLL(2)-radius*sind(eLRData.CAL-imrotation)],'b','LineWidth', 1)
                PlotEllipseDataL=ellipserim(eLRData.ellipse{1},linspace(0,2*pi,1000)');
                PlotEllipseDataR=ellipserim(eLRData.ellipse{2},linspace(0,2*pi,1000)');
                plot(PlotEllipseDataL(:,1),PlotEllipseDataL(:,2),'c','LineWidth',1)
                plot(PlotEllipseDataR(:,1),PlotEllipseDataR(:,2),'m','LineWidth',1)
            end
            if n_frame == 1
                figure(3); subplot(2,2,2)
                imagesc(mov{end})
                hold on
                colormap gray
                %             title('elliptical all frames')
                axis ij
                axis equal
                xlim([1 camwidth]);
                ylim([1 camheight]);
                ylabel('z (pixel)');
                xlabel('x (pixel)');
            end
            if n_frame ==1 || mod(n_frame,plotint)==0
                figure(3); subplot(2,2,2)
                t=linspace(-3,3);
                plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r','LineWidth',1)
                radius=60;
                imrotation=atand((y0R-y0L)/(x0R-x0L));
                ellipserimlr=@(e,t) ones(size(t))*[e.X0_in,e.Y0_in]+cos(t)*[cos(-e.phi),sin(-e.phi)]*e.a+sin(t)*[-sin(-e.phi),cos(-e.phi)]*e.b;
                lw=2;
                PlotEllipseDataLR=ellipserimlr(eTotalData.ellipse{1},linspace(0,2*pi,1000)');
                plot(PlotEllipseDataLR(:,1),PlotEllipseDataLR(:,2),'color',cc(ii,:),'LineWidth',1.5)
            end
            if n_frame ==1 || mod(n_frame,plotint)==0
                figure(3); subplot(2,2,3)
                hold on
                t=linspace(-3,3);
                plot((xy0RN(1)-xy0LN(1))*t+xy0LN(1),(xy0RN(2)-xy0LN(2))*t+xy0LN(2),'color',cc(ii,:),'LineWidth',1); hold on
                radius=60;
                imrotation=atand((y0R-y0L)/(x0R-x0L));
                ellipserimlr=@(e,t) ones(size(t))*[0,0]+cos(t)*[1,0]*e.a(n_frame)+sin(t)*[0,1]*e.c(n_frame);
                PlotEllipseDataLR=ellipserimlr(eTotalOutput,linspace(0,2*pi,1000)');
                plot(PlotEllipseDataLR(:,1),PlotEllipseDataLR(:,2),'color',cc(ii,:),'LineWidth',1)
                scatter(tababcdefg.xL(n_frame,1),tababcdefg.zL(n_frame,1))
                scatter(tababcdefg.xR(n_frame,1),tababcdefg.zR(n_frame,1))
                axis ij
                axis equal
                xlim([min(PlotEllipseDataLR(:,1))-50 max(PlotEllipseDataLR(:,1))+50])
                ylim([min(PlotEllipseDataLR(:,2))-50 max(PlotEllipseDataLR(:,2))+50])
                ylabel('z (pixel)');
                xlabel('x (pixel)');
                axis square
                box on
            end

            if Qcrop
                if n_frame ==1 || mod(n_frame,plotint)==0
                    figure(3); subplot(2,2,4)
                    if n_frame ==1
                        imagesc(mov{1}); hold all
                        t=linspace(-3,3);
                        plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'b','LineWidth',1)
                        plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L-zpixelcutoff,'m','LineWidth',1)
                        plot([xi xi xf xf xi], [yf yi yi yf yf],'r','LineWidth',1) % plot crop box of the drop
                        axis ij
                        axis equal
                        hold on
                        colormap gray
                        title('edgecutoff')
                        xlim([1 camwidth]);
                        ylim([1 camheight]);
                    end
                    edgecutoff=max([y0L y0R])-zpixelcutoff;
                    scatter(EdgeCellLR{n_frame,1}.x(EdgeCellLR{n_frame,1}.y<edgecutoff),EdgeCellLR{n_frame,1}.y(EdgeCellLR{n_frame,1}.y<edgecutoff),2)
                end
            end
        end

        %% Poly
        if polycal
            if n_frame ==1
                warning('off')
                [PolyData,PolyDataCoef]=fit(edgeLR.x(edgeLR.y<min([y0L y0R])-10),edgeLR.y(edgeLR.y<min([y0L y0R])-10),'a*x^5+b*x^4+c*x^3+d*x^2+e*x+f','tolx',1e-9);
                warning('on')
            end
            [PolyData,PolyDataCoef]=fit(edgeLR.x(edgeLR.y<min([y0L y0R])-10),edgeLR.y(edgeLR.y<min([y0L y0R])-10),'a*x^5+b*x^4+c*x^3+d*x^2+e*x+f',...
                'startpoint',[PolyData.a PolyData.b PolyData.c PolyData.d PolyData.e PolyData.f]);
            PolyData=polyfit(edgeLR.x(edgeLR.y<min([y0L y0R])-10),edgeLR.y(edgeLR.y<min([y0L y0R])-10),5)

            if n_frame ==1
                figure(100)
                imagesc(mov{1}); hold all
                axis ij
                axis equal
                colormap gray
                xlim([1 camwidth]);
                ylim([1 camheight]);
            end
            if n_frame ==1 || mod(n_frame,plotint)==0
                figure(100)
                h=plot(PolyData,edgeLR.x(edgeLR.y<min([y0L y0R])-10),edgeLR.y(edgeLR.y<min([y0L y0R])-10)); hold all
                set(h(1),'color',cc(ii,:),'LineWidth',0.1)
                set(h(2),'color',cc(ii,:),'LineWidth',1.2)
                legend off
            end
        end
        %%
        if mod(n_frame,plotint)==0 || n_frame==1
            if ii==25
                ii=0;
            else
                ii=ii+1;
            end
        end
    end
end

avcabegin=[];
for n=1:camframetot
    if timest(1,n)<=contactcalsec
        avcabegin(1,end+1)=timest(1,n);
    end
end

[eLRCAD,Dposit,~]=deleteoutliers(mean(eLRCA,2));

if Qrotstage
    fprintf('Contact angle (1s sep. L/R) = %1.2f\n',mean(mean(eCA(1:size(avcabegin',1),1:2))))
    fprintf('Contact angle for the first frame (whole L+R) = %1.2f\n',mean(mean(eLRCA(1:size(avcabegin',1),1:2))))
else
    fprintf('Contact angle (1s sep. L/R) = %1.2f, (L-R=%1.2f)\n',mean(mean(eCA(1:size(avcabegin',1),1:2))),-diff(mean(eCA(1:size(avcabegin',1),1:2))))
    fprintf('Contact angle (all frame, whole L+R, delete outliers) = %1.2f\n',(mean(eLRCAD(1:size(avcabegin',1),1:1))))
end

%%
%---------------------------------------------------------------------------
%        5.       Plot result for video analysis
%---------------------------------------------------------------------------

%%
if polycal
    figure
    p1=plot(timest,PolyCAS(:,1),'r--');hold on
    p2=plot(timest,PolyCAS(:,2),'r-');
end
if  polycal
    figure
    p11=plot(timest,edropvolumespherepoly,'b-');
    hYLabel1=ylabel('Volume (pixel^3)');
    hXLabel=xlabel('Time (s)');
    legend on
    legend('poly')
end

%% PLOTS and Illus.

timesec=timest';
if ellipcal
    % contact new
    eLavg= mean(eLRCAD,2);
    figure;
    timeplot=timesec;
    timeplot(Dposit,:)=[];
    p4=plot(timeplot,eLavg,'r-','linewidth',1.2);
    ylabel('Contact angle (degree)');
    xlabel('Time (s)');
    xlim([0 camduration])
    ylim([0 180])
    set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.02,0.02])
    set(gcf,'windowstyle','normal')
    axis square
    box on

    %Pic gray scale
    figure;
    imagesc(mov{1})
    axis ij
    axis equal
    hold on
    colormap gray
    xlim([1 camwidth]);
    ylim([1 camheight]);
    xlim([0 camwidth]);
    ylim([0 camheight]);
    set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.02,0.02])
    set(gcf,'windowstyle','normal')
    xlabel('x (pixel)')
    ylabel('y (pixel)')

    %Rescale baseline
    figure;
    imagesc(mov{n_frame})
    colormap gray
    axis ij
    axis equal
    xlim([1 camwidth]);
    ylim([1 camheight]);
    hold on
    title('elliptical')
    t=linspace(-3,3);
    plot((x0R-x0L)*t+x0L,(y0R-y0L)*t+y0L,'r--','LineWidth',2)
    radius=60;
    imrotation=atand((y0R-y0L)/(x0R-x0L));
    ellipserim=@(e,t) ones(size(t))*[e.X0_in,e.Y0_in]+cos(t)*[cos(-e.phi),sin(-e.phi)]*e.a+sin(t)*[-sin(-e.phi),cos(-e.phi)]*e.b;
    lw=2;
    plot([eTotalData.TLR(1),eTotalData.TLR(1)-radius*cosd(eTotalData.CAR+imrotation)],[eTotalData.TLR(2),eTotalData.TLR(2)-radius*sind(eTotalData.CAR+imrotation)],'r','LineWidth', 1)
    plot([eTotalData.TLL(1),eTotalData.TLL(1)+radius*cosd(eTotalData.CAL-imrotation)],[eTotalData.TLL(2),eTotalData.TLL(2)-radius*sind(eTotalData.CAL-imrotation)],'b','LineWidth', 1)
    PlotEllipseDataL=ellipserim(eTotalData.ellipse{1},linspace(0,2*pi,1000)');
    PlotEllipseDataR=ellipserim(eTotalData.ellipse{2},linspace(0,2*pi,1000)');
    plot(PlotEllipseDataL(:,1),PlotEllipseDataL(:,2),'c','LineWidth',1)
    plot(PlotEllipseDataR(:,1),PlotEllipseDataR(:,2),'m','LineWidth',1)

    %Scatter base of every frame
    movSum=zeros(camheight,camwidth);
    for i=1:camframetot
        movSum=movSum+double(mov{i});
    end
    movSum=movSum/camframetot;
    figure(7);
    imagesc(movSum)
    axis ij
    axis equal
    hold on
    colormap gray
    xlim([1 camwidth]);
    ylim([1 camheight]);
    xlim([0 camwidth]);
    ylim([0 camheight]);
    set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.02,0.02])
    set(gcf,'windowstyle','normal')
    xlabel('x (pixel)')
    ylabel('z (pixel)')
    scatter(eTotalDataBaselineL(:,1),eTotalDataBaselineL(:,2),20,'r','filled')
    scatter(eTotalDataBaselineR(:,1),eTotalDataBaselineR(:,2),20,'r','filled')
end

save([pathnamm  fnam '-CASV.mat'])