function f=segmentImage(imageFileName)
        
    realImage = imread(imageFileName);    
    %figure;imshow(realImage);    
    %mainImage = imresize(mainImage, 0.5);
    
    
    binaryImage = im2bw(realImage,graythresh(realImage));    
    %figure;imshow(binaryImage);
    
    mainImage = imcomplement(binaryImage);
    %figure;imshow(mainImage);
       
    [row, column]=size(mainImage);
    
    setM(11);
    setH(13);
    zoneWidthPercent=0.05;    
    
    verticalZoneWidth=ceil(column*zoneWidthPercent);            
    verticalZone = cell(1/zoneWidthPercent,1) ;
    
    %diving into zone
    zoneNo=1;
    for i = 1:verticalZoneWidth:column        
        verticalZone{zoneNo}=mainImage(1:row,i:min(i+verticalZoneWidth,column));           
        zoneNo=zoneNo+1;
    end
    zoneNo=zoneNo-1;
    
    %figure;imshow(verticalZone{1});        
    validZone = cell(zoneNo,1);        
    validThreshold=50;    
    
    for i=1:1:zoneNo
        if(getWhitePixelDensityInTheZone(verticalZone{i})<validThreshold)
            validZone{i}=0;
        else 
            validZone{i}=1;
        end        
    end
    
    PR = cell(zoneNo,1);        
    for i=1:1:zoneNo
        PR{i}=projectionProfile(verticalZone{i});        
    end
        
    
    SPR=projectionProfileAfterSmoothing(PR,validZone);
    
    %{
    figure;hold on;plot(PR{2},'b');    
    plot(SPR{2},'r');
    hold off;
    %}
  
    delSPR=derivative(SPR);  
    %delSPR=derivative(PR);
        
    %figure;plot(delSPR{2},'b');    
    
    crudeSegmentation=findCrudeVerticalSegmentationline(delSPR,validZone);       
    grayImage=totalCrudeSegementedImage(crudeSegmentation,verticalZone,validZone);
    %figure;imshow(grayImage);
    
    
    refinedSegementation=findRefinedSegmentation(crudeSegmentation, verticalZone,validZone);
    grayImage=totalCrudeSegementedImage(refinedSegementation,verticalZone,validZone);
    %figure;imshow(grayImage);
    
    trickSegmentation=trickRefinedSegmentation(crudeSegmentation,verticalZone,validZone);
    grayImage=totalCrudeSegementedImage(trickSegmentation,verticalZone,validZone);
    %figure;imshow(grayImage);
    
    %Text line seperator drawing algorithm
    
    %Association phase       
    
    finalSeperator=finalSeperation(trickSegmentation,validZone);    
    grayImage=drawSeperator(finalSeperator,verticalZone,validZone);
    %figure;imshow(grayImage);
    
    %trickSegmentation{3}
    %finalSeperator{3}
    
    refinedSeperator=refineSeperator(finalSeperator,verticalZone,validZone);    
    grayImage=drawSeperator(refinedSeperator,verticalZone,validZone);
    %figure;imshow(grayImage);
    
    [validSeperator, validcolumn]=getValidIndex(refinedSeperator,validZone);
    grayImage=drawSeperator(refinedSeperator,verticalZone,validZone);
    %figure;imshow(grayImage);
    
    
    %validSeperator
    [pixeltoline,lineNo]=assignPixelToLines(validSeperator, validcolumn,verticalZone);    
    
    lines=linebylineAnalysis(pixeltoline,lineNo,validSeperator, validcolumn,verticalZone);
    
    
    for i=1:1:lineNo
        imshow(lines{i});
        pause(1);
    end
    
    
    linewords=segmentWord(lines,lineNo);
    
    %linewords{1}
    
    
    for i=1:1:1
       lineImage=linewords{i};
       
       for j=1:1:size(lineImage,2);
           image=lineImage{j}; 
           imshow(image);
           pause(1);
       end
    end
    
    %}
    
    %pixeltoline=assignPixeltotextline
        
    %lines=lineDrawing(trickSegmentation,verticalZone,validZone);
    
   % close all;    
    %crudeSegmentation{3}
    %refinedSegementation{3}
  
    %grayImage=drawCrudeSegmentation(crudeSegmentation{3},verticalZone{3});
    %figure;imshow(grayImage);    
    %grayImage=drawCrudeSegmentation(refinedSegementation{3},verticalZone{3});
    %figure;imshow(grayImage);    
   % close all;
   
   f=0;
end

function [lineWords]=segmentWord(lines,lineNo)

    lineWords=cell(lineNo,1);
    
    for i=1:1:lineNo
        line=lines{i};
        words=extractWord(line);
        lineWords{i}=words;
    end

    lineWords;
end

function [words]=extractWord(line)
        
    %imshow(line);
   %{ 
    line=[
        0 1 1 0 0 0 1;
        0 0 0 1 0 0 0;
        1 1 0 1 0 0 0;
        0 1 0 0 0 0 1;];
    %}
    cc=findAllConnectedComponents(line);                
    word=groupCC(cc);
    
    [r,c]=size(line);
    
    words=wordtoimage(r,c,word);
       
    
    
    
end

function [image]=wordtoimage(r,c,word)
    [rw,cw]=size(word);    
    image=cell(rw,cw);
    
    for i=1:1:cw        
       image{i}=makeImage(r,c,word{i}) ;
    end
    
    
end

function [image]=makeImage(r,c,word)
    word=sortrows(word,[-2]);
    maxcol=word(1,2);
    
    word=sortrows(word,[2]);
    mincol=word(1,2);
    
    
    
    
    image=zeros(r,maxcol-mincol);
    
    for i=1:1:size(word,1)
        x=word(i,1);
        y=word(i,2);
        image(x,y-mincol+1)=255;
    end
    
    
   
    
end


function word=groupCC(cc)
    
    word=[];
    
    [r,c]=size(cc);
        
    index=1;
    
    word{index}=cc{1};
    
    for i=2:1:c;
        %dist=checkConntectivity(word{index},cc{i});        
        
        dist=checkConntectivity(cc{i-1},cc{i});        
                          
        if dist<15
            word{index}=[word{index};cc{i}];
        else
            index=index+1;
            word{index}=cc{i};
        end             
        
    end
        
end

function mindist=checkConntectivity(ci,cj)
 
    x=[1 -1];
    y=[2 -2];
    
    mindist=10000;
    
    for i=x
        for j=y
           ci=sortrows(ci,[i j]);
           corneri=ci(1,:);
           for k=x
               for l=y
                   cj=sortrows(cj,[k l]);
                   cornerj=cj(1,:);
                   
                   dx=corneri(1)-cornerj(1);
                   dy=corneri(2)-cornerj(2);
                   
                   dist=sqrt(dx*dx+dy*dy);
                   
                   if(dist<mindist)
                       mindist=dist;
                   end
               end
           end
                       
        end
    end
    
    
    
   

end


function [surrogate,cc]=DFS(surrogate,i,j,cc)
    x=[-1,-1,-1, 0, 0, 1, 1, 1];
    y=[-1, 0, 1,-1, 1,-1, 0, 1];
    
    [r,c]=size(surrogate);
    
    if i<1 || i>r || j<1 || j>c
       return 
    end
    
    if surrogate(i,j)==0
        return;            
    end
    
    if surrogate(i,j)==1
        cc=[cc;i j];
        surrogate(i,j)=-1;
        for k=x
            for l=y
               [surrogate,cc]=DFS(surrogate,i+k,j+l,cc); 
            end
        end
    end
    return;
end

function cc=findAllConnectedComponents(line)
    
    [r,c]=size(line);
    
    surrogate=line;
    cc=[];
    ccno=1;
    
    
    for j=1:1:c
        for i=1:1:r
            if surrogate(i,j)==1
                ccf=[];
                [surrogate,ccf]=DFS(surrogate,i,j,ccf);
                cc{ccno}=ccf;
                ccno=ccno+1;
            end
        end
    end    
    
    %surrogate
    %{
    for i=1:1:size(cc,2)        
        cc{i}
    end
    %}

end

function lines=linebylineAnalysis(pixeltoline,lineNo,validSeperator, validcolumn,verticalZone)
    
    [rs,cs]=size(validSeperator);
    
    lines=cell(lineNo,1);
    for i=1:1:lineNo
       min=100000000;
        max=-1;       
       for j=1:1:rs
            seperator=validSeperator{j};
                                   
            lmin=seperator{i};
            lmax=seperator{i+1};         
          
            if min>lmin
                min=lmin;
            end
            if max<lmax
                max=lmax;
            end
       end     
       
       
       line=extractline(min,max,pixeltoline,i,validcolumn);
       lines{i}=line;
    end    
end

function line=extractline(lmin,lmax,pixeltoline,lineNo,validcolumn)
    %pixeltoline
     
    
    [rp,rc]=size(pixeltoline);    
    [rd,cd]=size(pixeltoline{1});
    
    rowsize=lmax-lmin;
    columnsize=rp*cd;
    
    line=zeros(rowsize,columnsize);
    
    r=1;
    c=1;
    
    for i=lmin:1:lmax
        c=1;
        for j=1:1:rp
            currentseg=pixeltoline{j};
            for k=1:1:cd
               if currentseg(i,k)==lineNo
                   line(r,c)=1;
               end
               c=c+1; 
            end
        end
        
        r=r+1;
    end
    %}
    %line=[];
end

function [assignment,lineNo]=assignPixelToLines(validSeperator, validcolumn,verticalZone)
    [r,c]=size(validcolumn);               
    towhichline=cell(c,1);
    lineNo=size(validSeperator{validcolumn{1}},2);
    
    for i=1:1:c        
        zone=verticalZone{validcolumn{i}};
        seperator=validSeperator{i};
       
        % error('stop');
        %zone=[0 0 1 0 0 0;0 1 1 1 0 0;0 0 0 0 1 0;]
       
        lp=size(seperator,2);
        lineNo=lp;
        [rz,cz]=size(zone);
        linezone=zeros(rz,cz);
        
        
        for j=1:1:lp-1
           line1=[seperator{j} seperator{j+1}];           
           for k=line1(1):1:line1(2)
                for p=1:1:cz
                    if zone(k,p)==1
                        linezone(k,p)=j;
                    end
                end
           end
        end
        
        towhichline{i}=linezone;
        
    end
    
    %towhichline{validcolumn{1}};
    lineNo=lineNo-1;
    assignment=towhichline;
end

function [seperator,validIndex]=getValidIndex(refinedSeperator,validZone)
    validIndex=[];
    indexNo=1;
    for i=1:1:size(validZone,1)       
        if validZone{i}==1                        
            validIndex{indexNo}=i;
            indexNo=indexNo+1;
        end
    end           
    indexNo=indexNo-1;
    seperator=cell(indexNo,1);
    %{
        for i=1:1:indexNo
            seperator{i}=refinedSeperator{validIndex{i}};
        end
    %}
    
    lineNo=size(refinedSeperator{validIndex{1}},2);
    
    for i=1:1:indexNo
    %for i=4:1:4
        seperator{i}=filterSeperator(refinedSeperator{validIndex{i}},lineNo);
    end
    
end

function seperator=filterSeperator(refinedSeperator,lineNo)
    seperator=refinedSeperator;
        
    
    
    if size(refinedSeperator,2)==lineNo
        return
    elseif size(refinedSeperator,2)<lineNo
        %error('change image');        
        return
    else 
        while size(refinedSeperator,2)~=lineNo
            
            value=abs(refinedSeperator{1}-refinedSeperator{size(refinedSeperator,2)});
            index=-1;

            for i=1:1:size(refinedSeperator,2)-1
                if abs(refinedSeperator{i}-refinedSeperator{i+1})<value
                    value=abs(refinedSeperator{i}-refinedSeperator{i+1});
                    index=i;
                end
            end
            candidate1=index
            candidate2=index+1
            finalCandidate=index;
            if candidate1==1 
                finalCandidate=2;
            elseif candidate2==size(refinedSeperator,2)
                finalCandidate=candidate1;
            else
                d1=abs(refinedSeperator{candidate1}-refinedSeperator{candidate1-1});
                d2=abs(refinedSeperator{candidate1}-refinedSeperator{candidate2+1});
                fd1=max(d1,d2);
                d1=abs(refinedSeperator{candidate2}-refinedSeperator{candidate1-1});
                d2=abs(refinedSeperator{candidate2}-refinedSeperator{candidate2+1});
                fd2=max(d1,d2);

                if fd1>fd2
                    finalCandidate=candidate1;
                else 
                    finalCandidate=candidate2;
                end            
            end

            refinedSeperator(finalCandidate)=[];
        end
    end
    
    seperator=refinedSeperator;
    
end

function refinedSeperator=refineSeperator(finalSeperator,verticalZone,validZone)
    textLines=cell(size(validZone));
    validIndex=[];
    indexNo=1;
    for i=1:1:size(validZone,1)       
        if validZone{i}==1                        
            validIndex{indexNo}=i;
            indexNo=indexNo+1;
        end
    end           
    indexNo=indexNo-1;
    
    p1RefinedSeperator=finalSeperator;    
    p1RefinedSeperator{validIndex{1}}=finalSeperator{validIndex{1}};
    
    
    
    for in=2:1:indexNo                
   % for in=15:1:15                
        i1=validIndex{in};
        i0=validIndex{in-1};        
        p1seperator=update(in,i1,i0,finalSeperator,verticalZone);
        %p1seperator                       
        p1RefinedSeperator{in}=num2cell(p1seperator);        
    end
    
    %p1RefinedSeperator{15}
    
            
    p2RefinedSeperator=p1RefinedSeperator;
    p2RefinedSeperator{validIndex{indexNo}}=p1RefinedSeperator{validIndex{indexNo}};
    
    
    
    for in=indexNo-1:-1:1                
        i0=validIndex{in+1};
        i1=validIndex{in};         
       
        p2seperator=update(in,i1,i0,p1RefinedSeperator,verticalZone);
        %p1seperator                        
        p2RefinedSeperator{in}=num2cell(p2seperator);        
    end
    
    refinedSeperator=p2RefinedSeperator;
end

function p1seperator=update(in,i1,i0,finalSeperator,verticalZone)
        
        %i1=validIndex{in};
        %i0=validIndex{in-1};        
        sepi1=num2cell(sort(cell2mat(finalSeperator{i1})));
        sepi0=num2cell(sort(cell2mat(finalSeperator{i0})));      
        
        %sepi1={27, 100, 170}        
        %sepi0={26, 81, 125, 169}
        
        closestSepi1=[];
        
        for i=1:1:size(sepi1,2)
            closestSepi1{i}=[sepi1{i} 0];
        end
                
        
        for i=1:1:size(sepi0,2)
            i0val=sepi0{i};
            [index,value]=findClosest(i0val,sepi1);
            scell=closestSepi1{index};
            scell(2)=scell(2)+1;
            closestSepi1{index}=[scell,i0val];
        end
        %refine phase
        closestSepi1{4};
        
        
        p1seperator=[];
        
        for i=1:1:size(sepi1,2)
            scell=closestSepi1{i};
            
            if scell(2)==1 || i==1 || i==size(sepi1,2)
                p1seperator=[p1seperator,scell(1)];
            elseif scell(2)>1
                for l=1:1:scell(2)
                    ith=2+l;
                    
                    if scell(ith)<scell(1)
                        partition=[sepi1{i} sepi1{i-1}];
                    else
                        partition=[sepi1{i+1} sepi1{i}];
                    end
                    npart=getnewSeperator(partition,scell(ith),verticalZone{in});
                    p1seperator=[p1seperator,npart];
                end
            end            
        end
end

function npart=getnewSeperator(partition,ith,zone)
      %zone=vzone(partition(1):partition(2),:);
     % partition
      [r,c]=size(zone);
      min=10000;
      k=-1;
      for i=partition(2):1:partition(1)
         row=zone(i,:);
         pm=(sum(row(:))/c);
         dm=abs(ith-i)/(partition(1)-partition(2));
         qm=(pm+1)*(dm+1);
         if qm<min            
             min=qm;
             k=i;
         end
      end
      
      npart=k;
end

function f=evalate()
    f=0;
end
    
function [index,mainvalue]=findClosest(i0val,sepi1);
    value=abs(i0val-sepi1{1});
    mainvalue=sepi1{1};
    index=1;
    for i=1:1:size(sepi1,2)
        if abs(i0val-sepi1{i})<value
            value=abs(i0val-sepi1{i});
            mainvalue=sepi1{i};
            index=i;
        end
    end
end



function grayImage=drawSeperator(finalSeperator,verticalZone,validZone)
    grayImage=[];
    for i=1:1:size(validZone,1)
        if validZone{i}==1
            grayImage=horzcat(grayImage,zeros(size(grayImage,1),1),drawCrudeSeperator(	finalSeperator{i},verticalZone{i}));        
        end
    end
    
end

function grayImage=drawCrudeSeperator(segmentation,zone)
    grayImage = 255 * uint8(zone);    
    for i=1:1:size(segmentation,2)       
       for j=1:1:size(zone,2)
         grayImage(segmentation{i},j)=255;
       end
    end        
end

function finalSeperator=finalSeperation(trickSegmentation,validZone)
    
    finalSeperator=cell(size(trickSegmentation));

    for i=1:1:size(validZone,1)       
        if validZone{i}==1                        
            segment=trickSegmentation{i};
            %seperator=cell(size(segment),1);
            index=1;
            for j=1:1:size(segment,1)
                if segment(j,3)==-1
                    seperator{index}=round((segment(j,2)+segment(j,1))/2);
                    index=index+1;
                end
            end
            finalSeperator{i}=seperator;
        else
            finalSeperator{i}=[];
        end
    end    
end


function textLines=lineDrawing(trickSegmentation,verticalZone,validZone)
    textLines=cell(size(validZone));
    validIndex=[];
    indexNo=1;
    for i=1:1:size(validZone,1)       
        if validZone{i}==1                        
            validIndex{indexNo}=i;
            indexNo=indexNo+1;
        end
    end           
    indexNo=indexNo-1;
    %{
    mainSeperators=cell(indexNo,1);    
    for k=1:1:size(validIndex,2)
        i=validIndex{k};
        tricksegement=trickSegmentation{i}
        
        [rs,cs]=size(tricksegement);
        
        seperator=cell(rs+1,1);        
        
        seperator{1}=tricksegement(1,1);
        for j=1:1:rs
            seperator{j+1}=tricksegement(j,2);
        end
        
        mainSeperators{i}=seperator;
    end        
       %}
    
    
end


function refinedSegmentation=trickRefinedSegmentation(crudeSegmentation,verticalZone,validZone)
    
    refinedSegmentation=cell(size(crudeSegmentation));
    
   for i=1:1:size(validZone,1)    
   %for i=3:1:3
        if validZone{i}==1                        
            refinedSegmentation{i}=trick(crudeSegmentation{i},verticalZone{i});
        else
            refinedSegmentation{i}=[];
        end
    end           
end

function Z=trick(Y,yZone)
    X=Y;     
    for i=2:size(Y,1)-2
        zone=yZone(Y(i,1):Y(i,2),:);      
        [r,c]=size(zone);
        value=sum(zone(:));        
        
        if Y(i-1,3)==-1 && value<10            
            X(i,3)=-1;
        elseif Y(i,3)==1 && abs(Y(i,2)-Y(i,1))<7 && (abs(Y(i+1,2)-Y(i+1,1))<7) && (abs(Y(i+2,2)-Y(i+2,1))>7)
            Y(i+1,3)=1;
        elseif Y(i,3)==-1 && Y(i-1,3)==1 && abs(Y(i,2)-Y(i,1))<7
            X(i,3)=1;
        elseif Y(i,3)==1 && (Y(i-1,3)==1 || Y(i+1,3)==1) && abs(Y(i,2)-Y(i,1))<7
            X(i,3)=1;
        end
        
    end
    
    Z(1,1:3)=X(1,1:3);
    index=1;
    for i=2:1:size(Y,1)  
        
        if X(i,3)== Z(index,3)
            Z(index,2)=X(i,2);            
        else            
            index=index+1;
            Z(index,1:3)=X(i,1:3);
        end
    end

end

function refinedSegmentation=findRefinedSegmentation(crudeSegmentation, verticalZone,validZone)
    %states
    C0=1; %text
    C1=-1; %gapRegion
    
    %start probabilities/initial probabilites
    phi0=0.5;
    phi1=0.5;
    
    [m0,m1]=getMeanHeight(crudeSegmentation,validZone); %m0=mean height of text region m1=mean height of gap region
     
    %m1=m0;
    %disp(transition(1,-1,[15 15],15))
    %transition can be calculated don't worry
    
    %emission probability
    [mu0,sigma0, mu1, sigma1, num0,num1]=getemissionProbablityParameter(crudeSegmentation,validZone,verticalZone,m0,m1);
    
    %bluh=crudeSegmentation{3};
    %emission([mu0,sigma0, mu1, sigma1],-1,bluh(3,1:3),verticalZone{3});
    
    refinedSegmentation=cell(size(crudeSegmentation));
    
   for i=1:1:size(validZone,1)
   % for i=3:1:3
        if validZone{i}==1            
            S=[C0,C1];
            Y=crudeSegmentation{i};
            Yzone=verticalZone{i};
            A=[m0,m1];            
            B=[mu0,sigma0,mu1,sigma1,num0,num1];
            phi=[0.5,0.5];
            refinedSegmentation{i}=VITERBI(S,phi,Y,Yzone,A,B);
        else
            refinedSegmentation{i}=[];
        end
    end        
    %}
end

function X=VITERBI(S,phi,Y,Yzone,A,B)
    K=size(S,2);
    T=size(Y,1);
    
    T1=cell(K,T);
    T2=cell(K,T);
    
    stateIterator=[1,2];    
    
    for i=stateIterator
        T1{i,1}=phi(i)*emission(B,S(i),Y(1,1:3),Yzone);
        T2{i,1}=2;
    end
    
    for i=2:1:T
        for j=stateIterator
            max=-1;
            best=-1;
            for k=stateIterator
                value=T1{k,i-1}*transition(S(k),S(j),A,abs(Y(i,2)-Y(i,1)))*emission(B,S(j),Y(i,1:3),Yzone);
                if(max<value)
                    max=value;
                    best=k;
                end
            end
            if best==-1
                disp('something went wrong at viterbi');
               % exit(0);
            end
            T1{j,i}=max;
            T2{j,i}=best;
        end
    end
    
    Zt=-1;
    max=-1;
    
    for i=stateIterator
        if T1{i,T}>max
            Zt=i;
            max=T1{i,T};
        end
    end
    
    X=Y;    
    %X(T,3)=S(Zt);
    for i=T:-1:1
        Zt=T2{Zt,i};
        X(i,3)=S(Zt);
    end
    
end

function p=transition(Ci,Cj,A,H)
    if Ci==-1 && Cj==-1
        p=exp(-H/A(2));
    elseif Ci==-1 && Cj==1
        p=1-exp(-H/A(2));
    elseif Ci==1 && Cj==-1
        p=1-exp(-H/A(1));
    else 
        p=exp(-H/A(1));
    end    
   % disp('Transition');
   % Ci
   % Cj
   % H
   % p
end

function p=emission(B,S,Y,Yzone)
   % Y
    value=Yzone(Y(1,1):Y(1,2),:);
    [r,c]=size(value);
    
   % disp('emission');
 %   S
 %   Y
    
    if S==-1
        X=(sum(value(:)))/(r*c);                
        if X==0 
            p=1;
        else            
            Xi=log(X);
           %[Xl,Xu]=closest(Xi,B(6));
           %p=normcdf((Xi+Xu)/2,B(3),B(4))-normcdf((Xi+Xl)/2,B(3),B(4));
           p=normpdf(Xi,B(3),B(4))/normpdf(B(3),B(3),B(4));
        end
    else
        X=(sum(value(:)))/(r*c);
        if X==0
            p=0;
        else                       
           Xi=log(X);
           %[Xl,Xu]=closest(Xi,B(5));
           %p=normcdf(Xi+(Xu-Xi)/2,B(1),B(2))-normcdf(Xi-(Xi-Xl)/2,B(1),B(2));           
           p=normpdf(Xi,B(1),B(2))/normpdf(B(1),B(1),B(2));
        end
    end    
  %  p
end

function [Xl,Xu]=closest(X,num)
    num=sort(num);
        Xl=X;
        Xu=X;               
    for i=1:1:size(num,2)
        if X>num(i)
            Xl=num(i);
        else
            Xu=num(i);
            break;
        end
    end
end

function p=getNormValue(x,mu,sigma)
    p1=1/(sigma*sqrt(2*pi));
    p2=exp(-((x-mu)*(x-mu))/(2*sigma*sigma));
    
    p=p1*p2;

end


function [mu0,sigma0, mu1, sigma1,num0,num1]=getemissionProbablityParameter(crudeSegmentation,validZone,verticalZone,m0,m1)
    
    num0=[];
    num1=[];
    
    for i=1:1:size(validZone,1)
        if validZone{i}==1
            segment= crudeSegmentation{i};
            zone=verticalZone{i};
            
            for j=1:1:size(segment,1)
                value=zone(segment(j,1):segment(j,2),:);
                [r,c]=size(value);
                if segment(j,3)==-1    
                    density=(sum(value(:)))/(r*c);
                    if r>=m1*.2 && density>0
                        num1=[num1,log(density)];                                                
                    end
                else
                    density=(sum(value(:)))/(r*c);
                    if r>=m0*.2 && density>0
                        num0=[num0,log(density)];
                    end
                end                
            end
        end            
    end            
    mu0=mean(num0);
    mu1=mean(num1);
    sigma0=std(num0);
    sigma1=std(num1);    
    %figure;hist([num0,num1],40);
    
end

function [m0,m1]=getMeanHeight(crudeSegementation,validZone)
    num0=[];
    num1=[];
    
    for i=1:1:size(validZone,1)
        if validZone{i}==1
            segment= crudeSegementation{i};
            
            for j=1:1:size(segment,1)
                if abs(segment(j,2)-segment(j,1))>5
                if segment(j,3)==-1                    
                    num1=[num1,abs(segment(j,2)-segment(j,1))];                                        
                else
                    num0=[num0,abs(segment(j,2)-segment(j,1))];
                end                
                end
            end
        end            
    end
    
    m0=mean(num0);
    m1=mean(num1);
    
end

function grayImage=totalCrudeSegementedImage(crudeSegmentation,verticalZone,validZone)
    grayImage=[];
    for i=1:1:size(validZone,1)
        if validZone{i}==1
            grayImage=horzcat(grayImage,zeros(size(grayImage,1),1),drawCrudeSegmentation(crudeSegmentation{i},verticalZone{i}));        
        end
    end
end
function grayImage=drawCrudeSegmentation(segmentation,zone)
    grayImage = 255 * uint8(zone);    
    for i=1:1:size(segmentation,1)       
       for j=segmentation(i,1):1:segmentation(i,2)
          if segmentation(i,3)==-1
              color=64;
          else
              color=256-64;
          end
          
          for k=1:1:size(grayImage,2)
              if ~(segmentation(i,3)==1 && grayImage(j,k)==255)                               
                grayImage(j,k)=color;
              end
              
          end
       end
       
        
    end
    
    
end

function crude=findCrudeVerticalSegmentationline(delSPR,validZone)
    crude=cell(size(delSPR));
    for i=1:1:size(delSPR,1)
        if validZone{i}==1
            crude{i}=findMinMax(delSPR{i});
        else
            crude{i}=[];
        end 
    end
end


%-1 for blank 1 for text region
function crudei=findMinMax(delSPRi)
    crudei=zeros(1,3);
    %crudei(1:1,1:3)=[0,0,-1];    
    maxmin=-1;
    start=1;
    flag=-1;
    index=1;
    
    for i=1:1:size(delSPRi,1)
        if flag==-1
            if delSPRi(i)>=maxmin
                maxmin=delSPRi(i);
            else
                crudei(index,1:3)=[start, i, flag];
                flag=1;
                start=i;
                index=index+1;
            end
            
        else
            if delSPRi(i)<=maxmin
                maxmin=delSPRi(i);
            else
                crudei(index,1:3)=[start,i,flag];
                flag=-1;
                start=i;
                index=index+1;
            end   
        end
    end
end


function setM(m)
  global M;
  M = m;
end

function m = getGlobalM
    global M
    m = M;
end

function setH(h)
  global H;
  H = h;
end

function h = getGlobalH
    global H;
    h = H;
end
    
function d=getWhitePixelDensityInTheZone(image)    
    d = sum(image(:));
end

function d=projectionProfile(image)
    d = sum(image,2);    
end

function SPR=projectionProfileAfterSmoothing(PR,validZone)
   [r,c]=size(PR);
   SPR=cell(r,c);
   %M=4;
   M=getGlobalM();
   for i=1:1:r
       sum=0;
       for j=-M:1:M
            if (i+j)>0 && (i+j)<=r    
                sum=sum+validZone{i+j}*getWeight(j,M)*PR{i+j};
            end
       end 
       SPR{i}=sum;
   end
end

function w=getWeight(i,M)

p1=exp((-3*abs(i))/(M+1));
p2=0;

    for k=-M:1:M
       p2=p2+exp((-3*abs(k))/(M+1)); 
    end
    w=p1/p2;
end


function delSPR=derivative(SPR)
    
    delSPR=cell(size(SPR));
    for i=1:1:size(SPR,1)       
       delSPR{i}=delSPRi(SPR{i}); 
    end        
end

function SPRj=delSPRi(SPR)    
    [r,c]=size(SPR);
    SPRj=zeros(r,c);
    %h=13;
    h=getGlobalH();
    p1=1/(h*(h+1));
    for j=1:1:r
       p2=0;
       for k=1:1:h
           if (j+k)<=r && (j-k)>0                             
             p2=p2+k*(SPR(j+k)-SPR(j-k));
           end
       end
       SPRj(j)=p1*p2;
    end
end