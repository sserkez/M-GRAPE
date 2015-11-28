function dist=distcut(dist,xlim,pxlim,ylim,pylim,tlim,glim,mode)

%dist=beam;
sizedist=numel(dist.T);

index.o=false(sizedist,1);
index.T1=index.o;
index.T2=index.o;
index.G1=index.o;
index.G2=index.o;
index.X1=index.o;
index.X2=index.o;
index.PX1=index.o;
index.PX2=index.o;
index.Y1=index.o;
index.Y2=index.o;
index.PY1=index.o;
index.PY2=index.o;
% 
%  index.X1=dist.X<=-1e-4;
%  index.X2=dist.X>=1e-4;
%  index.PX1=dist.PX<=-1.5e-5;
%  index.PX2=dist.PX>=1.5e-5;
%  
%  index.Y1=dist.Y<=-1e-4;
%  index.Y2=dist.Y>=1e-4;
%  index.PY1=dist.PY<=-1.5e-5;
%  index.PY2=dist.PY>=1.5e-5;
if numel(xlim)==2
 index.X1=dist.X<=xlim(1);
 index.X2=dist.X>=xlim(2);
end
if numel(pxlim)==2
 index.PX1=dist.PX<=pxlim(1);
 index.PX2=dist.PX>=pxlim(2);
end
if numel(ylim)==2
 index.Y1=dist.Y<=ylim(1);
 index.Y2=dist.Y>=ylim(2);
end
if numel(pylim)==2
 index.PY1=dist.PY<=pylim(1);
 index.PY2=dist.PY>=pylim(2);
end
if numel(tlim)==2
 index.T1=dist.T<=tlim(1);
 index.T2=dist.T>=tlim(2);
end
if numel(glim)==2
 index.G1=dist.G<=glim(1);
 index.G2=dist.G>=glim(2);
end


index.f=index.T1|index.T2|index.X1|index.X2|index.PX1|index.PX2|index.Y1|index.Y2|index.PY1|index.PY2|index.G1|index.G2;
if strcmp(mode,'leave')
elseif strcmp(mode,'delete')
    index.f=~index.f;
else
    error('cut mode not defined: "leave" or "delete"')
end
    %disp(' ');
    disp([num2str(sum(index.f)./sizedist.*100),'% particles is cut off']);
    %disp('% is cut off');
    %%
    dist.T(index.f)=[];
    dist.PX(index.f)=[];
    dist.X(index.f)=[];
    dist.Y(index.f)=[];
    dist.PY(index.f)=[];
    dist.G(index.f)=[];

    dist.charge;
    dist.charge=dist.charge./sizedist.*sum(~index.f);

end



%clear index sizedist

