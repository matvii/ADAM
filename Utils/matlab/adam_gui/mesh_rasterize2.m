function  I= mesh_rasterize2(tlist,vlist,visible,d,I_size )
%Based on patch_rasterize by Andrea.Tagliasacchi@epfl.ch
%p(1),p(2) image size in pixels p(1) is vertical, p(2) horizontal
%d pixel size, d(1) is vertical, d(2) horizontal
Vs = vlist(:,1:2);
Fs = tlist;
%keyboard
d=fliplr(d);
dist=vlist(:,3);
% Build triangulation data structure (for barycentric)
Vss=[Vs(:,1)/d(1),Vs(:,2)/d(2)];
TR = triangulation(Fs,Vss);

% Top-right vertex gives image size

%I_size = fliplr( ceil( max(Vs)./d) );
%zero location
z=[I_size(1)/2+1,I_size(2)/2+1];

% outside of triangulation you have NAN
I = nan( I_size );
Idist=inf(I_size);

% for every face
for ti=1:size(Fs,1)
    if visible(ti)==0
        continue;
    end
  
    F = Fs(ti,:);
    
    % get this triangle's bounding box   
    vs = Vs(F,:);
    pmin = floor( min(vs)./d );
    pmax = ceil( max(vs)./d );
 
    % implicit coordinates for triangle
    [XX,YY] = meshgrid( pmin(1):pmax(1), pmin(2):pmax(2) );
    XXYY = [XX(:),YY(:)];
    TI = ti*ones(size(XXYY,1),1);
    %XXYY=[XXYY(:,1)*d(1) XXYY(:,2)*d(2)]; %Convert to vlist coordinate frame
    % get barycentric & check interior
    bary = TR.cartesianToBarycentric(TI,XXYY); %Convert to vlist coordinate frame
    inside = ( all( bary>=0, 2 ) & all( bary<=1, 2 ) );
    
    % set face index inside the image
   % keyboard
    idxs = sub2ind(I_size',XXYY(inside,2)+z(1),XXYY(inside,1)+z(2));
    Vd=dist(F);
    Pd=bary*Vd;
    Pdi=Pd(inside);
    idc=Idist(idxs)>Pdi;
    Idist(idxs(idc))=Pdi(idc);
   % if ~isempty(idxs)
    %    keyboard
   % end
    I(idxs(idc)) = ti;
end
I(isnan(I))=0;



end

