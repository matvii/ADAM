%facets - logical list of fixed facets
nfac=size(tlist,1);
nvert=size(vlist,1);
invis=false(1,nfac);
invis(facets)=true;
%keyboard
tol=min(max(vlist(:,1))-min(vlist(:,1)))*1e-4;
figure;
P = trisurf(tlist,vlist(:,1),vlist(:,2),vlist(:,3)); axis equal
%P.FaceColor       = 'flat';
fcolor=repmat([1,0,0],nfac,1);
fcolor(invis,1)=0;
fcolor(invis,2)=1;
fcolor(invis,3)=0;
P.FaceVertexCData=fcolor;
P.UserData=invis;
fun = @(obj,hif,face)set(obj,'FaceVertexCData',[obj.FaceVertexCData(1:(face-1),:); double([not(obj.FaceVertexCData(face,1)) not(obj.FaceVertexCData(face,2)) 0]); obj.FaceVertexCData((face+1):end,:)]);
P.ButtonDownFcn = @(obj,hit) patchfacefcn(obj,hit,fun,'face',tol);

fixfacets=P.UserData;
invis=unique(reshape(tlist(fixfacets,:),[],1)); %fixed vertices