function patchfacefcn(obj,hit,fun,selection,dLim)
% callback function for patch face selection
%
% PATCHFACEFCN(obj,hit,callbackfun,selection)
%
% The function can find the index of the face in a patch object which was
% clicked on by the mouse. The function should be used as a callback
% function for the ButtonDownFcn event of patch object and it will call a
% user defined function with the same arguments as the ButtonDownFcn call,
% but adding an extra argument, the face index. Thus the user defined
% callback function will have the following header:
%
%           callbackfun(obj,hit,faceIndex)
%
% The function can detect if the mouse click was on a face or on an edge of
% the patch object.
%
% Input:
%
% obj           The patch object that calls the function.
% hit           Hit object that contains the point where the object was hit.
% callbackfun   User function that will be called in case of a click event
%               on obj. It should have the following header:
%                   callbackfun(obj,hit,faceIndex)
%               where face Index contains the index of the face that was
%               clicked on, it can contain a single index or more depending
%               on the selection type.
% selection     String defines three diferent selection criteria when the
%               callbackfun() function will be called:
%                   'face'  The callbackfun() will be triggered if a face
%                           was clicked (numel(faceIndex)==1).
%                   'edge'  The callbackfun() will be triggered if an edge
%                           is clicked (numel(faceIndex)==2).
%                   'all'   The callbackfun() will be triggered for both
%                           faces and edges.
%
%
% Example: 
%   The color of any face of the red triangulated icosahedron will be
%   changed to green if clicked on.
%
% mesh = icomesh(1);
% V = mesh.Points;
% F = mesh.ConnectivityList;
% hPatch = patch('Faces',F,'Vertices',V,'FaceColor','r','EdgeColor','none');
% axis equal
% box on
% view(3)
% camlight('right')
% hPatch.FaceColor       = 'flat';
% hPatch.FaceVertexCData = repmat([1 0 0],[size(F,1) 1]);
% fun = @(obj,hif,face)set(obj,'FaceVertexCData',[obj.FaceVertexCData(1:(face-1),:); [0 1 0]; obj.FaceVertexCData((face+1):end,:)]);
% hPatch.ButtonDownFcn = @(obj,hit)swplot.patchfacefcn(obj,hit,fun,'face');
%
% See also PATCH.
%

% point where we hit the surface
P = hit.IntersectionPoint;
%disp('In function')
nF = size(obj.Faces,1);

% vertices of patch
V = obj.Vertices;

if size(V,2) ==2
    flat = true;
    % for 2D patch
    P = P(1:2);
else
    flat = false;
end

if size(obj.Faces,2) ~=3
    error('patchfacefcn:WrongObject','The pathcfacefcn() callback function works only for triangulated surfaces!');
end

% precision for finding planes of faces
if nargin < 5
    dLim = 1e-7;
end

% shift the origin to the first vertex of every triangle
E1 = V(obj.Faces(:,2),:)-V(obj.Faces(:,1),:);
E2 = V(obj.Faces(:,3),:)-V(obj.Faces(:,1),:);
D  = bsxfun(@minus,P,V(obj.Faces(:,1),:));

if ~flat
    % check whether the plane of any face contains the point
    det = sum(cross(D,E1,2).*E2,2);
    
    % find points within the plane
    pIdx = find(abs(det)<dLim);
    
    % determine barycentric coordinates
    bCoord = zeros(numel(pIdx),2);
    for ii = 1:numel(pIdx)
        bCoord(ii,:) = ([E1(pIdx(ii),:)' E2(pIdx(ii),:)']\D(pIdx(ii),:)')';
    end
    
    % find the right face(s)
    fIdx = pIdx(all(bCoord>=0 & bCoord<=1 & sum(bCoord,2)<=1,2));
else
    % flat patch, all triangles are in the plane
    % determine barycentric coordinates
    bCoord = zeros(nF,2);
    for ii = 1:nF
        bCoord(ii,:) = ([E1(ii,:)' E2(ii,:)']\D(ii,:)')';
    end
    
    % find the right face(s)
    fIdx = find(all(bCoord>=0 & bCoord<=1 & sum(bCoord,2)<=1,2));
end
%sprintf('fIdx is %d\n',fIdx)
switch selection
    case 'all'
        % include clicks on edges (numel(fIdx)>1)
        if numel(fIdx) > 0
            fun(obj,hit,fIdx);
        end
    case 'face'
        % only trigger for faces
        if numel(fIdx) == 1
            fun(obj,hit,fIdx);
         %   disp('face')
            obj.UserData(fIdx)=not(obj.UserData(fIdx));
        end
    case 'edge'
        % only trigger for faces
        if numel(fIdx) == 2
            fun(obj,hit,fIdx);
            disp('edge')
        end
    otherwise
        error('pathcfacefcn:WrongInput','The pathcfacefcn() callback has only two modes: ''body'' and ''face''!');
end

end