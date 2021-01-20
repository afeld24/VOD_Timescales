%Input center coordinate (either lat or lon) and this linearly interpolates
%the 4 side and 4 corner coordinates
function [RightBound,LeftBound,LowerBound,UpperBound,...
   TopLeftCorner,TopRightCorner,BottomLeftCorner,BottomRightCorner]...
                = InterpEdgesCorners(Xmat)
Mid = Xmat(:,1:end-1)+diff(Xmat,1,2)/2;
RightBound = cat(2,Mid,NaN(size(Xmat,1),1));
LeftBound = cat(2,NaN(size(Xmat,1),1),Mid);

Mid = Xmat(1:end-1,:)+diff(Xmat,1,1)/2;
LowerBound = cat(1,Mid,NaN(1,size(Xmat,2)));
UpperBound = cat(1,NaN(1,size(Xmat,2)),Mid);

Up = UpperBound(:,1:end-1)+diff(UpperBound,1,2)/2;
TopLeftCorner = cat(2,NaN(size(Xmat,1),1),Up);
TopRightCorner = cat(2,Up,NaN(size(Xmat,1),1));

Down = LowerBound(:,1:end-1)+diff(LowerBound,1,2)/2;
BottomLeftCorner = cat(2,NaN(size(Xmat,1),1),Down);
BottomRightCorner = cat(2,Down,NaN(size(Xmat,1),1));
end