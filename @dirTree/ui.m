function h = ui(obj)
cpy = copy(obj);

h = figure;
leftbox = uicontrol('parent', h, 'style', 'listbox', 'string', obj.lsnames, 'max',2, 'min',0);
rightbox = uicontrol('parent', h, 'style', 'listbox', 'string', obj.lsnames, 'max',2, 'min',0);
commitbtn = uicontrol('parent', h, 'style', 'pushbutton', 'string', '==>');
finishbtn = uicontrol('parent', h, 'style', 'pushbutton', 'string', 'Done!');

set(h, 'resizefcn', {@resizecallback, leftbox, rightbox, commitbtn, finishbtn});
h.UserData.lvl = 0; % directory level saved in UserData
h.UserData.obj = cpy;
set(commitbtn, 'callback', {@branch, leftbox});
set(leftbox, 'callback', {@update, leftbox, rightbox});
resizecallback(h, [], leftbox, rightbox, commitbtn, finishbtn);
branch(leftbox, [], leftbox);
update(leftbox, [], leftbox, rightbox);

function plotls(obj, box, lvl)
box.Value = 1;
ls = obj.lsnames(lvl);
set(box, 'string', unique(ls));

function branch(h, ~, box)
h.Parent.UserData.lvl = h.Parent.UserData.lvl + 1;
h.Parent.UserData.ls = h.Parent.UserData.obj.lsnames(h.Parent.UserData.lvl);
plotls(h.Parent.UserData.obj, box, h.Parent.UserData.lvl);

function update(h, ~, source, target)
temp = copy(h.Parent.UserData.obj);
idx = setxor(1:length(source.String), source.Value);
[~,idx] = intersect(h.Parent.UserData.ls, source.String(idx));
temp.children(idx) = [];
plotls(temp, target, h.Parent.UserData.lvl + 1);

% function finish(h, obj)
% obj = 

function resizecallback(h, ~, leftbox, rightbox, commitbtn, finishbtn)
figsize = h.Position;
lcoor = [.05 .1 .4 .8] .* repmat(figsize(3:4), [1 2]);
rcoor = [.55 .1 .4 .8] .* repmat(figsize(3:4), [1 2]);
bcoor = [.475 .475 .05 .05] .* repmat(figsize(3:4), [1 2]);
fcoor = [.475 .1 .05 .05] .* repmat(figsize(3:4), [1 2]);
leftbox.Position = round(lcoor);
rightbox.Position = round(rcoor);
commitbtn.Position = round(bcoor);
finishbtn.Position = round(fcoor);