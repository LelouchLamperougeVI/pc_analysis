function s = parseXML(filename, targets)
% Extract attributes and data from target fields

try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end
try
    s = crawl(tree, targets);
catch
    error('Unable to parse XML file %s.',filename);
end


function [s, name] = crawl(node, targets)
name = char(node.getNodeName);
if any( cellfun(@(x) strcmpi(name, x), targets) )
    s = struct(struct('Attributes', struct(), 'Data', ''));
    if any(strcmp(methods(node), 'getData'))
        s.Data = char(node.getData);
    end
    if node.hasAttributes
        attributes = node.getAttributes;
        for count = 1:attributes.getLength
            s.Attributes.(char(attributes.item(count-1).getName)) = char(attributes.item(count-1).getValue);
        end
    end
else
    s = [];
end

if node.hasChildNodes
    children = node.getChildNodes;
    for ii = 1:children.getLength
        [rets, retn] = crawl(children.item(ii-1), targets);
        if ~isempty(rets)
            s.(retn) = rets;
        end
    end
end