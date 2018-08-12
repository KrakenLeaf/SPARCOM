function alignEquals
% Select text from the editor and run this function to align equals sign
% You may add a shortcut to this function in the Matlab toolbar
%
% Alessandro Masullo - 2014

% Get the selected text in the active editor
aed = matlab.desktop.editor.getActive;
sel = aed.Selection;
startSel = matlab.desktop.editor.positionInLineToIndex(aed, sel(1), sel(2));
endSel = matlab.desktop.editor.positionInLineToIndex(aed, sel(3), sel(4));

% Whole text
mainText = aed.Text;
PreText = mainText(1:startSel-1);
PosText = mainText(endSel+1:end); 
% Selected text
SelText = mainText(startSel:min(endSel,numel(mainText)));
% Split selected text in lines
Lines = strsplit(SelText,char(10),'CollapseDelimiters',false);

% Text parsing
for i = 1:numel(Lines)
    if ~isempty(strfind(Lines{i},'=')) && ...                   % Check if contains equals
            isempty(regexp(Lines{i},'^\s*\%.*','match','once')) % Check for comments
        % Remove the spaces before and after the first equal sign
        Lines{i} = regexprep(Lines{i},'\s*=\s*','=','once');
        % Splits the line in pre-equal and post-equal
        bf = strsplit(Lines{i},'=');
        pre{i} = bf{1};
        % Select everything after the equal sign
        pos{i} = regexp(Lines{i},'=.*','match','once');
        pos{i} = pos{i}(2:end); % Without the equal sign
        % Strtrim to avoid indentation spaces
        len(i) = length(strtrim(regexprep(pre{i},'\s*=\s*','=','once')));
    else
        % Skip this line
        pre{i} = [];
    end
end
max_len = max(len);

% Align the text
AlignedText = '';
for i = 1:numel(Lines)
    if ~isempty(pre{i})
        % Repeat space to fill the gap
        pre{i} = [pre{i},repmat(' ',1,max_len-len(i))];
        AlignedText = [AlignedText, pre{i},' = ',pos{i},char(10)];
    else
        % Skip this line
        AlignedText = [AlignedText, Lines{i},char(10)];
    end
end

% Replace the text. end-1 is to delete the last \n in AlignedText
aed.Text = [PreText AlignedText(1:end-1) PosText];
% Re-select the text
[sel(3), sel(4)] = matlab.desktop.editor.indexToPositionInLine(aed, numel([PreText AlignedText(1:end-1)]));
aed.Selection = sel;
