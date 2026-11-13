function added = gif_ensure_loop(file)
% GIF_ENSURE_LOOP  make sure a GIF file loops forever: it must contain the NETSCAPE 2.0 application extension
%   with loop count 0. If it is missing it is inserted right after the global colour table (the place every
%   viewer expects it) and the file is rewritten. Returns true if the block had to be added.
%   (Some viewers - e.g. the Windows photo viewer - ignore it and play a GIF once anyway; browsers, GitHub,
%   PowerPoint, VLC and most image viewers loop.)
fid = fopen(file,'rb'); b = fread(fid, Inf, 'uint8=>uint8').'; fclose(fid);
added = false;
if numel(b) < 14 || ~strcmp(char(b(1:3)),'GIF'), error('gif_ensure_loop: %s is not a GIF file', file); end
sig = uint8('NETSCAPE2.0');
hit = strfind(char(b), char(sig));
if ~isempty(hit)
    % the sub-block is 03 01 lo hi (loop count, 0 = forever)
    p = hit(1) + numel(sig);
    if p+3 <= numel(b) && b(p) == 3 && b(p+1) == 1 && b(p+2) == 0 && b(p+3) == 0, return, end
    if p+3 <= numel(b), b(p+2:p+3) = 0; end                       % force loop count 0 (= forever)
    fid = fopen(file,'wb'); fwrite(fid,b,'uint8'); fclose(fid); added = true; return
end
flags = double(b(11));
pos = 13;                                                       % 6 header + 7 logical screen descriptor bytes
if bitand(flags,128), pos = pos + 3*2^(bitand(flags,7)+1); end    % global colour table
blk = [uint8([33 255 11]) sig uint8([3 1 0 0 0])];
b = [b(1:pos) blk b(pos+1:end)];
fid = fopen(file,'wb'); fwrite(fid,b,'uint8'); fclose(fid);
added = true;
end
