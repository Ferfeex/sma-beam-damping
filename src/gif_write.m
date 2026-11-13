function gif_write(file, frames, map, fps)
% GIF_WRITE  write an animated GIF that LOOPS FOREVER (pure Octave/MATLAB, no external tool).
%   gif_write(file, frames, map, fps)   frames: cell array of indexed images (uint8, 0-based indices into map),
%   map: n x 3 colour map (0..1), fps: frames per second (GIF delays are multiples of 1/100 s).
%   The first imwrite carries 'LoopCount', Inf; the following frames are appended (Octave's 4-D multi-frame
%   imwrite repeats the first frame - checked with Octave 10.3 - while 'append' works in Octave and MATLAB).
%   Afterwards gif_ensure_loop checks that the file really contains the NETSCAPE 2.0 loop block and inserts
%   it if the imwrite back-end (GraphicsMagick build, MATLAB version) left it out.
d = fileparts(file);
if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end
if exist(file,'file'), delete(file); end
delay = round(100/fps)/100;
uselzw = true;      % Octave (GraphicsMagick) writes GIFs UNCOMPRESSED unless asked: 'Compression','lzw' cuts the size ~10x
for j = 1:numel(frames)
    if j == 1, args = {'LoopCount', Inf}; else, args = {'WriteMode','append'}; end
    args = [{'DelayTime', delay} args];
    if uselzw
        try
            imwrite(frames{j}, map, file, 'gif', args{:}, 'Compression', 'lzw');
        catch
            uselzw = false;
            if j > 1 && exist(file,'file') == 0, error('gif_write: could not write %s', file); end
        end
    end
    if ~uselzw, imwrite(frames{j}, map, file, 'gif', args{:}); end
end
gif_ensure_loop(file);
end
