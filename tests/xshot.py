"""Screenshot helper for tests/test_gui.m (Linux/X11 only): xshot.py DISPLAY out.png "X=.. Y=.. WIDTH=.. HEIGHT=.."."""
import sys, re
from PIL import ImageGrab
disp, out, geom = sys.argv[1], sys.argv[2], sys.argv[3]
g = dict((k, int(v)) for k, v in re.findall(r'(X|Y|WIDTH|HEIGHT)=(-?\d+)', geom))
im = ImageGrab.grab(xdisplay=disp)
if 'WIDTH' in g:
    im = im.crop((max(g['X'], 0), max(g['Y'], 0), g['X'] + g['WIDTH'], g['Y'] + g['HEIGHT']))
im.save(out)
