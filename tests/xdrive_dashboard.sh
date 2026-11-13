#!/bin/bash
# Helper of tests/test_gui.m (Linux, X11, xdotool): operates the controls of a dashboard window with real mouse clicks
# while the animation loop is running:  xdrive_dashboard.sh DISPLAY shotfile.png flagfile
# Sequence: faster (+) twice, Pause, Play, switch to the third tab, screenshot, Stop; then waits for the flag file and
# presses Play once more (the loop is restarted from the button).
# Click positions are fractions of the Octave Qt window (menu + toolbar above the figure area).
export DISPLAY=$1
for i in $(seq 1 120); do WID=$(xdotool search --name "results" 2>/dev/null | head -1); [ -n "$WID" ] && break; sleep 0.5; done
[ -z "$WID" ] && exit 1
sleep 3
eval $(xdotool getwindowgeometry --shell $WID)
clk() { xdotool mousemove $(python3 -c "print(int($X+$1*$WIDTH), int($Y+$2*$HEIGHT))") click 1; sleep 0.8; }
BAR=0.939; TABS=0.0955
clk 0.589 $BAR     # +
clk 0.589 $BAR     # +
clk 0.062 $BAR     # Pause
clk 0.062 $BAR     # Play
clk 0.565 $TABS    # third tab
sleep 1.5
/usr/bin/python3 "$(dirname "$0")/xshot.py" $DISPLAY "$2" "X=$X Y=$Y WIDTH=$WIDTH HEIGHT=$HEIGHT"
clk 0.928 $BAR     # Stop
for i in $(seq 1 120); do [ -f "$3" ] && break; sleep 0.5; done
sleep 1
clk 0.062 $BAR     # Play: restarts the loop from the callback
