"""Timestamped progress log shared by the engine, the command line and the web pages.

The engine reports work by calling `progress(message)`.  A message that starts with `[k/T]` opens step k of
T (it drives the progress bar); other messages are detail lines belonging to the step above them.
"""
import re
import time

_STEP = re.compile(r"^\[(\d+)/(\d+)\]")


class StepLog:
    """Callable progress sink that keeps every message with the seconds elapsed since it was created.

    `sinks` are called as sink(line, fraction) for every message, e.g. to redraw a progress bar.
    """

    def __init__(self, *sinks):
        self.t0 = time.time()
        self.lines = []
        self.fraction = 0.0
        self.steps = 0
        self.sinks = list(sinks)

    def elapsed(self):
        return time.time() - self.t0

    def __call__(self, msg):
        m = _STEP.match(msg)
        if m:
            self.steps += 1
            self.fraction = min(1.0, (int(m.group(1)) - 1) / int(m.group(2)))   # earlier steps are done
        line = f"[{self.elapsed():6.1f}s] {msg}"
        self.lines.append(line)
        for sink in self.sinks:
            sink(line, self.fraction)

    def finish(self, msg):
        """Record the closing message and set the bar to full."""
        self.fraction = 1.0
        self(msg)

    def text(self, last=None):
        return "\n".join(self.lines[-last:] if last else self.lines)
