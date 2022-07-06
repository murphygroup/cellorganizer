'''
Utility for debugging from Stack Overflow.

Taraz Buck
2021-10-10
'''


# https://stackoverflow.com/questions/13174412/python-start-interactive-debugger-when-exception-would-be-otherwise-thrown

## {{{ http://code.activestate.com/recipes/65287/ (r5)
# code snippet, to be included in 'sitecustomize.py'
import sys

def info(type, value, tb):
   if hasattr(sys, 'ps1') or not sys.stderr.isatty():
      # we are in interactive mode or we don't have a tty-like
      # device, so we call the default hook
      sys.__excepthook__(type, value, tb)
   else:
      import traceback, pdb
      # we are NOT in interactive mode, print the exception...
      traceback.print_exception(type, value, tb)
      print
      # ...then start the debugger in post-mortem mode.
      pdb.pm()

sys.excepthook = info
## end of http://code.activestate.com/recipes/65287/ }}}
