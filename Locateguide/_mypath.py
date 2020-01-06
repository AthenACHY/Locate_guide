import os, sys
thisdir = os.path.dirname(__file__)
testdir=os.path.join(thisdir, 'Locateguide/bin/test')

if testdir not in sys.path:
      sys.path.insert(0, testdir)

print testdir
