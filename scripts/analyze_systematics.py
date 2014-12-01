#! /usr/bin/env python

from ROOT import *
import sys
from operator import itemgetter, attrgetter, methodcaller

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')

mt_central = 0
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   if tree.syst == 'Central':
      mt_central = tree.mt

systarray = []
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   systarray.append( (str(tree.syst), tree.mt-mt_central) )

outformat = '   %-25r: %+3.2f'

print( "\n*** Central Value = %5.2f ***\n" % mt_central )
print( "JES Correlation Groups:\n" )
for syst,val in systarray:
   if "CorrelationGroup" in syst:
      print( outformat % (syst.replace("CorrelationGroup",""), val) )

print ''
for syst,val in systarray:
   if "Total" in syst:
      print( outformat % (syst, val) )

print( "\nOther systematics:\n" )
for syst,val in systarray:
   if not "CorrelationGroup" in syst and not "Central" in syst and not "Total" in syst:
      print( outformat % (syst, val) )

print ''
