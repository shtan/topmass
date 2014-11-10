#! /usr/bin/env python

from ROOT import *
import sys

file = TFile( 'fitresults.root' )
tree = file.Get('FitResults')

mt_central = 0
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   if tree.syst == 'Central':
      mt_central = tree.mt

print '*** Central Value = '+str(mt_central)+' ***'
for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   print tree.syst+': '+str(tree.mt-mt_central)
