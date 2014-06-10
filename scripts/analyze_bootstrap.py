#! /usr/bin/env python

from ROOT import *
import sys

file = TFile( sys.argv[1] )
tree = file.Get('FitResults')

countmt = [0]*8
meanmt = [0]*8
varmt = [0]*8

mcmasses = [161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5]

for i in range( tree.GetEntries() ):
   tree.GetEntry(i)
   iter=-1
   if tree.mcmass == 161.5:
      iter=0
   elif tree.mcmass == 163.5:
      iter=1
   elif tree.mcmass == 166.5:
      iter=2
   elif tree.mcmass == 169.5:
      iter=3
   elif tree.mcmass == 172.5:
      iter=4
   elif tree.mcmass == 175.5:
      iter=5
   elif tree.mcmass == 178.5:
      iter=6
   elif tree.mcmass == 181.5:
      iter=7
   countmt[iter] += 1
   meanmt[iter] += tree.mt
   varmt[iter] += (tree.mt)**2

for iter in range(8):
   if countmt[iter] != 0:
      meanmt[iter] = meanmt[iter]/countmt[iter]
      varmt[iter] = varmt[iter]/countmt[iter] - meanmt[iter]**2

gresults = TGraphErrors()
chi2=0
for i in range(8):
   gresults.SetPoint(i, mcmasses[i], meanmt[i]-mcmasses[i])
   gresults.SetPointError(i, 0.0, sqrt(varmt[i]))
   chi = 0
   if varmt[i] != 0:
      chi = (meanmt[i]-mcmasses[i])/sqrt(varmt[i])
   chi2 += chi*chi
   print str(mcmasses[i])+': '+str(meanmt[i])+' +- '+str(sqrt(varmt[i]))

print 'chi2 = '+str(chi2)
gresults.SetMarkerStyle(20)

fline = TF1('fline','[0]*x',150,200)
fline.SetParameter(0,0.0)
fline.SetLineStyle(7)

canvas = TCanvas('canvas','canvas',800,800)
canvas.cd()

gresults.SetMaximum(1.5)
gresults.SetMinimum(-1.5)
gresults.Draw('AEP')
fline.Draw('same')

raw_input('')
