#! /usr/bin/env python

from ROOT import *
import sys

# dataset: draw random numbers from unit gaussian

rand = TRandom3(1)

data = []
for n in range(10000):
   data.append( rand.Gaus() )

# sweep over 'mean' parameter
Ls = []
for p in range(100):
   ll = 0
   for n in range(len(data)):
      ll -= log( TMath.Gaus(data[n],-0.025+1.0*p/2000,1.0,true) )

   Ls.append( ll )

# plot
gL = TGraph()
for p in range(len(Ls)):
   gL.SetPoint(p, -0.025+1.0*p/2000, Ls[p]-min(Ls))

gL.SetMarkerStyle(20)

gL.Draw('AP')

raw_input('')
