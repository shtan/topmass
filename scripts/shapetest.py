#! /usr/bin/env python

# fit mbl distribution to parameterized shape

from ROOT import *
import sys
import collections

f = TFile( 'plotsDataMC.root' )

hm = collections.OrderedDict()
hm['161.5'] = f.Get( 'all/mbl/hmbl_ttbar161_signal' )
hm['163.5'] = f.Get( 'all/mbl/hmbl_ttbar163_signal' )
hm['166.5'] = f.Get( 'all/mbl/hmbl_ttbar166_signal' )
hm['169.5'] = f.Get( 'all/mbl/hmbl_ttbar169_signal' )
hm['172.5'] = f.Get( 'all/mbl/hmbl_ttbar172_signal' )
hm['175.5'] = f.Get( 'all/mbl/hmbl_ttbar175_signal' )
hm['178.5'] = f.Get( 'all/mbl/hmbl_ttbar178_signal' )
hm['181.5'] = f.Get( 'all/mbl/hmbl_ttbar181_signal' )

callmasses = TCanvas('callmasses','all masses',800,800)
callmasses.cd()

lallmasses = TLegend(0.65,0.5,0.85,0.8)
lallmasses.SetFillStyle(0)
lallmasses.SetBorderSize(0)

numparams = 9
func = TF1('func','gaus(0)+gaus(3)+landau(6)',0,300)
func.SetParameter(0, 9.44076e-03)
func.SetParameter(1, 1.11430e+02)
func.SetParameter(2, 2.02132e+01)
func.SetParameter(3, 8.69241e-03)
func.SetParameter(4, 7.39658e+01)
func.SetParameter(5, 1.97367e+01)
func.SetParameter(6, 1.24551e-02)
func.SetParameter(7, 3.79318e+01)
func.SetParameter(8, 7.60931e+00)

"""
hm['172.5'].Scale( 1.0/hm['172.5'].Integral('width') )
hm['172.5'].Draw()
#hm['172.5'].Fit(func, 'EM')
#func.Draw('same')
callmasses.Draw()
raw_input('')
sys.exit()
"""

gpar = []
for i in range(numparams):
   gpar.append(TGraphErrors())
   gpar[i].SetTitle('parameter '+str(i))
   gpar[i].GetXaxis().SetTitle('Top Mass (GeV)')
   gpar[i].GetYaxis().SetTitle('par value')
   gpar[i].SetMarkerStyle(20)


count = 0
hm['161.5'].Draw('HIST')
for name, hist in hm.iteritems():
   hist.SetLineColor(count+1)
   hist.Scale( 1.0/hist.Integral('width') )
   hist.Draw('same HIST')
   lallmasses.AddEntry( hist, name, 'l' )

   hist.Fit(func, 'EM')
   for i in range(numparams):
      gpar[i].SetPoint(count, float(name), func.GetParameter(i))
      gpar[i].SetPointError(count, 0.0, func.GetParError(i))

   count=count+1


lallmasses.Draw( 'same' )
callmasses.Draw()

fout = TFile( 'output.root', 'RECREATE' )
fout.cd()

fsmooth = TF1('fsmooth','pol1',160,183)

gStyle.SetOptFit(1)
for i in range(numparams):
   ctemp = TCanvas('cpar'+str(i),'cpar'+str(i),800,800)
   ctemp.cd()
   gpar[i].Fit(fsmooth)
   gpar[i].Draw('AP')
   ctemp.Write()
   gpar[i].GetXaxis().SetLimits(160,200)
   ctemp.Print('parfits/'+ctemp.GetName()+'.pdf')

fout.Close()

#hm['172.5'].SetLineColor(1)
#hm['172.5'].Draw('HIST')

#raw_input('')


