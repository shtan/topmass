#! /usr/bin/env python

# fit mbl distribution to parameterized shape

from ROOT import *
import sys
import collections

f = TFile( 'plotsDataMC.root' )

hm = collections.OrderedDict()
hm['161.5'] = f.Get( 'masses/mbl_ttbar_mass161' )
hm['163.5'] = f.Get( 'masses/mbl_ttbar_mass163' )
hm['166.5'] = f.Get( 'masses/mbl_ttbar_mass166' )
hm['169.5'] = f.Get( 'masses/mbl_ttbar_mass169' )
hm['172.5'] = f.Get( 'masses/mbl_ttbar' )
hm['175.5'] = f.Get( 'masses/mbl_ttbar_mass175' )
hm['178.5'] = f.Get( 'masses/mbl_ttbar_mass178' )
hm['181.5'] = f.Get( 'masses/mbl_ttbar_mass181' )

callmasses = TCanvas('callmasses','all masses',800,800)
callmasses.cd()

lallmasses = TLegend(0.65,0.5,0.85,0.8)
lallmasses.SetFillStyle(0)
lallmasses.SetBorderSize(0)

numparams = 10
func = TF1('func','gaus(0)+exp(-[3]*x)*pol5(4)',0,250)
func.SetParameter(0,-1.33E03 )
func.SetParameter(1, 1.51E02 )
func.SetParameter(2, 2.62E01 )
func.SetParameter(3, 2.05E-02)
func.SetParameter(4,-1.37E01 )
func.SetParameter(5, 8.19E00 )
func.SetParameter(6,-1.64E00 )
func.SetParameter(7, 1.16E-01)
func.SetParameter(8,-9.20E-04)
func.SetParameter(9, 1.90E-06)

gpar = []
for i in range(numparams):
   gpar.append(TGraph())
   gpar[i].SetTitle('parameter '+str(i))
   gpar[i].GetXaxis().SetTitle('Top Mass (GeV)')
   gpar[i].GetYaxis().SetTitle('par value')
   gpar[i].SetMarkerStyle(20)


count = 0
hm['161.5'].Draw('HIST')
for name, hist in hm.iteritems():
   hist.SetLineColor(count+1)
   hist.Scale( hm['172.5'].Integral('width')/hist.Integral('width') )
   hist.Draw('same HIST')
   lallmasses.AddEntry( hist, name, 'l' )

   hist.Fit(func, 'EM')
   for i in range(numparams):
      gpar[i].SetPoint(count, float(name), func.GetParameter(i))

   count=count+1


#callmasses.Draw()
#lallmasses.Draw( 'same' )

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


