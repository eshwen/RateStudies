#!/usr/bin/python 
from ROOT import *
from array import *
from math import fabs,sqrt

DIR= '/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias6Feb2017HighPU_ibx0_-1Events/'
FILEIN = 'L1total.root'
FILEOUT = 'L1Plots.root'

FILEIN = DIR+FILEIN
FILEOUT = DIR+FILEOUT

class Plots:
	def __init__(self, name):
		self.name = name
		self.plots = {} 
		
	def addPlot(self, plotName, nbins, xmin, xmax,xlabel = '',ylabel = ''):
		plot = TH1F ('plot_'+plotName+'_'+self.name, 'plot_'+plotName+'_'+self.name, nbins, xmin, xmax)
		plot.Sumw2(True)       
		plot.SetTitle('')
		plot.GetXaxis().SetTitle(xlabel)
		plot.GetYaxis().SetTitle(ylabel)       
		self.plots[plotName] = plot
  
	def add2DPlot (self, plotName, xbins, xmin, xmax,ybins,ymin,ymax,xlabel = '',ylabel = ''):
		plot2d = TH2F ('plot_'+name+'_'+self.name, 'plot_'+name+'_'+self.name, xbins, xmin, xmax, ybins, ymin,ymax)
		plot2d.GetXaxis().SetTitle(xlabel)
		plot2d.GetYaxis().SetTitle(ylabel)
		plot2d.Sumw2(True) 
		plot2d.Draw("colz")
		self.plots[name] = plot2d
		
	def saveToFile(self, tFile):
		tFile.cd()
		for x in self.plots:
			self.plots[x].Write()



	def getPlot (self, name, isPass= False): # isPass: True: for _pass_, False for _tot_, 
		return self.plots[name]



class effPlots:
	def __init__(self, name):
		self.name = name
		self.plots = {} 
		self.effplots = {} ## contains efficiency plots
       
 
	def addPlot(self, plotName, nbins, xmin, xmax):
		plot_pass = TH1F ('plot_pass_'+plotName+'_'+self.name, 'plot_pass_'+plotName+'_'+self.name, nbins, xmin, xmax)
		plot_tot  = TH1F ('plot_tot_'+plotName+'_'+self.name, 'plot_tot_'+plotName+'_'+self.name, nbins, xmin, xmax)
		self.plots[plotName] = [plot_pass, plot_tot]
            ## call sumw2
		for x in self.plots[plotName]:
			x.Sumw2(True)

	def add2DPlot (self, plotName, xbins, xmin, xmax,ybins,ymin,ymax):
		plot_pass = TH2F ('plot_pass_'+plotName+'_'+self.name, 'plot_pass_'+plotName+'_'+self.name, xbins, xmin, xmax,ybins,ymin,ymax)
		plot_tot  = TH2F ('plot_tot_'+plotName+'_'+self.name, 'plot_tot_'+plotName+'_'+self.name, xbins, xmin, xmax,ybins,ymin,ymax)
		self.plots[plotName] = [plot_pass, plot_tot]
            ## call sumw2
		for x in self.plots[plotName]:
			x.Sumw2(True)
            
            
	def makeEffPlot(self, plotName,xlabel,ylabel):
		# hwPass = self.plots[plotName][0].Clone(self.name+'_eff_'+plotName)
		# hwTot = self.plots[plotName][2].Clone(self.name+'_tot_'+plotName)
		# hwPass.Add(self.plots[plotName][1], -1)
		# hwTot.Add(self.plots[plotName][3], -1)
		# hwPass.Divide(hwTot)
		# self.effplots[plotName] = hwPass
		hwPass = self.plots[plotName][0].Clone(self.name+'_pass_'+plotName)
		hwTot = self.plots[plotName][1].Clone(self.name+'_tot_'+plotName)
		grEff = TGraphAsymmErrors()
		grEff.SetName(self.name+'_'+plotName)
		grEff.BayesDivide(hwPass, hwTot)
		grEff.SetMarkerStyle(8)
		grEff.GetXaxis().SetTitle(xlabel)
		grEff.GetYaxis().SetTitle(ylabel)
		self.effplots[plotName] = grEff


	def make2DEffPlot(self, plotName,xlabel,ylabel):
        # hwPass = self.plots[plotName][0].Clone(self.name+'_eff_'+plotName)
        # hwTot = self.plots[plotName][2].Clone(self.name+'_tot_'+plotName)
        # hwPass.Add(self.plots[plotName][1], -1)
        # hwTot.Add(self.plots[plotName][3], -1)
        # hwPass.Divide(hwTot)
        # self.effplots[plotName] = hwPass
		hEff = self.plots[plotName][0].Clone(self.name+'_pass_'+plotName)
		hTot = self.plots[plotName][1].Clone(self.name+'_tot_'+plotName)
		
		hEff.SetName(self.name+'_'+plotName)
		hEff.Divide(hTot)
		
		hEff.GetXaxis().SetTitle(xlabel)
		hEff.GetYaxis().SetTitle(ylabel)
		self.effplots[plotName] = hEff
		


   

            
	def saveToFile(self, tFile):
		tFile.cd()
		for x in self.plots:
			for y in self.plots[x]:
				y.Write()
		for x in self.effplots:
			self.effplots[x].Write()
	    
	    
	def getPlot (self, name, isPass= False): # isPass: True: for _pass_, False for _tot_, 
		if isPass:
			return self.plots[name][0]
		
		if not isPass:
			return self.plots[name][1]
            
	def getEffPlot (self, name): 
		return self.effplots[name]
    






Plots = Plots("Plots")
effPlots = effPlots("effPlots")


Plots.addPlot('DeltaR_taujet',20,0,0.5,'#Delta R(L1#tau,L1jet)','Events')
Plots.addPlot('DeltaR_jetjet',50,0,5,'#Delta R(L1jet_{i},L1jet_{k})','Events')




#### 
fIn   = TFile.Open(FILEIN)
if not fIn.IsZombie(): print ('File '+FILEIN+' opened')
tIn   = fIn.Get('L1Tree/L1Tree')
nEvt  = tIn.GetEntries()


if (nEvt > 100000): nEvt = 100000 

for ev in range(0, nEvt):
    tIn.GetEntry(ev)
    if (ev%10000==0): print 'Entry {0} / {1}'.format(ev,nEvt) 
    #variables

    jetN = tIn.Stage2jetsNumber
    tauN = tIn.Stage2tausNumber
    mass = float(0.)        
#facciamo una bella funzione per questa cosa?
    dR_taujet = []
    dR_jetjet = []
    for iJet in range(0,jetN):

        jetEt = tIn.stage2_jetEt.at(iJet)   
        jetEta = tIn.stage2_jetEta.at(iJet)
        jetPhi = tIn.stage2_jetPhi.at(iJet)
        jet = TLorentzVector() 
        jet.SetPtEtaPhiM(
            jetEt,
            jetEta,
            jetPhi,
            mass
            )
        for iTau in range(0,tauN):
            tauEt = tIn.stage2_tauEt.at(iTau)   
            tauEta = tIn.stage2_tauEta.at(iTau)
            tauPhi = tIn.stage2_tauPhi.at(iTau)
            tau = TLorentzVector() 
            tau.SetPtEtaPhiM(
                tauEt,
                tauEta,
                tauPhi,
                mass
                )
            dR_taujet.append(tau.DeltaR(jet))
	for kJet in range(0,jetN):
		if not(kJet==iJet):
			jet2Et = tIn.stage2_jetEt.at(kJet)   
			jet2Eta = tIn.stage2_jetEta.at(kJet)
			jet2Phi = tIn.stage2_jetPhi.at(kJet)
			jet2 = TLorentzVector() 
			jet2.SetPtEtaPhiM(
				jet2Et,
				jet2Eta,
				jet2Phi,
				mass
				)
			dR_jetjet.append(jet2.DeltaR(jet))
    

    dR_taujet = sorted(dR_taujet)
    dR_jetjet = sorted(dR_jetjet)

   
#conditions
    


#fill plots        
    if len(dR_taujet)>0: Plots.getPlot('DeltaR_taujet',False).Fill(dR_taujet[0])
    if len(dR_jetjet)>0: Plots.getPlot('DeltaR_jetjet',False).Fill(dR_jetjet[0])
    

fIn.Close()





fOut = TFile (FILEOUT, "recreate")

Plots.saveToFile(fOut)
if not fOut.IsZombie(): print ('Plots saved in '+FILEOUT)









