#!/usr/bin/python
from ROOT import *
import math


directory = '/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6245/';

outputfile = 'VBF_620_110_35_PU_Un_vs_uGT'
drawPU = True
drawLS = False
EmuVsUn = False
compareTrg = False
seeds = ['rateL1_VBF_620_90_30_2017_PU.root',
         'rateL1_VBF_620_100_35_2017_PU.root',
         'rateL1_VBF_620_110_35_2017_PU.root',
         'rateL1_VBF_620_115_35_2017_PU.root',
         'rateL1_VBF_620_115_40_2017_PU.root'
]

#seeds = ['rateL1_VBF_620_110_35_2017_PU.root',
#         'Emu_rateL1_VBF_620_110_35_2017_PU_trg.root',
#]

colors = [kBlue, kMagenta, kSpring, kRed, kCyan]

nBunches =1
scale = (nBunches*11245.6)
unit = "Hz"

        
def getTresholds(string,trglabel):
    if EmuVsUn:
        if (not string.startswith('Emu')):
            temp, proc, mjj , lead, sub,year, end= string.split("_")
            unpack = "Un"
            newstring = (proc,mjj,lead,sub,unpack)
        else:
            emu, temp,proc, mjj , lead, sub,year, end= string.split("_")
            newstring = (proc,mjj,lead,sub,emu)
    else:
        temp, proc, mjj , lead, sub,year, end = string.split("_")
        newstring = (proc,mjj,lead,sub)
    if compareTrg and not trglabel:
        unpack = "Un"
        newstring = (proc,mjj,lead,sub,unpack)
    if trglabel:
        ugt = "uGT"
        newstring = (proc,mjj,lead,sub,ugt)
    dash = "_"
    thresholds = dash.join(newstring)
    return thresholds
    
def DrawAll(ratePlots,canvas,ymin,ymax,pu,ls):
    for i,plot in enumerate(ratePlots):
        ratePlots[i].SetMarkerSize(1.5)
        ratePlots[i].SetMarkerStyle(20)
        ratePlots[i].SetMarkerColor(colors[i])
        if i == 0:
            ratePlots[i].Draw("AP")
            if pu: ratePlots[i].GetXaxis().SetTitle("PileUp")
            if ls: ratePlots[i].GetXaxis().SetTitle("LS")
            ratePlots[i].GetYaxis().SetRangeUser(ymin, ymax)
            ratePlots[i].GetYaxis().SetTitle("Rate (nBunches = %d) [%s]" % (nBunches, unit))
#            ratePlots[i].GetYaxis().SetTitleOffset(1.)
        else:
            ratePlots[i].Draw("P")
            
        canvas.Update()

def DrawRate(infile,pu,ls,trg):
    if pu:
        if trg:
            VBF_Pass=infile.Get("VBF_PassTrg")
        else:
            VBF_Pass=infile.Get("VBF_Pass")
        nEventsPass = infile.Get("nEventsPass")
    if ls:
        if trg:
            VBF_Pass=infile.Get("VBF_PassLSTrg")
        else:
            VBF_Pass=infile.Get("VBF_PassLS")
        nEventsPass = infile.Get("nEventsPassLS") 
    

    
    VBF_Pass.Sumw2(True)
    nEventsPass.Sumw2(True)
    VBF_Pass.SetBinContent(0,0)
    nEventsPass.SetBinContent(0,0)
    if ls:VBF_Pass.Rebin(100)
    if ls and not trg:nEventsPass.Rebin(100)
    if pu:VBF_Pass.Rebin(2)
    if pu and not trg:nEventsPass.Rebin(2)
    plot = TGraphAsymmErrors()
    plot.Divide(VBF_Pass,nEventsPass)
    for point in range(0,plot.GetN()):
        plot.GetY()[point] *= scale
        errorlow = plot.GetErrorYlow(point)
        errorhigh = plot.GetErrorYhigh(point)
        plot.SetPointError(point,0,0,errorlow*scale,errorhigh*scale)
        
    return plot

def Frame(gPad,width=2):
    gPad.Update()
    gPad.RedrawAxis()
    l = TLine()
    l.SetLineWidth(width)
    # top
    l.DrawLine(gPad.GetUxmin(), gPad.GetUymax(), gPad.GetUxmax(), gPad.GetUymax())
    # right
    l.DrawLine(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax())
    # bottom
    l.DrawLine(gPad.GetUxmin(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymin())
    # left
    l.DrawLine(gPad.GetUxmin(), gPad.GetUymin(), gPad.GetUxmin(), gPad.GetUymax())
    

def MakePlot(pu,ls,trg):
    ratePlots = []
    inFile = []
    canvas = TCanvas("canvas","canvas",600,600)
    legsize = 0.08*len(seeds)
    if pu:
        coord_leg =[0.11,0.89,0.5,0.9-legsize]
    else:
        coord_leg =[0.5,0.89,0.89,0.9-legsize]
    leg = TLegend(coord_leg[0],coord_leg[1],coord_leg[2],coord_leg[3])
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    for i, seed in enumerate(seeds):
        print seed
        legentry = getTresholds(seed,False)
        if trg: legentrytrg = getTresholds(seed, True)
        inFile.append(TFile.Open(directory+seed))
        ratePlot = TGraphAsymmErrors()
        ratePlot = DrawRate(inFile[i],pu,ls,False)
        ratePlots.append(ratePlot)
        leg.AddEntry(ratePlot, legentry, "p")
        if trg and i==0:
            ratePlot = DrawRate(inFile[i],pu,ls,True)
            ratePlots.append(ratePlot)
            leg.AddEntry(ratePlot, legentrytrg, "p")    
    ymax = max([plot.GetMaximum() for i,plot in enumerate(ratePlots)])
    ymin = min([plot.GetMinimum() for i,plot in enumerate(ratePlots)])    
    TGaxis.SetMaxDigits(2);
    DrawAll(ratePlots,canvas,ymin,ymax,pu,ls)

    
    label = TLatex()
    l47 = TLine(47, ymin, 47, ymax)
    l47.SetLineColor(2)
    l47.SetLineWidth(2)
    l47.Draw()
    l55 = TLine(55, ymin, 55, ymax)
    l55.SetLineColor(2)
    l55.SetLineWidth(2)
    l55.Draw()
    label.SetTextSize(0.03)
    label.SetTextAlign(33)
    label.SetTextFont(42)
    label.SetTextAngle(90)
    label.SetTextColor(2)
    if pu: label.DrawLatex(47+0.5,ymax-(ymax-ymin)/50,"PU = 47")
    if pu: label.DrawLatex(55+0.5,ymax-(ymax-ymin)/50,"PU = 55")
    #label.DrawLatex(60+0.5,ymax-1,"PU = 60")
    leg.Draw("Same")
    Frame(gPad)
    canvas.Update()
    raw_input()
    if(pu):
        print "saved as "+directory+outputfile+"_PU.pdf"
        canvas.SaveAs(directory+outputfile+"_PU.pdf")
    if(ls):
        print "saved as "+directory+outputfile+"_LS.pdf"
        canvas.SaveAs(directory+outputfile+"_LS.pdf")



if drawPU: MakePlot(pu = True, ls = False, trg = compareTrg)
if drawLS: MakePlot(pu = False, ls = True, trg = compareTrg)
