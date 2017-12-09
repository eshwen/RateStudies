#!/usr/bin/python
from ROOT import *
import math
import ctypes
#lib = ctypes.cdll.LoadLibrary('/in/vrtime/mahesh/blue/rnd/software/test/test.so')

directory = '/data_CMS/cms/amendola/Rate2017_VBF/HighPUFill/'

outputfile = 'VBFinclPlots/VBFincl'
isVBFTau = False
drawPU = False
drawLS = False
if not drawPU: drawLS = True 
EmuVsUn = False
compareTrg = False
compareRejTT28 = True
makeRatioPlots = True

seeds = [
    'VBFincl_620_90_30',
    'VBFincl_620_100_35',
    'VBFincl_620_110_35',
    'VBFincl_620_115_35',
    'VBFincl_620_115_40'
]
seedsVBFtau = [
    'VBFtau_450_35_45',
    'VBFtau_450_40_45',
    'VBFtau_450_45_45',
    'VBFtau_450_45_50',
]

if isVBFTau:
    seeds = list(seedsVBFtau)


    
colors = [
    kBlue,
    kMagenta,
    kRed,
    kSpring,
    kCyan
]

if compareRejTT28: colors = [kBlue , kSpring, kRed]

hexcolors = ["#3f35ff",
             "#e01a4f",
             "#00e2df",
             "#f9d800",
             "#2bd11f",
             "#4cb944"
]


nBunches =1
scale = (nBunches*11245.6)
unit = "Hz"

        
def getTresholds(string,trglabel,rejlabel):
    if EmuVsUn:
        if (not string.startswith('Emu')):
            temp, proc, mjj , lead, sub,year, end= string.split("_")
            unpack = "Un"
            newstring = (proc,mjj,lead,sub,unpack)
        else:
            emu, temp,proc, mjj , lead, sub,year, end= string.split("_")
            newstring = (proc,mjj,lead,sub,emu)
    
    else:
        proc, mjj , lead, sub = string.split("_")
        newstring = (proc,mjj,lead,sub)
    if (compareTrg and not trglabel and not rejlabel):
        unpack = "Un"
        newstring = (proc,mjj,lead,sub,unpack)
    if trglabel:
        ugt = "uGT"
        newstring = (proc,mjj,lead,sub,ugt)

    if rejlabel:
        rej = "w/o TT28"
        newstring = (proc,mjj,lead,sub,rej)
    dash = "_"
    thresholds = dash.join(newstring)
    return thresholds
    
def DrawAll(Plots,canvas,ymin,ymax,pu,ls,makeRatio):
    for i,plot in enumerate(Plots):
    
        Plots[i].SetMarkerSize(1.5)
        Plots[i].SetMarkerStyle(20)
        Plots[i].SetMarkerColor(colors[i])
        if makeRatio:   Plots[i].SetMarkerColor(colors[i])
        #color=TColor()
        #Plots[i].SetMarkerColor(color.GetColor(colors[i][0],colors[i][1],colors[i][2]))
        if (i == 0 and not makeRatio) or (i==1 and makeRatio):
            Plots[i].Draw("AP")
            if pu: Plots[i].GetXaxis().SetTitle("PileUp")
            if ls: Plots[i].GetXaxis().SetTitle("LS")
            Plots[i].GetYaxis().SetRangeUser(ymin, ymax)
            Plots[i].GetYaxis().SetTitle("Rate (nBunches = %d) [%s]" % (nBunches, unit))
            if makeRatio:
                Plots[i].GetYaxis().SetTitle("Rate reduction")
                Plots[i].GetYaxis().SetTitleOffset(1.2)
        elif (not makeRatio) or (i>0 and makeRatio):
            Plots[i].Draw("P")
        canvas.Update()

        
def DrawRatio(ratePlots):
    Plots = []
    for i, plot in enumerate(ratePlots):
        plots = TGraphAsymmErrors()
        if i > 0:
            #plot = ratePlots[1].Clone()
            for point in range(0,ratePlots[1].GetN()):
                if ratePlots[0].GetY()[point] > 0:
                    value = ratePlots[i].GetY()[point]
                    value*= 1./ratePlots[0].GetY()[point]
                   
                    errorlow = ratePlots[i].GetErrorYlow(point)/ratePlots[0].GetY()[point]
                    errorhigh = ratePlots[i].GetErrorYhigh(point)/ratePlots[0].GetY()[point]
                    xvalue = ratePlots[i].GetX()[point]
                   
                    plots.SetPoint(point,xvalue,value)
                    plots.SetPointError(point,0,0,errorlow,errorhigh)
        Plots.append(plots)
    return Plots

            
            
        
def DrawRate(infile,pu,ls,trg,reject,reject60,makeRatio):
    if pu:
        if trg:
            VBF_Pass=infile.Get("VBF_PassTrg")
        elif reject:
            VBF_Pass=infile.Get("VBF_Pass_rej")
        elif reject60:
            print infile
            VBF_Pass=infile.Get("VBF_Pass_rej60")
            print VBF_Pass
        else:
            VBF_Pass=infile.Get("VBF_Pass")
        nEventsPass = infile.Get("nEventsPass")
    if ls:
        if trg:
            VBF_Pass=infile.Get("VBF_PassLSTrg")
        elif reject:
            VBF_Pass=infile.Get("VBF_Pass_rejLS")
        elif reject60:
            VBF_Pass=infile.Get("VBF_Pass_rej60LS")
        else:
            VBF_Pass=infile.Get("VBF_PassLS")
        nEventsPass = infile.Get("nEventsPassLS") 
        
    VBF_Pass.Sumw2(True)
    nEventsPass.Sumw2(True)
    VBF_Pass.SetBinContent(0,0)
    nEventsPass.SetBinContent(0,0)
    if ls:VBF_Pass.Rebin(100)
    if ls and not trg and not reject and not reject60: nEventsPass.Rebin(100)

    if pu:VBF_Pass.Rebin(2)
    if pu and not trg and not reject and not reject60: nEventsPass.Rebin(2)
    print "rebin"
    if makeRatio:
        VBF_Pass.Rebin(2)
        if not reject and not reject60: nEventsPass.Rebin(2)

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
    

def MakePlot(pu,ls,trg,reject,makeRatio):
    ratePlots = []
    inFile = []
    canvas = TCanvas("canvas","canvas",600,600)



    for i, seed in enumerate(seeds):
        print seed
        filename = "VBF_ratePU.root"
        if isVBFTau:
            filename = "VBFtau_ratePU.root"
        inFile.append(TFile.Open(directory+seed+"/"+filename))
        ratePlot = TGraphAsymmErrors()
        ratePlot = DrawRate(inFile[i],pu,ls,False,False,False,makeRatio)
        ratePlots.append(ratePlot)
        
        if trg and i==0:
            ratePlot = DrawRate(inFile[i],pu,ls,True,False,False,makeRatio)
            ratePlots.append(ratePlot)
        
        if reject and i==0:
            ratePlot = DrawRate(inFile[i],pu,ls,False,True,False,makeRatio)
            ratePlots.append(ratePlot)
            ratePlot = DrawRate(inFile[i],pu,ls,False,False,True,makeRatio)
            ratePlots.append(ratePlot)
            break

            
    ymax = max([plot.GetMaximum() for i,plot in enumerate(ratePlots)])
    ymin = min([plot.GetMinimum() for i,plot in enumerate(ratePlots)])

    
    TGaxis.SetMaxDigits(2);
    if makeRatio :
        ratioPlots=list(DrawRatio(ratePlots))
        DrawAll(ratioPlots,canvas,0,1.5,pu,ls,makeRatio)
        xmin = ratioPlots[1].GetXaxis().GetXmin()
        xmax = ratioPlots[1].GetXaxis().GetXmax() 
        l1 = TLine(xmin, 1, xmax, 1)
        l1.SetLineColor(colors[0])
        l1.SetLineWidth(2)
        l1.Draw()
    else:
        DrawAll(ratePlots,canvas,ymin,ymax,pu,ls,makeRatio)

    
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
   
    if makeRatio:     Plots = list(ratioPlots)
    else:  Plots = list(ratePlots)


    legysize = 0.04*len(Plots)+0.01*(len(Plots)-1)+0.02
    legxsize = 0.4
    if reject: legxsize = 0.6
    if pu:
        coord_leg =[0.11,0.89,0.11+legxsize,0.9-legysize]
    else:
        coord_leg =[0.89-legxsize,0.89,0.89,0.9-legysize]
    if makeRatio:
        coord_leg =[0.8-legxsize,0.89,0.8,0.9-legysize]
        if not pu:
            coord_leg =[0.11,0.89,0.11+legxsize,0.9-legysize]
    leg = TLegend(coord_leg[0],coord_leg[1],coord_leg[2],coord_leg[3])
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    for i, seed in enumerate(seeds):
        
        legentry = getTresholds(seed,False,False)
        if trg: legentrytrg = getTresholds(seed, True,False)
        if reject: legentryrej = getTresholds(seed, False,True)
        
        
        if makeRatio and not reject:
            if i==0: leg.AddEntry(l1, legentry, "l")
            else:  leg.AddEntry(Plots[i], legentry, "p")
        else:
            if (not reject): leg.AddEntry(Plots[i], legentry, "p")
            if trg and i==0:
        
                leg.AddEntry(Plots[i], legentrytrg, "p")
            if reject and i==0:
                if makeRatio: leg.AddEntry(l1, legentry, "l")
                if not makeRatio: leg.AddEntry(Plots[0], legentry, "p")
                leg.AddEntry(Plots[1], legentryrej, "p")
                legentryrej = legentryrej + ' (jets pT < 60 GeV)'
                leg.AddEntry(Plots[2], legentryrej, "p")
                break



    leg.Draw("Same")
    Frame(gPad)
    canvas.Update()
    raw_input()
    filename = outputfile
    if(makeRatio):  filename=outputfile+"_wrtBaseline"
    if(reject): filename = filename + "_noTT28"
    if(pu):
        print "saved as "+directory+filename+"_PU.pdf"
        canvas.SaveAs(directory+filename+"_PU.pdf")
        canvas.SaveAs(directory+filename+"_PU.png")
    if(ls):
        print "saved as "+directory+filename+"_LS.pdf"
        canvas.SaveAs(directory+filename+"_LS.pdf")
        canvas.SaveAs(directory+filename+"_LS.png")
        



if drawPU: MakePlot(pu = True, ls = False, trg = compareTrg,reject = compareRejTT28, makeRatio=makeRatioPlots)
if drawLS: MakePlot(pu = False, ls = True, trg = compareTrg,reject = compareRejTT28, makeRatio=makeRatioPlots)
