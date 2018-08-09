import ROOT

def plotcomp(TH1, TH2, xtitle, ytitle, th_x, mA, ebeam):
   plot = ROOT.TCanvas("Plot comparison"+th_x,"Plot comparison",1) 
   div_line = 0.35
   pad1 = ROOT.TPad("pad1","pad1",0,div_line,1.0,1.0)
   pad1.SetBottomMargin(0.)
   pad1.SetRightMargin(0.03)
   pad1.Draw()
   pad1.cd()
   TH1.SetTitle("")
   TH1.SetStats(0)
   TH1.SetLineColor(2)
   TH1.Scale(1/TH1.GetEntries())
   TH1.Draw()
   TH2.Scale(1/TH2.GetEntries())
   TH2.Draw("same")
   TH1.GetYaxis().SetTitle(ytitle)
   TH1.GetYaxis().SetTitleSize(20)
   TH1.GetYaxis().SetTitleFont(43)
   TH1.GetYaxis().SetTitleOffset(1.05)
   TH1.GetYaxis().SetLabelFont(43)
   TH1.GetYaxis().SetLabelSize(15)
   TH1.SetMinimum(0.000000011)
   if th_x == "x":
      TH1.SetMaximum(0.1)
   if th_x == "theta":
      TH1.SetMaximum(1.0)
   pad1.SetLogy()
   if th_x == "theta":
      leg_x0 = 0.6
      leg_y0 = 0.6
   if th_x == "x":
      leg_x0 = 0.2
      leg_y0 = 0.3
   leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.15,leg_y0+0.15)
   leg.SetBorderSize(0)
   leg.AddEntry(TH1, TH1.GetName(),"l")
   leg.AddEntry(TH2, TH2.GetName(),"l")
   leg.Draw()
   ROOT.gStyle.SetTextFont(43)
   label = ROOT.TPaveText(0.4,0.9,0.99,0.96,"brNDC")
   label.SetTextSize(20)
   label.AddText("A' mass = "+str(mA)+" GeV, Beam energy = "+str(ebeam)+" GeV")
   label.SetFillColor(0)
   label.SetBorderSize(0)
   label.Draw() 
   plot.cd()
   pad2 = ROOT.TPad("pad2","pad2",0.0, 0.06, 1.0, div_line)
   pad2.SetTopMargin(0.0)
   pad2.SetBottomMargin(0.1)
   pad2.SetRightMargin(0.03)
   pad2.Draw()
   pad2.cd()
   ratio = TH2.Clone("h3")
   if th_x == "theta":
      pad2.SetLogy()
      ratio.SetMinimum(0.11)
      ratio.SetMaximum(99)
   if th_x == "x":
      ratio.SetMinimum(0.05)
      ratio.SetMaximum(1.95)
   ratio.Sumw2()
   ratio.SetStats(0)
   ratio.Divide(TH1)
#   ratio.SetMarkerStyle(2)
   ratio.SetLineColor(1)
   ratio.Draw("HIST")

   ratio.SetLineWidth(2)
   TH1.SetLineWidth(2)
   TH2.SetLineWidth(2)
   ratio.SetTitle("")
   ratio.GetYaxis().SetTitle("Geant/Mg")
   ratio.GetYaxis().CenterTitle()
   ratio.GetYaxis().SetTitleSize(20)
   ratio.GetYaxis().SetTitleFont(43)
   ratio.GetYaxis().SetTitleOffset(1.05)
   ratio.GetYaxis().SetLabelFont(43)
   ratio.GetYaxis().SetLabelSize(15)

   #ratio.GetXaxis().SetTitle(xtitle)
   #ratio.GetXaxis().SetTitleSize(20)
   #ratio.GetXaxis().SetTitleFont(43)
   #ratio.GetXaxis().SetTitleOffset(3.)
   ratio.GetXaxis().SetLabelFont(43)
   ratio.GetXaxis().SetLabelSize(15)
   ratio.GetXaxis().SetTickSize(0.07)
  
   plot.cd()
   xlabel = ROOT.TPaveText(0.5,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(20)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   if th_x == "theta":
      img.WriteImage("newersimtheta_map_"+str(mA)+"_ebeam_"+str(ebeam)+".png")
   if th_x == "x":
      img.WriteImage("newersimx_map_"+str(mA)+"_ebeam_"+str(ebeam)+".png")
   return plot
