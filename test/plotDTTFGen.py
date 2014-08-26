from ROOT import *


def drawth1(name,rfile1,rfile2):
    canvas=TCanvas(name+' Legacy vs Phase1  ',name+' Legacy vs Phase1  ')
    h1=rfile1.Get('L1ITMuDTTFPlotter/'+name)
    h2=rfile2.Get('L1ITMuDTTFPlotter/'+name)
    h1.SetLineColor(kRed)
    h2.SetLineColor(kBlue)
    h1.GetXaxis().SetRangeUser(-3,3)
    h2.GetXaxis().SetRangeUser(-3,3)
    h2.Draw()
    h1.Draw('same')
    canvas.cd().SetLogy()
    return [canvas,h1,h2]

def drawgraph(canvas,name,rfile,color,opt):
    h1=rfile.Get('L1ITMuDTTFPlotter/'+name+'Num')
    h2=rfile.Get('L1ITMuDTTFPlotter/'+name)
    #print h1,h2
    graph=TGraphAsymmErrors(h1,h2)
    graph.SetLineColor(color)
    canvas.cd().SetGridy()
    canvas.cd().SetGridx()
    graph.GetYaxis().SetRangeUser(0,1)
    graph.DrawClone(opt)

def drawgraphpoker(name,rfile1,rfile2):
    canvas=TCanvas(name+' Legacy vs Phase1  ',name+' Legacy vs Phase1  ')
    bx0=''
    h11=rfile1.Get('L1ITMuDTTFPlotter/'+name+'Num'+bx0)
    h12=rfile1.Get('L1ITMuDTTFPlotter/'+name)
    h21=rfile2.Get('L1ITMuDTTFPlotter/'+name+'Num'+bx0)
    h22=rfile2.Get('L1ITMuDTTFPlotter/'+name)
    bx0='Bx0'
    h110=rfile1.Get('L1ITMuDTTFPlotter/'+name+'Num'+bx0)
    h120=rfile1.Get('L1ITMuDTTFPlotter/'+name)
    h210=rfile2.Get('L1ITMuDTTFPlotter/'+name+'Num'+bx0)
    h220=rfile2.Get('L1ITMuDTTFPlotter/'+name)
    graph1=TGraphAsymmErrors(h11,h12)
    graph1.SetLineColor(kRed)
    graph2=TGraphAsymmErrors(h21,h22)
    graph2.SetLineColor(kBlue)
    graph3=TGraphAsymmErrors(h110,h120)
    graph3.SetLineColor(kMagenta)
    graph4=TGraphAsymmErrors(h210,h220)
    graph4.SetLineColor(kCyan+2)
    canvas.cd().SetGridy()
    canvas.cd().SetGridx()
    graph1.GetYaxis().SetRangeUser(0,1)
    graph2.GetYaxis().SetRangeUser(0,1)
    graph3.GetYaxis().SetRangeUser(0,1)
    graph4.GetYaxis().SetRangeUser(0,1)
    graph1.Draw('ap')
    graph2.Draw('p')
    graph3.Draw('p')
    graph4.Draw('p')
    canvas.Update()
    return [canvas,graph1,graph2,graph3,graph4,h11,h12,h21,h22,h110,h120,h210,h220]
def drawth2projXratio(file1,name1,file2,name2,min=0,max=1, color=kRed):
    canvas=TCanvas(name1+' div  '+name2+' proj X',name1+' div '+name2+' proj X')
    h1=file1.Get('L1ITMuDTTFPlotter/'+name1)
    h2=file2.Get('L1ITMuDTTFPlotter/'+name2)
    h11=h1.ProjectionX()
    h21=h2.ProjectionX()
    canvas.cd().SetGridy()
    canvas.cd().SetGridx()
    graph1=TGraphAsymmErrors(h11,h21)
    graph1.GetYaxis().SetRangeUser(min,max)
    graph1.GetXaxis().SetRangeUser(0,1.2) #only eta positive
    graph1.SetLineColor(color)
    graph1.Draw('ap')
    return [canvas,graph1,h11,h21,h1,h2]

def drawth2ratio(file1,name1,file2,name2,min=0,max=1):
    canvas=TCanvas(name1+' div '+name2,name1+' div '+name2)
    h1=TH2D()
    h1=file1.Get('L1ITMuDTTFPlotter/'+name1)
    h2=file2.Get('L1ITMuDTTFPlotter/'+name2)
    h1.Divide(h2)
    h1.SetName(name1+'DIV'+name2)
    h1.SetTitle(name1+'DIV'+name2)
    h1.SetMaximum(max)
    h1.SetMinimum(min)
    h1.GetXaxis().SetRangeUser(0,1.2)
    h1.Draw('colz')
    return [canvas,h1]

def drawth1ratiographN(file,num,den,min=0,max=1,colors=kRed):
    canvas=TCanvas()
    hn=[]
    hd=[]
    g=[]
    for j,n in enumerate(num):
        hn.append(file.Get('L1ITMuDTTFPlotter/'+n))
        hd.append(file.Get('L1ITMuDTTFPlotter/'+den[j]))
        g.append(TGraphAsymmErrors(hn[-1],hd[-1]))
        g[-1].SetLineColor(colors[j])
        canvas.cd().SetGridx()
        canvas.cd().SetGridy()
        g[-1].GetYaxis().SetRangeUser(min,max)
        if (j==0):
            g[-1].Draw('ap')
        else:
            g[-1].Draw('p')
    canvas.Update()
    return [canvas,g,hn,hd]

def drawgraphpair(name,rfile1,rfile2,bx0):
    canvas=TCanvas(name+' Legacy vs Phase1 '+bx0,name+' Legacy vs Phase1 '+bx0)
    h11=rfile1.Get('L1ITMuDTTFPlotter/'+name+'Num'+bx0)
    h12=rfile1.Get('L1ITMuDTTFPlotter/'+name)
    h21=rfile2.Get('L1ITMuDTTFPlotter/'+name+'Num'+bx0)
    h22=rfile2.Get('L1ITMuDTTFPlotter/'+name)
    #print h11,h12,h21,h22
    graph1=TGraphAsymmErrors(h11,h12)
    graph1.SetLineColor(kRed)
    graph2=TGraphAsymmErrors(h21,h22)
    graph2.SetLineColor(kBlue)
    canvas.cd().SetGridy()
    canvas.cd().SetGridx()
    graph1.GetYaxis().SetRangeUser(0,1)
    graph2.GetYaxis().SetRangeUser(0,1)
    graph1.Draw('ap')
    graph2.Draw('p')
    canvas.Update()
    return [canvas,graph1,graph2,h11,h12,h21,h22]

def drawth2(name, file):
    canvas=TCanvas('c'+name,'c'+name)
    h=file.Get('L1ITMuDTTFPlotter/'+name)
    h.Draw('colz')
    return [canvas, h]

plots=['phiGen','etaGen','ptGen']
phase1_file=TFile('phase1L1ITMuonBarrelDTTFStudies.root','read')

cc=[]
drawall='yes'
if drawall=='yes':
    for p in plots:
        cc.append(drawth1ratiographN(phase1_file,[p+'NumBx0',p+'NumBx0T2'],[p,p],0,1,[kRed,kBlue]))
cc.append(drawth2projXratio(phase1_file,'ptGenVsEtaGenNumBx0T2NotT1',phase1_file,'ptGenVsEtaGen',0,0.05))
cc.append(drawth2projXratio(phase1_file,'ptGenVsEtaGenNumBx0T1NotT2',phase1_file,'ptGenVsEtaGen',0,0.05))
cc.append(drawth2('qualT2vsQualT1',phase1_file))
cc.append(drawth2('ptT2vsPtT1',phase1_file))
cc.append(drawth2('dptT2vsDptT1',phase1_file))
         
for c in cc:
    c[0].cd().Print('dttfplots/'+c[0].GetTitle()+'.png')
raw_input()
