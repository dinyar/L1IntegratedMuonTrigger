from ROOT import *
import os

"""
_plotDccRes_

Plots comparation between old and new residuals in phi and phi-bending

"""
max_qual = 16

def names( plot, base, wheel, sector, station, quality ):
    """
    _names_

    Compose output names for plots

    """

    phi = {}
    localdir = "%s/Wheel%d/Sector%d/Station%d/%s/quality_%d/%s"%( base, wheel,
                                                                  sector,
                                                                  station,
                                                                  plot,
                                                                  quality,
                                                                  plot )

    newRes = "%s_new_W%d_Sec%d_St%d_quality_%d"%( localdir, wheel, sector,
                                                  station, quality )
    oldRes = "%s_old_W%d_Sec%d_St%d_quality_%d"%( localdir, wheel, sector,
                                                  station, quality )

    name = "%s_W%d_Sec%d_St%d_quality_%d"%( plot, wheel, sector,
                                            station, quality )
    phi["new"] = newRes
    phi["old"] = oldRes
    phi["name"] = name
    return phi

def phiNames( base, wheel, sector, station, quality ):
    """
    _phiNames_

    Compose phi output name

    """

    return names( "DCC_PhiResidual", base, wheel, sector, station, quality )



def phibNames( base, wheel, sector, station, quality ):
    """
    _phiNames_

    Compose phi-bending output name

    """

    return names( "DCC_PhibResidual", base, wheel, sector, station, quality )



def draw( old, new, name, outdir, write=False ):
    """
    _draw_

    build canvas with all cosmetics

    """

    ### build canvas & cosmetics
    canvas = TCanvas( 'Legacy vs Phase1 '+name,
                      'Legacy vs Phase1 '+name,
                      200, 10, 700, 500 )
    canvas.cd()
    old.SetLineColor( ROOT.kBlue )
    new.SetLineColor( ROOT.kRed )
    #old.SetLineWidth( 2 )
    #new.SetLineWidth( 2 )
    old.SetTitle(name)

    ### coordinates for legend and stats
    X1NDC = .77
    X2NDC = .98
    Y1NDC = .80
    Y2NDC = .90

    ### Legend
    leg = TLegend( X1NDC, Y1NDC, X2NDC, Y2NDC );
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.AddEntry( new, "new", "pl" );
    leg.AddEntry( old, "old", "pl" );

    ### stats for new: need a Draw to make stats accessible
    new.Draw()
    canvas.Update()
    sb3 = new.FindObject("stats")
    Y2NDC = Y1NDC
    Y1NDC = Y2NDC - .12
    sb3.SetX1NDC( X1NDC )
    sb3.SetX2NDC( X2NDC )
    sb3.SetY1NDC( Y1NDC )
    sb3.SetY2NDC( Y2NDC )

    ### stats for old: need a Draw to make stats accessible
    old.Draw()
    canvas.Update()
    sb2 = old.FindObject("stats")
    Y2NDC = Y1NDC
    Y1NDC = Y2NDC - .12
    sb2.SetX1NDC( X1NDC )
    sb2.SetX2NDC( X2NDC )
    sb2.SetY1NDC( Y1NDC )
    sb2.SetY2NDC( Y2NDC )

    ### Draw everything
    old.Draw()
    new.Draw("same")
    sb3.Draw("same")
    sb2.Draw("same")
    leg.Draw()

    canvas.SetLogy()
    canvas.Update()
    if write : 
        canvas.SaveAs( outdir + '/' + name + ".png" )

def build( hfile, base, outdir, write=False ):
    """
    _build_

    fills up arrays out of the plot in the rootfiles,
    combine info and call draw()

    """

    oldPhi = []
    newPhi = []
    oldPhib = []
    newPhib = []
    sector = 1
    station = 1
    windex = 0
    for wheel in range( -2, 3 ) :
        oldPhi.append( [] )
        newPhi.append( [] )
        oldPhib.append( [] )
        newPhib.append( [] )
        for quality in range( 0, max_qual ) :
            phi = phiNames( base, wheel, sector, station, quality )
            phib = phibNames( base, wheel, sector, station, quality )
            oldPhi[windex].append( hfile.Get( phi["old"] ) )
            newPhi[windex].append( hfile.Get( phi["new"] ) )
            oldPhib[windex].append( hfile.Get( phib["old"] ) )
            newPhib[windex].append( hfile.Get( phib["new"] ) )
        windex += 1

    oldPhiSum = oldPhi[0][0].Clone()
    newPhiSum = newPhi[0][0].Clone()

    oldPhibSum = oldPhib[0][0].Clone()
    newPhibSum = newPhib[0][0].Clone()

    qualO = [ 0 ] * max_qual
    qualN = [ 0 ] * max_qual

    for wheel in range( 0, windex ) :
        for quality in range( 0, max_qual ) :
            if wheel == 0 and quality == 0 : continue
            oldPhiSum.Add( oldPhi[wheel][quality] )
            newPhiSum.Add( newPhi[wheel][quality] )
            oldPhibSum.Add( oldPhib[wheel][quality] )
            newPhibSum.Add( newPhib[wheel][quality] )
            qualO[quality] = qualO[quality] + oldPhi[wheel][quality].Integral()
            qualN[quality] = qualN[quality] + newPhi[wheel][quality].Integral()

    for quality in range( 0, max_qual ) :
        print quality, qualO[quality], qualN[quality]

    qualOtot = 0
    qualNtot = 0
    for quality in range( 2, max_qual ) :
        qualOtot = qualOtot + qualO[quality]
        qualNtot = qualNtot + qualN[quality]
    print qualOtot, qualNtot

    ## Draw&Save Phi Residuals
    name = "DCC_PhiResidual_Sec%d_St%d"%( sector, station )
    draw( oldPhiSum, newPhiSum, name, outdir, write )
    
    ## Draw&Save Phi Bending Residuals
    name = "DCC_PhibResidual_Sec%d_St%d"%( sector, station )
    draw( oldPhibSum, newPhibSum, name, outdir, write )

############ split Hi and Ho

#4+5+6 13+14+15
#2 HI 3+5+7+9
#3 HO 4+6+8+10

    oldPhiHi = oldPhi[0][2].Clone()
    oldPhiHo = oldPhi[0][3].Clone()
    
    newPhiHi = newPhi[0][3].Clone()
    newPhiHi.Add( newPhi[0][5] )
    newPhiHi.Add( newPhi[0][7] )
    newPhiHi.Add( newPhi[0][9] )
    
    newPhiHo = newPhi[0][4].Clone()
    newPhiHo.Add( newPhi[0][6] )
    newPhiHo.Add( newPhi[0][8] )
    newPhiHo.Add( newPhi[0][10] )

    oldPhiUncHi = oldPhi[0][2].Clone()
    oldPhiUncHi.SetTitle( "Phi HI uncorrelated compare (old 2, new 3)" )
    newPhiUncHi = newPhi[0][3].Clone()
    newPhiUncHi.SetTitle( "Phi HI uncorrelated compare (old 2, new 3)" )
    oldPhiUncHo = oldPhi[0][3].Clone()
    oldPhiUncHo.SetTitle( "Phi HO uncorrelated compare (old 3, new 4)" )
    newPhiUncHo = newPhi[0][4].Clone()
    newPhiUncHo.SetTitle( "Phi HO uncorrelated compare (old 3, new 4)" )
    
    for wheel in range( 1, windex ) :
        oldPhiHi.Add( oldPhi[wheel][2] )
        oldPhiHo.Add( oldPhi[wheel][3] )

        newPhiHi.Add( newPhi[wheel][3] )
        newPhiHi.Add( newPhi[wheel][5] )
        newPhiHi.Add( newPhi[wheel][7] )
        newPhiHi.Add( newPhi[wheel][9] )
        newPhiHo.Add( newPhi[wheel][4] )
        newPhiHo.Add( newPhi[wheel][6] )
        newPhiHo.Add( newPhi[wheel][8] )
        newPhiHo.Add( newPhi[wheel][10] )

        oldPhiUncHi.Add( oldPhi[wheel][2] )
        newPhiUncHi.Add( newPhi[wheel][3] )
        oldPhiUncHo.Add( oldPhi[wheel][3] )
        newPhiUncHo.Add( newPhi[wheel][4] )

    
    ## Draw&Save Phi Residuals for HI
    name = "DCC_PhiResidual_Hi_Sec%d_St%d"%( sector, station )
    draw( oldPhiHi, newPhiHi, name, outdir, write )
    
    ## Draw&Save Phi Residuals for HO
    name = "DCC_PhiResidual_Ho_Sec%d_St%d"%( sector, station )
    draw( oldPhiHo, newPhiHo, name, outdir, write )

    # Draw&Save Phi Residuals for HI
    name = "Phi uncorrelated compare (old 3, new 4) Hi_Sec%d_St%d"%( sector, station )
    draw( oldPhiUncHi, newPhiUncHi, name, outdir, write )
    print name, oldPhiUncHi.GetBinContent(1), newPhiUncHi.GetBinContent(1)
    print name, oldPhiUncHi.GetBinContent(401), newPhiUncHi.GetBinContent(401)

    # Draw&Save Phi Residuals for HO
    name = "Phi uncorrelated compare (old 3, new 4) Ho_Sec%d_St%d"%( sector, station )
    draw( oldPhiUncHo, newPhiUncHo, name, outdir, write )
    print name, oldPhiUncHo.GetBinContent(1), newPhiUncHo.GetBinContent(1)
    print name, oldPhiUncHo.GetBinContent(401), newPhiUncHo.GetBinContent(401)
    
############

    oldPhibHi = oldPhib[0][2].Clone()
    oldPhibHo = oldPhib[0][3].Clone()

    newPhibHi = newPhib[0][3].Clone()
    newPhibHi.Add( newPhib[0][5] )
    newPhibHi.Add( newPhib[0][7] )
    newPhibHi.Add( newPhib[0][9] )
    
    newPhibHo = newPhib[0][4].Clone()
    newPhibHo.Add( newPhib[0][6] )
    newPhibHo.Add( newPhib[0][8] )
    newPhibHo.Add( newPhib[0][10] )
    

    oldPhibUncHi = oldPhib[0][2].Clone()
    oldPhibUncHi.SetTitle( "Phib HI uncorrelated compare (old 2, new 3)" )
    newPhibUncHi = newPhib[0][3].Clone()
    newPhibUncHi.SetTitle( "Phib HI uncorrelated compare (old 2, new 3)" )
    oldPhibUncHo = oldPhib[0][3].Clone()
    oldPhibUncHo.SetTitle( "Phib HO uncorrelated compare (old 3, new 4)" )
    newPhibUncHo = newPhib[0][4].Clone()
    newPhibUncHo.SetTitle( "Phib HO uncorrelated compare (old 3, new 4)" )

    for wheel in range( 1, windex ) :
        oldPhibHi.Add( oldPhib[wheel][2] )
        oldPhibHo.Add( oldPhib[wheel][3] )

        newPhibHi.Add( newPhib[wheel][3] )
        newPhibHi.Add( newPhib[wheel][5] )
        newPhibHi.Add( newPhib[wheel][7] )
        newPhibHi.Add( newPhib[wheel][9] )
        newPhibHo.Add( newPhib[wheel][4] )
        newPhibHo.Add( newPhib[wheel][6] )
        newPhibHo.Add( newPhib[wheel][8] )
        newPhibHo.Add( newPhib[wheel][10] )


        oldPhibUncHi.Add( oldPhib[wheel][2] )
        newPhibUncHi.Add( newPhib[wheel][3] )
        oldPhibUncHo.Add( oldPhib[wheel][3] )
        newPhibUncHo.Add( newPhib[wheel][4] )
    
    # Draw&Save Phi Bending Residuals for HI
    name = "DCC_PhibResidual_Hi_Sec%d_St%d"%( sector, station )
    draw( oldPhibHi, newPhibHi, name, outdir, write )
    
    # Draw&Save Phi Bending Residuals for HO
    name = "DCC_PhibResidual_Ho_Sec%d_St%d"%( sector, station )
    draw( oldPhibHo, newPhibHo, name, outdir, write )


    # Draw&Save Phi Bending Residuals for HI
    name = "Phib uncorrelated compare (old 3, new 4) Hi_Sec%d_St%d"%( sector, station )
    draw( oldPhibUncHi, newPhibUncHi, name, outdir, write )
    print name, oldPhibUncHi.GetBinContent(1), newPhibUncHi.GetBinContent(1)
    print name, oldPhibUncHi.GetBinContent(401), newPhibUncHi.GetBinContent(401)

    # Draw&Save Phi Bending Residuals for HO
    name = "Phib uncorrelated compare (old 3, new 4) Ho_Sec%d_St%d"%( sector, station )
    draw( oldPhibUncHo, newPhibUncHo, name, outdir, write )
    print name, oldPhibUncHo.GetBinContent(1), newPhibUncHo.GetBinContent(1)
    print name, oldPhibUncHo.GetBinContent(401), newPhibUncHo.GetBinContent(401)

    if not write :
        raw_input()


#####################
### MAIN
#####################
def main() :
    """
    __main__

    """
    # write = False
    write = True
    hfile = TFile('/tmp/gcodispo/L1ITMuonBarrelResolutionComparePlots.root', 'read')
    outdir = "dccplots/"
    os.system( "mkdir -p " + outdir );
    gROOT.Reset()

    base = "L1ITMuResolutionPlotter"
    build( hfile, base, outdir, write )

#####################
if __name__ == '__main__':

    main()

