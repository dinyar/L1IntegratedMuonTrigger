from ROOT import *
import os

"""
_plotDccRes_

Plots comparation between old and new residuals in phi and phi-bending

"""

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
    _phiNamess_

    Compose phi output name

    """

    return names( "DCC_PhiResidual", base, wheel, sector, station, quality )



def phibNames( base, wheel, sector, station, quality ):
    """
    _phiNamess_

    Compose phi-bending output name

    """

    return names( "DCC_PhibResidual", base, wheel, sector, station, quality )



def draw( old, new, name, outdir, write=False ):
    """
    _draw_

    build canvas with all cosmetics

    """
    canvas = TCanvas( 'Legacy vs Phase1 '+name,
                      'Legacy vs Phase1 '+name,
                      200, 10, 700, 500 )
    canvas.cd()
    old.SetLineColor( ROOT.kBlue )
    new.SetLineColor( ROOT.kRed )
    old.SetTitle(name)

    X1NDC = .77
    X2NDC = .98
    Y1NDC = .80
    Y2NDC = .90

    #leg = TLegend(.7, 0.75, 0.9, 0.9  );
    leg = TLegend( X1NDC, Y1NDC, X2NDC, Y2NDC );
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.AddEntry( new, "new", "pl" );
    leg.AddEntry( old, "old", "pl" );

    old.Draw()
    canvas.Update()
    sb2 = old.FindObject("stats")
    Y2NDC = Y1NDC
    Y1NDC = Y2NDC - .12
    sb2.SetX1NDC( X1NDC )
    sb2.SetX2NDC( X2NDC )
    sb2.SetY1NDC( Y1NDC )
    sb2.SetY2NDC( Y2NDC )

    new.Draw()
    canvas.Update()
    sb3 = new.FindObject("stats")
    Y2NDC = Y1NDC
    Y1NDC = Y2NDC - .12
    sb3.SetX1NDC( X1NDC )
    sb3.SetX2NDC( X2NDC )
    sb3.SetY1NDC( Y1NDC )
    sb3.SetY2NDC( Y2NDC )

    print sb2, sb3
    #old.SetLineWidth( 2 )
    #new.SetLineWidth( 2 )
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
        for quality in range( 0, 7 ) :
            phi = phiNames( base, wheel, sector, station, quality )
            phib = phibNames( base, wheel, sector, station, quality )
            oldPhi[windex].append( hfile.Get( phi["old"] ) )
            newPhi[windex].append( hfile.Get( phi["new"] ) )
            oldPhib[windex].append( hfile.Get( phib["old"] ) )
            newPhib[windex].append( hfile.Get( phib["new"] ) )
        windex += 1
    print windex

    oldPhiSum = oldPhi[0][0].Clone()
    newPhiSum = newPhi[0][0].Clone()

    oldPhibSum = oldPhib[0][0].Clone()
    newPhibSum = newPhib[0][0].Clone()

    for wheel in range( 0, windex ) :
        for quality in range( 0, 7 ) :
            if wheel == 0 and quality == 0 : continue
            oldPhiSum.Add( oldPhi[wheel][quality] )
            newPhiSum.Add( newPhi[wheel][quality] )
            oldPhibSum.Add( oldPhib[wheel][quality] )
            newPhibSum.Add( newPhib[wheel][quality] )

    ## Draw&Save Phi Residuals
    name = "DCC_PhiResidual_Sec%d_St%d"%( sector, station )
    draw( oldPhiSum, newPhiSum, name, outdir, write )
    
    ## Draw&Save Phi Bending Residuals
    name = "DCC_PhibResidual_Sec%d_St%d"%( sector, station )
    draw( oldPhibSum, newPhibSum, name, outdir, write )

############ split Hi and Ho
    oldPhiHi = oldPhi[0][2].Clone()
    oldPhiHo = oldPhi[0][3].Clone()
    
    newPhiHi = newPhi[0][0].Clone()
    newPhiHi.Add( oldPhi[0][2] )
    
    newPhiHo = newPhi[0][1].Clone()
    newPhiHo.Add( newPhi[0][3] )
    
    for wheel in range( 1, windex ) :
        oldPhiHi.Add( oldPhi[wheel][2] )
        oldPhiHo.Add( oldPhi[wheel][3] )
        newPhiHi.Add( newPhi[wheel][0] )
        newPhiHi.Add( newPhi[wheel][2] )
        newPhiHo.Add( newPhi[wheel][1] )
        newPhiHo.Add( newPhi[wheel][3] )

############
    oldPhibHi = oldPhib[0][2].Clone()
    oldPhibHo = oldPhib[0][3].Clone()
    
    newPhibHi = newPhib[0][0].Clone()
    newPhibHi.Add( oldPhib[0][2] )
    
    newPhibHo = newPhib[0][1].Clone()
    newPhibHo.Add( newPhib[0][3] )
    
    for wheel in range( 1, windex ) :
        oldPhibHi.Add( oldPhib[wheel][2] )
        oldPhibHo.Add( oldPhib[wheel][3] )
        newPhibHi.Add( newPhib[wheel][0] )
        newPhibHi.Add( newPhib[wheel][2] )
        newPhibHo.Add( newPhib[wheel][1] )
        newPhibHo.Add( newPhib[wheel][3] )
    
    
    ## Draw&Save Phi Residuals for HI
    name = "DCC_PhiResidual_Hi_Sec%d_St%d"%( sector, station )
    draw( oldPhibHi, newPhibHi, name, outdir, write )

    ## Draw&Save Phi Residuals for HO
    name = "DCC_PhiResidual_Ho_Sec%d_St%d"%( sector, station )
    draw( oldPhibHo, newPhibHo, name, outdir, write )

    # Draw&Save Phi Bending Residuals for HI
    name = "DCC_PhibResidual_Hi_Sec%d_St%d"%( sector, station )
    draw( oldPhibHi, newPhibHi, name, outdir, write )

    # Draw&Save Phi Bending Residuals for HO
    name = "DCC_PhibResidual_Ho_Sec%d_St%d"%( sector, station )
    draw( oldPhibHo, newPhibHo, name, outdir, write )

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
    hfile = TFile('L1ITMuonBarrelResolutionComparePlots.root', 'read')
    outdir = "dccplots/"
    os.system( "mkdir -p " + outdir );
    gROOT.Reset()

    base = "L1ITMuResolutionPlotter"
    build( hfile, base, outdir, write )

#####################
if __name__ == '__main__':

    main()

