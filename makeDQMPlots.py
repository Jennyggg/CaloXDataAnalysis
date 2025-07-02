import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels
from utils.utils import number2string
from utils.html_generator import generate_html
from utils.validateMap import DrawFERSBoards, DrawDRSBoards
from runNumber import runNumber
from TSwindow import peakTSWindow, deltaTSWindow 
import numpy as np
from DRSPeakTre import threBinsSci,threBinsCer, threColorsSci, threColorsCer, threBoardColorSci, threBoardColorCer

print("Start running script")
ROOT.gROOT.SetBatch(True)

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)
time_reference_channels = buildTimeReferenceChannels(run=runNumber)
hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)
hodo_pos_channels = buildHodoPosChannels(run=runNumber)


rootdir = f"root/Run{runNumber}/"
outdir = f"plots/Run{runNumber}/"


def makeFERS1DPlots():
    plots = []

    infile_name = f"{rootdir}fers_all_channels_1D.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/FERS_1D"
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hist_C_name = f"hist_FERS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}"
            hist_S_name = f"hist_FERS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}"
            hist_C = infile.Get(hist_C_name)
            hist_S = infile.Get(hist_S_name)
            if not hist_C or not hist_S:
                print(
                    f"Warning: Histograms {hist_C_name} or {hist_S_name} not found in {infile_name}")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"Board: {FERSBoard.boardNo}")
            extraToDraw.AddText(f"Tower X: {iTowerX}")
            extraToDraw.AddText(f"Tower Y: {iTowerY}")
            extraToDraw.AddText(
                f"Cer Channel: {FERSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True).channelNo}")
            extraToDraw.AddText(
                f"Sci Channel: {FERSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=False).channelNo}")

            output_name = f"Energy_Board{boardNo}_iTowerX{sTowerX}_iTowerY{sTowerY}"
            DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)

            plots.append(output_name + ".png")

    output_html = f"html/Run{runNumber}/FERS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


# 2D FERS histograms, hg vs lg
def makeFERS2DPlots():
    plots = []
    outdir_plots = outdir + "/FERS_2D"
    infile_name = f"{rootdir}/fers_all_channels_2D.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_hg_vs_lg"
                hist = infile.Get(hist_name)

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {FERSBoard.boardNo}")
                extraToDraw.AddText(f"Tower X: {iTowerX}")
                extraToDraw.AddText(f"Tower Y: {iTowerY}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                output_name = f"FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_hg_vs_lg"
                DrawHistos([hist], f"", 0, 9000, "HG", 0, 1500, "LG",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True, extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"html/Run{runNumber}/FERS_2D/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# FERS output vs event
def trackFERSPlots():
    plots = []
    outdir_plots = outdir + "/FERS_vs_Event"
    infile_name = f"{rootdir}/fers_all_channels_2D_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {FERSBoard.boardNo}")
                extraToDraw.AddText(f"Tower X: {iTowerX}")
                extraToDraw.AddText(f"Tower Y: {iTowerY}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                nEvents = hist.GetXaxis().GetXmax()

                output_name = f"FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist], "", 0, nEvents, "Event", 1, 1e5, f"{var} Energy HG",
                           output_name,
                           dology=True, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"html/Run{runNumber}/FERS_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# 1D histograms for DRS variables
def makeDRS1DPlots():
    # PeakT distributions of DRS peaks passing certain thresholds
    plots = []
    infile_name = f"{rootdir}/drs_all_channels_1D.root"
    infile = ROOT.TFile(infile_name, "READ")
    hist_C_merge = [None]*(len(threBinsCer)-1)
    hist_S_merge = [None]*(len(threBinsSci)-1)
    hist_C_merge_board = {}
    hist_S_merge_board = {}

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        hist_C_merge_board[boardNo] = [None]*(len(threBinsCer)-1)
        hist_S_merge_board[boardNo] = [None]*(len(threBinsSci)-1)
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for ithr in range(len(threBinsCer)-1):
                hist_C_name = f"hist_DRS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}_peakT_amp{threBinsCer[ithr]}-{threBinsCer[ithr+1]}"
                hist_C = infile.Get(hist_C_name)
                if hist_C:
                    if hist_C_merge[ithr]:
                        hist_C_merge[ithr].Add(hist_C)
                    else:
                        hist_C_merge[ithr] = hist_C.Clone(f"hist_DRS_Cer_peakT_merge_amp{threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
                    if hist_C_merge_board[boardNo][ithr]:
                        hist_C_merge_board[boardNo][ithr].Add(hist_C)
                    else:
                        hist_C_merge_board[boardNo][ithr] = hist_C.Clone(f"hist_DRS_board{boardNo}_Cer_peakT_merge_amp{threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
            for ithr in range(len(threBinsSci)-1):
                hist_S_name = f"hist_DRS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}_peakT_amp{threBinsSci[ithr]}-{threBinsSci[ithr+1]}"
                hist_S = infile.Get(hist_S_name)
                if hist_S:
                    if hist_S_merge[ithr]:
                        hist_S_merge[ithr].Add(hist_S)
                    else:
                        hist_S_merge[ithr] = hist_S.Clone(f"hist_DRS_Sci_peakT_merge_amp{threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
                    if hist_S_merge_board[boardNo][ithr]:
                        hist_S_merge_board[boardNo][ithr].Add(hist_S)
                    else:
                        hist_S_merge_board[boardNo][ithr] = hist_S.Clone(f"hist_DRS_board{boardNo}_Sci_peakT_merge_amp{threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
    if None in  hist_C_merge or None in hist_S_merge:
        print(
            f"Warning: Histograms hist_DRS_Cer_peakT_merge or hist_DRS_Sci_peakT_merge has None")
    else:
        bin_centers = (np.arange(peakTSWindow[0],peakTSWindow[1])+np.arange(peakTSWindow[0]+1,peakTSWindow[1]+1))/2
        C_mean = []
        C_std = []
        S_mean = []
        S_std = []
        label_C_merge = []
        label_S_merge = []
        for ithr in range(len(threBinsCer)-1):
            hist_C_merge_window_y = np.array([hist_C_merge[ithr].GetBinContent(i) for i in range(peakTSWindow[0]+1,peakTSWindow[1]+1)])
            C_mean.append(np.average(bin_centers, weights = hist_C_merge_window_y))
            C_std.append(np.sqrt(np.average((bin_centers-C_mean[ithr])**2, weights = hist_C_merge_window_y)))
            hist_C_merge[ithr].Rebin(10)
            label_C_merge.append(f"Cer amp: {threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
        for ithr in range(len(threBinsSci)-1):
            hist_S_merge_window_y = np.array([hist_S_merge[ithr].GetBinContent(i) for i in range(peakTSWindow[0]+1,peakTSWindow[1]+1)])
            S_mean.append(np.average(bin_centers, weights = hist_S_merge_window_y))
            S_std.append(np.sqrt(np.average((bin_centers-S_mean[ithr])**2, weights = hist_S_merge_window_y)))
            hist_S_merge[ithr].Rebin(10)
            label_S_merge.append(f"Sci amp: {threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
        extraToDraw = ROOT.TPaveText(0.15, 0.65, 0.9, 0.9, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Merged peak TS for all broads")
        extraToDraw.AddText(f"TS window: [{peakTSWindow[0]}, {peakTSWindow[1]}]")
        extraToDraw.AddText(f"Sci amp: {threBinsSci}, Cer amp: {threBinsCer}")
        extraToDraw.AddText("Sci std: "+", ".join([f"{S_std[ithr]:.2f}" for ithr in range(len(threBinsSci)-1)]))
        extraToDraw.AddText("Cer std: "+", ".join([f"{C_std[ithr]:.2f}" for ithr in range(len(threBinsCer)-1)]))
        output_name = f"DRS_PeakT_Merged"
        outdir_plots = outdir + "/DRS_1D"
        DrawHistos(hist_S_merge + hist_C_merge, label_S_merge+label_C_merge, peakTSWindow[0], peakTSWindow[1], "peak TS", 0.5, 1e6, "Counts",
            output_name,
            dology=True, drawoptions="HIST", mycolors=threColorsSci+threColorsCer, addOverflow=True, extraToDraw=extraToDraw,
            legendPos=(0.18, 0.5, 0.90, 0.65),
            legendNCols = len(hist_C_merge),
            outdir=outdir_plots)
        plots.append(output_name + ".png")

    output_html = f"html/Run{runNumber}/DRS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


# PeakT distributions of DRS peaks passing certain thresholds
plots = []
infile_name = f"{rootdir}/drs_all_channels_1D.root"
infile = ROOT.TFile(infile_name, "READ")
hist_C_merge = None
hist_S_merge = None
for _, DRSBoard in DRSBoards.items():
    boardNo = DRSBoard.boardNo
    for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
        sTowerX = number2string(iTowerX)
        sTowerY = number2string(iTowerY)

        chan_Cer = DRSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True)
        chan_Sci = DRSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=False)

        hist_C_name = f"hist_DRS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}_peakT"
        hist_S_name = f"hist_DRS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}_peakT"
        hist_C = infile.Get(hist_C_name)
        hist_S = infile.Get(hist_S_name)
        if hist_C:
            if hist_C_merge:
                hist_C_merge.Add(hist_C)
            else:
                hist_C_merge = hist_C.Clone("hist_DRS_Cer_peakT_merge")
        if hist_S:
            if hist_S_merge:
                hist_S_merge.Add(hist_S)
            else:
                hist_S_merge = hist_S.Clone("hist_DRS_Sci_peakT_merge")
if not hist_C_merge or not hist_S_merge:
    print(
        f"Warning: Histograms hist_DRS_Cer_peakT_merge or hist_DRS_Sci_peakT_merge created")
else:
    bin_centers = (np.arange(peakTSWindow[0],peakTSWindow[1])+np.arange(peakTSWindow[0]+1,peakTSWindow[1]+1))/2
    hist_C_merge_window_y = np.array([hist_C_merge.GetBinContent(i) for i in range(peakTSWindow[0]+1,peakTSWindow[1]+1)])
    C_mean = np.average(bin_centers, weights = hist_C_merge_window_y)
    C_std = np.sqrt(np.average((bin_centers-C_mean)**2, weights = hist_C_merge_window_y))
    hist_S_merge_window_y = np.array([hist_S_merge.GetBinContent(i) for i in range(peakTSWindow[0]+1,peakTSWindow[1]+1)])
    S_mean = np.average(bin_centers, weights = hist_S_merge_window_y)
    S_std = np.sqrt(np.average((bin_centers-S_mean)**2, weights = hist_S_merge_window_y))
    hist_C_merge.Rebin(10)
    hist_S_merge.Rebin(10)
    extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Merged peak TS for all broads")
    extraToDraw.AddText(f"TS window: [{peakTSWindow[0]}, {peakTSWindow[1]}]")
    extraToDraw.AddText(f"Sci mean: {S_mean:.2f}, std: {S_std:.2f}")
    extraToDraw.AddText(f"Cer mean: {C_mean:.2f}, std: {C_std:.2f}")
    output_name = f"DRS_PeakT_Merged"
    outdir_plots = outdir + "/DRS_1D"
    DrawHistos([hist_C_merge,hist_S_merge], ["Cer","Sci"], peakTSWindow[0], peakTSWindow[1], "peak TS", 0.5, 1e6, "Counts",
        output_name,
        dology=True, drawoptions="HIST", mycolors=[2,4], addOverflow=True, extraToDraw=extraToDraw,
        legendPos=(0.70, 0.88, 0.90, 0.68),
        outdir=outdir_plots)
    plots.append(output_name + ".png")

    generate_html(plots, outdir_plots,
              output_html=f"html/Run{runNumber}/DRS_peakT/viewer.html")
    
    # DeltaT distributions of DRS peaks passing certain thresholds
    plots = []
    infile_name = f"{rootdir}/drs_all_channels_1D.root"
    infile = ROOT.TFile(infile_name, "READ")
    hist_C_merge = [None]*(len(threBinsCer)-1)
    hist_S_merge = [None]*(len(threBinsSci)-1)
    hist_C_merge_board = {}
    hist_S_merge_board = {}
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        hist_C_merge_board[boardNo] = [None]*(len(threBinsCer)-1)
        hist_S_merge_board[boardNo] = [None]*(len(threBinsSci)-1)
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for ithr in range(len(threBinsCer)-1):
                hist_C_name = f"hist_DRS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}_deltaT_amp{threBinsCer[ithr]}-{threBinsCer[ithr+1]}"
                hist_C = infile.Get(hist_C_name)
                if hist_C:
                    if hist_C_merge[ithr]:
                        hist_C_merge[ithr].Add(hist_C)
                    else:
                        hist_C_merge[ithr] = hist_C.Clone(f"hist_DRS_Cer_deltaT_merge_amp{threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
                    if hist_C_merge_board[boardNo][ithr]:
                        hist_C_merge_board[boardNo][ithr].Add(hist_C)
                    else:
                        hist_C_merge_board[boardNo][ithr] = hist_C.Clone(f"hist_DRS_board{boardNo}_Cer_deltaT_merge_amp{threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
            for ithr in range(len(threBinsSci)-1):
                hist_S_name = f"hist_DRS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}_deltaT_amp{threBinsSci[ithr]}-{threBinsSci[ithr+1]}"
                hist_S = infile.Get(hist_S_name)
                if hist_S:
                    if hist_S_merge[ithr]:
                        hist_S_merge[ithr].Add(hist_S)
                    else:
                        hist_S_merge[ithr] = hist_S.Clone(f"hist_DRS_Sci_deltaT_merge_amp{threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
                    if hist_S_merge_board[boardNo][ithr]:
                        hist_S_merge_board[boardNo][ithr].Add(hist_S)
                    else:
                        hist_S_merge_board[boardNo][ithr] = hist_S.Clone(f"hist_DRS_board{boardNo}_Sci_deltaT_merge_amp{threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
    
    if None in  hist_C_merge or None in hist_S_merge:
        print(
            f"Warning: Histograms hist_DRS_Cer_deltaT_merge or hist_DRS_Sci_deltaT_merge has None")
    else:
        bin_centers = (np.arange(deltaTSWindow[0]+1024,deltaTSWindow[1]+1024)+np.arange(deltaTSWindow[0]+1025,deltaTSWindow[1]+1025))/2
        C_mean = []
        C_std = []
        S_mean = []
        S_std = []
        label_C_merge = []
        label_S_merge = []
        for ithr in range(len(threBinsCer)-1):
            hist_C_merge_window_y = np.array([hist_C_merge[ithr].GetBinContent(i) for i in range(deltaTSWindow[0]+1025,deltaTSWindow[1]+1025)])
            C_mean.append(np.average(bin_centers, weights = hist_C_merge_window_y))
            C_std.append(np.sqrt(np.average((bin_centers-C_mean[ithr])**2, weights = hist_C_merge_window_y)))
            hist_C_merge[ithr].Rebin(10)
            label_C_merge.append(f"Cer amp: {threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
        for ithr in range(len(threBinsSci)-1):
            hist_S_merge_window_y = np.array([hist_S_merge[ithr].GetBinContent(i) for i in range(deltaTSWindow[0]+1025,deltaTSWindow[1]+1025)])
            S_mean.append(np.average(bin_centers, weights = hist_S_merge_window_y))
            S_std.append(np.sqrt(np.average((bin_centers-S_mean[ithr])**2, weights = hist_S_merge_window_y)))
            hist_S_merge[ithr].Rebin(10)
            label_S_merge.append(f"Sci amp: {threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
        extraToDraw = ROOT.TPaveText(0.15, 0.65, 0.9, 0.9, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Merged delta TS for all broads")
        extraToDraw.AddText(f"TS window: [{deltaTSWindow[0]+1024}, {deltaTSWindow[1]+1024}]")
        extraToDraw.AddText(f"Sci amp: {threBinsSci}, Cer amp: {threBinsCer}")
        extraToDraw.AddText("Sci std: "+", ".join([f"{S_std[ithr]:.2f}" for ithr in range(len(threBinsSci)-1)]))
        extraToDraw.AddText("Cer std: "+", ".join([f"{C_std[ithr]:.2f}" for ithr in range(len(threBinsCer)-1)]))
        output_name = f"DRS_DeltaT_Merged"
        outdir_plots = outdir + "/DRS_1D"
        DrawHistos(hist_S_merge + hist_C_merge, label_S_merge+label_C_merge, deltaTSWindow[0]+1024, deltaTSWindow[1]+1024, "peak TS - trigger TS + 1024", 0.5, 1e6, "Counts",
            output_name,
            dology=True, drawoptions="HIST", mycolors=threColorsSci+threColorsCer, addOverflow=True, extraToDraw=extraToDraw,
            legendPos=(0.18, 0.5, 0.90, 0.65),
            legendNCols = len(hist_C_merge),
            outdir=outdir_plots)
        plots.append(output_name + ".png")

    C_mean = []
    C_std = []
    S_mean = []
    S_std = []
    label_C_merge = []
    label_S_merge = []
    bin_centers = (np.arange(deltaTSWindow[0]+1024,deltaTSWindow[1]+1024)+np.arange(deltaTSWindow[0]+1025,deltaTSWindow[1]+1025))/2
    for boardNo in hist_C_merge_board.keys():
        print("hist_C_merge_board",hist_C_merge_board)
        print("hist_C_merge_board[boardNo]",hist_C_merge_board[boardNo])
        if None in hist_C_merge_board[boardNo] or None in hist_S_merge_board[boardNo]:
            print(
            f"Warning: Histograms hist_DRS_board{boardNo}_Cer_deltaT_merge or hist_DRS_board{boardNo}_Sci_deltaT_merge has None")
        else:
            for ithr in range(len(threBinsCer)-1):
                hist_C_merge_window_y = np.array([hist_C_merge_board[boardNo][ithr].GetBinContent(i) for i in range(deltaTSWindow[0]+1025,deltaTSWindow[1]+1025)])
                C_mean.append(np.average(bin_centers, weights = hist_C_merge_window_y))
                C_std.append(np.sqrt(np.average((bin_centers-C_mean[ithr])**2, weights = hist_C_merge_window_y)))
                hist_C_merge_board[boardNo][ithr].Rebin(10)
                label_C_merge.append(f"Board {boardNo}, Cer amp: {threBinsCer[ithr]}-{threBinsCer[ithr+1]}")
            for ithr in range(len(threBinsSci)-1):
                hist_S_merge_window_y = np.array([hist_S_merge_board[boardNo][ithr].GetBinContent(i) for i in range(deltaTSWindow[0]+1025,deltaTSWindow[1]+1025)])
                S_mean.append(np.average(bin_centers, weights = hist_S_merge_window_y))
                S_std.append(np.sqrt(np.average((bin_centers-S_mean[ithr])**2, weights = hist_S_merge_window_y)))
                hist_S_merge_board[boardNo][ithr].Rebin(10)
                label_S_merge.append(f"Board {boardNo}, Sci amp: {threBinsSci[ithr]}-{threBinsSci[ithr+1]}")
    extraToDraw = ROOT.TPaveText(0.15, 0.65, 0.9, 0.9, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Merged delta TS")
    extraToDraw.AddText(f"TS window: [{deltaTSWindow[0]+1024}, {deltaTSWindow[1]+1024}]")
    for ib, boardNo in enumerate(hist_C_merge_board.keys()):
        extraToDraw.AddText(f"Board {boardNo} Sci amp: {threBinsSci}, Cer amp: {threBinsCer}")
        print("ib",ib, "ib*(len(threBinsSci)-1)+ithr",ib*(len(threBinsSci)-1)+ithr)
        print("len(S_std)",len(S_std))
        extraToDraw.AddText("Sci std: "+", ".join([f"{S_std[ib*(len(threBinsSci)-1)+ithr]:.2f}" for ithr in range(len(threBinsSci)-1)]))
        extraToDraw.AddText("Cer std: "+", ".join([f"{C_std[ib*(len(threBinsCer)-1)+ithr]:.2f}" for ithr in range(len(threBinsCer)-1)]))
    output_name = f"DRS_DeltaT_board_Merged"
    outdir_plots = outdir + "/DRS_1D"
    hists_plot = []
    for boardNo in hist_S_merge_board.keys():
        for ithr in range(len(threBinsSci)-1):
            hists_plot.append(hist_S_merge_board[boardNo][ithr])
    for boardNo in hist_C_merge_board.keys():
        for ithr in range(len(threBinsCer)-1):
            hists_plot.append(hist_C_merge_board[boardNo][ithr])
    DrawHistos(hists_plot, label_S_merge+label_C_merge, deltaTSWindow[0]+1024, deltaTSWindow[1]+1024, "peak TS - trigger TS + 1024", 0.5, 1e6, "Counts",
        output_name,
        dology=True, drawoptions="HIST", mycolors=threBoardColorSci+threBoardColorCer, addOverflow=True, extraToDraw=extraToDraw,
        legendPos=(0.18, 0.5, 0.90, 0.65),
        legendNCols = len(hist_C_merge),
        outdir=outdir_plots)
    plots.append(output_name + ".png")

    generate_html(plots, outdir_plots,
                output_html=f"html/Run{runNumber}/DRS_deltaT/viewer.html")

# DRS vs TS
def makeDRS2DPlots(doSubtractMedian=False):
    suffix = ""
    ymin = -50
    ymax = 50
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -20
        ymax = 40
    plots = []
    outdir_plots = outdir + "/DRS_vs_TS"
    infile_name = f"{rootdir}/drs_all_channels_2D.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}{suffix}"
                hist = infile.Get(hist_name)
                output_name = f"DRS_{var}_vs_TS_{sTowerX}_{sTowerY}{suffix}"
                plots.append(output_name + ".png")

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                value_mean = hist.GetMean(2)

                extraToDraw = ROOT.TPaveText(0.20, 0.75, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"B: {DRSBoard.boardNo}, G: {chan.groupNo}, C: {chan.channelNo}")
                extraToDraw.AddText(f"iTowerX: {iTowerX}")
                extraToDraw.AddText(f"iTowerY: {iTowerY}")

                DrawHistos([hist], "", 0, 1024, "Time Slice", value_mean + ymin, value_mean + ymax, f"DRS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
    output_html = f"html/Run{runNumber}/DRS_vs_TS{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# DRS mean vs event
def trackDRSPlots():
    plots = []
    infile_name = f"{rootdir}/drs_all_channels_2D_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/DRS_vs_Event"
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"Board: {DRSBoard.boardNo}")
                extraToDraw.AddText(f"iTowerX: {iTowerX}")
                extraToDraw.AddText(f"iTowerY: {iTowerY}")
                extraToDraw.AddText(f"{var} Group: {chan.groupNo}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                nEvents = hist.GetXaxis().GetXmax()

                output_name = f"DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist], "", 0, nEvents, "Event", 1400, 2300, f"{var} Mean",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"html/Run{runNumber}/DRS_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# time reference
def compareTimeReferencePlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -2500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/time_reference_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/TimeReference"
    for chan_name in time_reference_channels:
        hist_name = f"hist_{chan_name}{suffix}"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"TimeReference_{chan_name}{suffix}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

    output_html = f"html/Run{runNumber}/TimeReference{suffix}/index.html"

    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# trigger
def compareHodoTriggerPlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -1500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/hodo_trigger_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/HodoTrigger"
    for chan_name in hodo_trigger_channels:
        hist_name = f"hist_{chan_name}{suffix}"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"HodoTrigger_{chan_name}{suffix}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

    output_html = f"html/Run{runNumber}/HodoTrigger{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html

# hodo position


def compareHodoPosPlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -1500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/hodo_pos_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/HodoPos"
    for board, channels in hodo_pos_channels.items():
        for chan_name in channels:
            hist_name = f"hist_{chan_name}{suffix}"
            hist = infile.Get(hist_name)
            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue
            extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"{chan_name}")
            output_name = f"HodoPos_{chan_name}{suffix}"
            DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                       extraToDraw=extraToDraw,
                       outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
            plots.append(output_name + ".png")

    output_html = f"html/Run{runNumber}/HodoPos{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


if __name__ == "__main__":
    output_htmls = {}

    # validate DRS and FERS boards
    output_htmls["fers mapping"] = DrawFERSBoards(run=runNumber)
    output_htmls["drs mapping"] = DrawDRSBoards(run=runNumber)

    output_htmls["fers 1D"] = makeFERS1DPlots()
    output_htmls["drs 1D"] = makeDRS1DPlots()
    output_htmls["drs 2D"] = makeDRS2DPlots(doSubtractMedian=True)

    output_htmls["time reference"] = compareTimeReferencePlots(True)
    output_htmls["hodo trigger"] = compareHodoTriggerPlots(True)
    output_htmls["hodo pos"] = compareHodoPosPlots(True)

    print("\n\n\n")
    print("*" * 30)
    for key, value in output_htmls.items():
        print(f"✅ {key} plots can be viewed at: {value}")

    print("All plots generated successfully.")
