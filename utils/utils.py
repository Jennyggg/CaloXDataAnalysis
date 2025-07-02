from utils.channel_map import mapDRSChannel2TriggerChannel
import ROOT
def number2string(n):
    s = str(n)
    return s.replace('-', 'm').replace('.', 'p')


def string2number(s):
    return float(s.replace('m', '-').replace('p', '.'))


def getDataFile(runNumber):
    runNum = str(runNumber)
    import json
    jsonFile = "data/datafiles.json"
    with open(jsonFile, 'r') as f:
        data = json.load(f)
    if runNum in data:
        return data[runNum]
    else:
        raise ValueError(f"Run number {runNum} not found in datafiles.json")


def getBranchStats(rdf, branches):
    stats = {
        br: {
            "mean": rdf.Mean(br),
            "min": rdf.Min(br),
            "max": rdf.Max(br)
        } for br in branches
    }
    return stats


def processDRSBoards(rdf):
    import re
    import ROOT
    # Get the list of all branch names
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"DRS.*Group.*Channel.*")
    drs_branches = [b for b in branches if pattern.search(b)]
    stats = getBranchStats(rdf, drs_branches)
    print("DRS branches statistics:")
    for br, res in stats.items():
        print(f"{br}: mean = {res['mean'].GetValue():.4f}, "
              f"min = {res['min'].GetValue():.4f}, "
              f"max = {res['max'].GetValue():.4f}")
        stats[br] = {
            "mean": res['mean'].GetValue(),
            "min": res['min'].GetValue(),
            "max": res['max'].GetValue()
        }

    # Create an array of indices for DRS outputs
    rdf = rdf.Define("TS", "FillIndices(1024)")

    # get the mean of DRS outputs per channel
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_median",
            f"compute_median({varname})"
        )
        rdf = rdf.Define(
            f"{varname}_subtractMedian",
            f"{varname} - {varname}_median"
        )

    return rdf


def filterPrefireEvents(rdf, TS=350):
    # use the hodo trigger to filter prefire events
    from utils.channel_map import buildHodoTriggerChannels
    trigger_name_top, trigger_name_bottom = buildHodoTriggerChannels()
    # index of the minimum value in the trigger channels
    rdf = rdf.Define(
        "TS_fired_up", f"ROOT::VecOps::ArgMin({trigger_name_top})")
    rdf = rdf.Define(
        "TS_fired_down", f"ROOT::VecOps::ArgMin({trigger_name_bottom})")

    rdf = rdf.Define(
        "NormalFired", f"(TS_fired_up >= {TS}) && (TS_fired_down >= {TS})")

    rdf_prefilter = rdf
    rdf = rdf.Filter("NormalFired == 1")

    return rdf, rdf_prefilter

def processDRSPeaks(rdf,drs_branches,trigger_channels):
    ROOT.gInterpreter.Declare("""
    #include "ROOT/RVec.hxx"
    #include <algorithm>
    float compute_peakAmp(ROOT::RVec<float> vec) {
        if (vec.empty()) return -9999;
        float peak_amp = *std::max_element(vec.begin(), vec.end());
        std::sort(vec.begin(), vec.end());
        size_t n = vec.size();
        if (n % 2 == 0){
            return peak_amp - (vec[n / 2 - 1] + vec[n / 2])/2;
        }
        else{
            return peak_amp - vec[n / 2];
        }
    }
    """)

    ROOT.gInterpreter.Declare("""
    #include "ROOT/RVec.hxx"
    #include <algorithm>
    size_t compute_peakT(ROOT::RVec<float> vec) {
        if (vec.empty()) return -9999;
        return std::distance(vec.begin(),std::max_element(vec.begin(), vec.end()));
    }
    """)
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_peakT",
            f"compute_peakT({varname})"
        )
        rdf = rdf.Define(
            f"{varname}_peakA",
            f"compute_peakAmp({varname})"
        )
    ROOT.gInterpreter.Declare("""
    #include "ROOT/RVec.hxx"
    #include <algorithm>
    size_t compute_triggerT(ROOT::RVec<float> vec) {
        if (vec.empty()) return -9999;
        size_t minT = std::distance(vec.begin(),std::min_element(vec.begin(), vec.end()));
        size_t maxT = std::distance(vec.begin(),std::max_element(vec.begin(), std::min_element(vec.begin(), vec.end())));
                             
        double half = (*std::min_element(vec.begin(), vec.end()) + *std::max_element(vec.begin(), std::min_element(vec.begin(), vec.end())))/2;
                          
        for(size_t iT = maxT; iT < minT; iT++){
            if( (vec[iT] > half) && (vec[iT+1] <= half) )         
                return (vec[iT] - half) > (half - vec[iT+1]) ? (iT+1) : iT ;      
        }
        return 0;          
    }
    """)

    for trigname in trigger_channels:
        rdf = rdf.Define(
            f"{trigname}_triggerT",
            f"compute_triggerT({trigname})"
        )
    for varname in drs_branches:
        triggername = mapDRSChannel2TriggerChannel(varname)
        print("DRS channel: ", varname, "mapped trigger: ",triggername)
        rdf = rdf.Define(
            f"{varname}_deltaT",
            f"{varname}_peakT - {triggername}_triggerT + 1024"
        )
    return rdf

