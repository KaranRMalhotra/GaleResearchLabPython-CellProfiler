import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.java
import pathlib
from collections import Counter
import pandas as pd
import numpy as np

cellprofiler_core.preferences.set_headless()
cellprofiler_core.utilities.java.start_java()
### start cellprofiler

Pipeline = cellprofiler_core.pipeline.Pipeline()
Pipeline.load("UpdatedSynapseIdentifier.cppipe")

SynapsePreFiles = list(pathlib.Path('.').absolute().glob('SynapseFiles/*.nd2'))
SynapseFiles = [file.as_uri() for file in SynapsePreFiles]

###adjustablevariables
synaptophysinobjectsize = 0
PSD95objectsize = 0
secondaryexpansionsize = 20
primaryexpansionsize = 0
maxsynaptophysinintensity = 0
maxpsd95intensity = 0

DoubleMaskAreaList = []
ReferenceParentList = []
ParentRingAroundNucleusList = []
ParentCombinedAreasList = []
ParentCombinedNumberList = []
CounterList = []
NucleusSynapseAverage = []
TotalSynapseList = []
AverageSynaptophysinPunctuaList =[]
AveragePSD95PunctuaList = []
namelist = []

minrange = 2
maxrange = 18
### setting for the area to be parsed for colocalization of the staining punctua

while minrange < 15:

    Pipeline.modules()[8].setting(5).set_value(secondaryexpansionsize)
    Pipeline.modules()[10].setting(3).set_value((minrange,maxrange))
    Pipeline.modules()[11].setting(3).set_value((minrange,maxrange))
    Pipeline.modules()[31].file_name_suffix.set_value("SynapseObjectNumber" + str(minrange) + "-" + str(maxrange))
    Pipeline.modules()[32].file_name_suffix.set_value("SynapseAreaSize" + str(minrange) + "-" + str(maxrange))
    Pipeline.modules()[33].file_name_suffix.set_value("PSD95 Unmasked" + str(minrange) + "-" + str(maxrange))
    Pipeline.modules()[34].file_name_suffix.set_value("Synaptophysin Unmasked" + str(minrange) + "-" + str(maxrange))
    Pipeline.modules()[35].file_name_suffix.set_value("PSD95 Masked" + str(minrange) + "-" + str(maxrange))
    Pipeline.modules()[36].file_name_suffix.set_value("Synaptophysin Masked" + str(minrange) + "-" + str(maxrange))
    ### setting cell profiler pipeline features

    Pipeline.read_file_list(SynapseFiles)
    output_measurements = Pipeline.run()

    print(secondaryexpansionsize)

    DoubleMaskAreaList = output_measurements.get_all_measurements("doublemaskedcombine2", "AreaShape_Area")
    ReferenceParentList = output_measurements.get_all_measurements("doublemaskedcombine2", "Parent_CombinedPSD95Synap")
    ParentRingAroundNucleusList = output_measurements.get_all_measurements("doublemaskedcombine2", "Parent_ringaroundnucleus40")
    ParentCombinedAreasList = output_measurements.get_all_measurements("CombinedPSD95Synap", "AreaShape_Area")
    ParentCombinedNumberList = output_measurements.get_all_measurements("CombinedPSD95Synap", "ObjectNumber")
    ParentSynaptophysinList = output_measurements.get_all_measurements("MaskedSynaptophysin", "Parent_ringaroundnucleus40")
    ParentPSD95List = output_measurements.get_all_measurements("MaskedPSD95", "Parent_ringaroundnucleus40")
    ### getting outputs from the cell profiler pipeline

    DoubleMaskAreaList = [list(x) for x in DoubleMaskAreaList]
    ReferenceParentList = [list(x) for x in ReferenceParentList]
    ParentRingAroundNucleusList = [list(x) for x in ParentRingAroundNucleusList]
    ParentCombinedAreasList = [list(x) for x in ParentCombinedAreasList]
    ParentCombinedNumberList = [list(x) for x in ParentCombinedNumberList]
    ###cell profiler outputs a list of coverslips (each one, which is a list of cells)

    print("DoubleMaskAreaList:" + str(DoubleMaskAreaList))
    print("ReferenceParentList:" + str(ReferenceParentList))
    print("ParentRingAroundNucleusList:" + str(ParentRingAroundNucleusList))
    print("ParentCombinedAreasList:" + str(ParentCombinedAreasList))

    arealist = []
    averageareaslist = []

    for parentcombinednumber, parentcombinedarea, referenceparentlist in zip(ParentCombinedNumberList, ParentCombinedAreasList, ReferenceParentList):
        for parentcombinednumber2, parentcombinedarea2 in zip(parentcombinednumber, parentcombinedarea):
            if parentcombinednumber2 not in referenceparentlist:
                print("removed")
                parentcombinedarea.remove(parentcombinedarea2)
    ### double loop for list of lists, removing cells lack any colocalization of both proteins

    for doublemaskarea, parentcombinedarea3, parentringnucleus in zip(DoubleMaskAreaList, ParentCombinedAreasList, ParentRingAroundNucleusList):
        for doublemaskarea2, parentcombinedarea4, parentringnucleus2 in zip(doublemaskarea, parentcombinedarea3, parentringnucleus):
            dividedarea = doublemaskarea2/parentcombinedarea4
            if (dividedarea < 0.50):
                print("failed")
                parentringnucleus.remove(parentringnucleus2)
            else:
                print("passed")
                arealist.append(doublemaskarea2)
        averagearea = np.average(arealist)
        print("the average area for this set is: " + str(averagearea))
        averageareaslist.append(averagearea)
        arealist.clear()
    ### measuring synaptic area, and removing cells that don't have at least 50% colocalization.

    totalaverage = np.average(averageareaslist)
    print("the total average area is: " + str(totalaverage))

    for PSD95_ringlist in ParentPSD95List:
        TotalPSD95Punctua = len(PSD95_ringlist)
        CellCountP = len(Counter(PSD95_ringlist).values())
        AveragePSD95Punctua_PerCell = TotalPSD95Punctua/CellCountP

        print("Psd95_ringlist: " + str(PSD95_ringlist))
        print("PSD_95 Punctua Per Cell List:" + str(Counter(PSD95_ringlist)))
        print("TotalPSD95Punctua: " + str(TotalPSD95Punctua))
        print('AveragePSD95Punctua_PerCell: ' + str(AveragePSD95Punctua_PerCell))

    for synaptophysin_ringlist in ParentSynaptophysinList:
        TotalSynaptophysinPunctua = len(synaptophysin_ringlist)
        CellCountS = len(Counter(synaptophysin_ringlist).values())
        AverageSynaptophysinPunctua_PerCell = TotalSynaptophysinPunctua / CellCountS

        print("Synaptophysin_ringlist: " + str(synaptophysin_ringlist))
        print("Synaptophysin Punctua Per Cell List:" + str(Counter(synaptophysin_ringlist)))
        print("TotalSynaptophysinPunctua:" + str(TotalSynaptophysinPunctua))
        print('AverageSynaptophysinPunctua_PerCell: ' + str(AverageSynaptophysinPunctua_PerCell))

    for synapse_ringlist in ParentRingAroundNucleusList:
        print("Synapse_ringlist:" + str(synapse_ringlist))
        TotalAmountofSynapses = len(synapse_ringlist) * 2
        NucleusSynapseCounts = len(Counter(synapse_ringlist).values())
        print("listofsynapsepernucleus:" + str(len(Counter(synapse_ringlist).keys())))
        print("Sum of values:" + str(2*(sum(Counter(synapse_ringlist).values()))))

        print("TotalAmountofSynapses:" + str(TotalAmountofSynapses))
        print("NucleusSynapseCounts:" + str(NucleusSynapseCounts))

        AverageSynapsePerNucleus = TotalAmountofSynapses/NucleusSynapseCounts
        print("AverageSynapsePerNucleus:" + str(AverageSynapsePerNucleus))
        ###returns number of synapses per ring (one nucleus)

        name = "SynapseData." + (str(minrange) + '-' + str(maxrange))
        print(name)

        namelist.append(name)
        CounterList.append(NucleusSynapseCounts)
        NucleusSynapseAverage.append(AverageSynapsePerNucleus)
        TotalSynapseList.append(TotalAmountofSynapses)
        AverageSynaptophysinPunctuaList.append(AverageSynaptophysinPunctua_PerCell)
        AveragePSD95PunctuaList.append(AveragePSD95Punctua_PerCell)

    minrange = minrange + 1

SynapseDataFrame = {"Range": namelist, "UniqueNucleiWithSynapseCount": CounterList, "AverageSynapsePerNucleus": NucleusSynapseAverage, "TotalAmountofSynapses": TotalSynapseList, 'TotalAmountofSynaptophysinPunctua': AverageSynaptophysinPunctuaList, 'TotalAmountofPSD95Punctua': AveragePSD95PunctuaList}
SynapseDataFrame = pd.DataFrame(SynapseDataFrame)
SynapseDataFrame.to_excel("SynapseDataSheetMaxRanges.xlsx")





