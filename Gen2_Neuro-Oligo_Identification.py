import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.java
import pathlib
import numpy as np
import pandas as pd

cellprofiler_core.preferences.set_headless()
cellprofiler_core.utilities.java.start_java()
### start cellprofiler

Pipeline = cellprofiler_core.pipeline.Pipeline()
Pipeline.load("Neuro_Oligo.cppipe")
### load pipline

threshold = 0.105
zikathreshold = 0.26
y = 0
z = 1
neuroninfected = 0.0
oligarchinfected = 0.0
numberofoligarch = 0
numberofneuron = 0
### set values

coverslipnamelist = []
percentneuronslist = []
percentoligarchslist = []
tempneuronslist = []
tempoligarchslist = []
thresholdlist = []
### create lists

S1CS1 = list(pathlib.Path('.').absolute().glob('230609_CS365_HFB_ZIKVtimecourseS1/*.nd2'))
### import all folders for analysis

FS1CS1 = [file.as_uri() for file in S1CS1]
filelist = [FS1CS1]
### create file lists

for i in filelist:
    threshold = 0.0
    print("cycle start")

    while threshold < .40:
        threshold = round(threshold, 2)

        Pipeline.read_file_list(i)
        output_measurements = Pipeline.run()

        possiblenucleilist = output_measurements.get_all_measurements("PossibleNuclei", "Number_Object_Number")
        neuronlist = output_measurements.get_all_measurements("NeuronCell", "Number_Object_Number")
        oligarchlist = output_measurements.get_all_measurements("OligarchCell", "Number_Object_Number")

        print(output_measurements.get_metadata_tags)
        print(output_measurements.get_names)

        for nuclei, neuron in zip(possiblenucleilist, neuronlist):
            percentneuron = len(neuron)/len(nuclei)
            percentneuronslist.append(percentneuron)
            tempneuronslist.append(percentneuron)
        for nuclei, oligarch in zip(possiblenucleilist, oligarchlist):
            percentoligarch = len(oligarch)/len(nuclei)
            percentoligarchslist.append(percentoligarch)
            tempoligarchslist.append(percentoligarch)

        n = len(possiblenucleilist)
        print(n)

        z = 1
        while z <= n:
            if z<=4:
                x = "1DPICS01"
            if 5 <= z <= 8:
                x = "1DPICS02"
            if 9 <= z <= 12:
                x = "2DPICS03"
            if 13 <= z <= 16:
                x = "3DPICS01"
            if 17 <= z <= 20:
                x = "3DPICS02"
            if 21 <= z <= 24:
                x = "4DPICS01"
            if 25 <= z <= 28:
                x = "4DPICS02"
            if z >= 29:
                x = "4DPICS03"

            name = x
            thresholdlist.append(threshold)
            coverslipnamelist.append(name)
            print(name)
            z = z + 1

        print(threshold)
        print(np.average(tempneuronslist))
        print(np.average(tempoligarchslist))

        tempneuronslist.clear()
        tempoligarchslist.clear()
        threshold = threshold + .01

        Pipeline.modules()[8].setting(16).set_value(threshold)
        Pipeline.modules()[9].setting(16).set_value(threshold)
        print("cycle finished")

    dataframe = {"coverslipname": coverslipnamelist, "thresholdlist": thresholdlist, "percentneuron": percentneuronslist,
                 "percentoligo": percentoligarchslist}
    df = pd.DataFrame(dataframe)
    df.to_excel("NeuronOligoCS365df" + ".xlsx")
    print("all done")
    ### create and export dataframe

cellprofiler_core.utilities.java.stop_java()
### terminate the virtual machine and program
