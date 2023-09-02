import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.java
import pathlib
import sys

import numpy as np
import pandas as pd
import os

cellprofiler_core.preferences.set_headless()
cellprofiler_core.utilities.java.start_java()
### start cellprofiler

Pipeline = cellprofiler_core.pipeline.Pipeline()
Pipeline.load("Astrocyte_NPC Cell.cppipe")
### load pipline

astrozikathreshold = 0.0
npczikathreshold = 0.0
astrothreshold = 0.065
neurothreshold = 0.035
zikathreshold = 0.0
expansion = 0
threshold = 0.105
imagenumber = 1
y = 0
numberofastro = 0
numberofnpc = 0
### set values

coverslipnamelista = []
coverslipnamelistn = []
percentastroslist = []
percentnpcslist = []
percentastrosinfectedlist = []
percentnpcsinfectedlist = []
### create lists

S2CS1 = list(pathlib.Path('.').absolute().glob('230609_CS365_HFB_ZIKVtimecourseS2/*.nd2'))
FS2CS1 = [file.as_uri() for file in S2CS1]
### load in files

filelist = [FS2CS1]

for i in filelist:
    while expansion < 100:
        astrothreshold = 0.100
        npcthreshold = 0.178

        Pipeline.read_file_list(i)
        ### read images

        Pipeline.modules()[9].setting(16).set_value(astrothreshold)
        Pipeline.modules()[10].setting(16).set_value(npcthreshold)
        Pipeline.modules()[11].setting(4).set_value(expansion)
        Pipeline.modules()[12].setting(4).set_value(expansion)
        output_measurements = Pipeline.run()
        ### run pipeline (output is default)

        astrozikalist = output_measurements.get_all_measurements("AstrocyteNuclei2", "Intensity_MeanIntensity_ZIKA")
        npczikalist = output_measurements.get_all_measurements("NPCNuclei2", "Intensity_MeanIntensity_ZIKA")
        ### intensity values

        astrozikalist2 = [list(x) for x in astrozikalist]
        npczikalist2 = [list(x) for x in npczikalist]

        for image in astrozikalist2:
            e = 1
            z = astrozikalist2.index(image) + 1
            if z <= 4:
                x = "1DPICS01"
            if 5 <= z <= 8:
                x = "1DPICS02"
            if 9 <= z <= 12:
                x = "2DPICS01"
            if 13 <= z <= 16:
                x = "2DPICS02"
            if 17 <= z <= 20:
                x = "2DPICS03"
            if 21 <= z <= 24:
                x = "3DPICS01"
            if 25 <= z <= 28:
                x = "3DPICS02"
            if z >= 29:
                x = "4DPICS01"
            while e <= len(image):
                name = str(z) + x + "_cell_" + str(threshold)
                coverslipnamelista.append(name)
                e = e + 1

        for image in npczikalist2:
            o = 1
            k = npczikalist2.index(image) + 1
            if k <= 4:
                x = "1DPICS01"
            if 5 <= k <= 8:
                x = "1DPICS02"
            if 9 <= k <= 12:
                x = "2DPICS01"
            if 13 <= k <= 16:
                x = "2DPICS02"
            if 17 <= k <= 20:
                x = "2DPICS03"
            if 21 <= k <= 24:
                x = "3DPICS01"
            if 25 <= k <= 28:
                x = "3DPICS02"
            if k >= 29:
                x = "4DPICS01"
            while o <= len(image):
                name2 = str(k) + x + "_cell_" + str(threshold)
                coverslipnamelistn.append(name2)
                o = o + 1

        astrozikalist = np.concatenate(astrozikalist).ravel()
        npczikalist = np.concatenate(npczikalist).ravel()

        dataframeastro = {"coverslipnameastros": coverslipnamelista, "astrozikaintensity": astrozikalist}
        dataframenpc = {"coverslipnamenpcs": coverslipnamelistn, "npczikaintensity": npczikalist}
        df = pd.DataFrame(dataframeastro)
        df2 = pd.DataFrame(dataframenpc)

        df.to_excel("astrodfzikavaluesupdatedv2" + str(expansion) + ".xlsx")
        df2.to_excel("npcdfzikavaluesupdatedv2" + str(expansion) + ".xlsx")

        coverslipnamelista.clear()
        coverslipnamelistn.clear()

        expansion = expansion + 5
        print("cyclefinished")

cellprofiler_core.utilities.java.stop_java()
### terminate the virtual machine and program
