#!/usr/bin/env python3

"""
This script codes for a GUI for the DADA2 pipeline.
It provides GUI access to the DADA2 scripts written by
Christoph Schmid.

edited 2018-07-12
"""

import tkinter as tk
import tkinter.filedialog as fd
import tkinter.messagebox
from pubsub import pub
import os
from pathlib import Path
import subprocess as sp
import csv

class sample(object):
    """This class is a data class to hold file references to
    forward and reverse sequence FASTQ files"""
    def __init__(self, name, forward, reverse, currPath):
        """Constructor for Sample Class"""
        self.name=str(name)
        currPath=os.path.normpath(currPath)
        if(os.path.isfile(os.path.join(currPath, forward)) &
               os.path.isfile(os.path.join(currPath, reverse))):
            self.forwardPath=os.path.join(currPath, forward)
            self.reversePath=os.path.join(currPath, reverse)
        else:
            raise FileNotFoundError

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if (isinstance(other, self.__class__)):
            return self.name == other.name
        elif (isinstance(other, str)):
            return self.name == other
        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

class selectFile(tk.Toplevel):

    def __init__(self, version):
        """Constructor Select Files Frame"""
        tk.Toplevel.__init__(self)
        #self.geometry('600x400')
        self.version=version
        self.title('File selection and quality plots (V. ' + self.version[0] + ')')
        self.protocol('WM_DELETE_WINDOW', self.onClose)
        self.sampleList=[]
        self.outDir=""

        self.initUI()

    def initUI(self):
        #HEAD FRAME with button and plot entry
        self.headFrame=tk.Frame(self)
        self.headFrame.pack(side=tk.TOP, fill=tk.X)
        #button for directory selection
        self.dirBtn=tk.Button(self.headFrame, text='Select input folder ...', command=self.chooseDir)
        self.dirBtn.pack(side=tk.LEFT, padx=10, pady=5)
        #label for plot number entry
        self.entryVal=tk.StringVar()
        self.entryVal.set("5")
        self.plotLabel=tk.Label(self.headFrame, text ='How many quality plots to produce?')
        self.plotLabel.pack(side=tk.LEFT, padx=10, pady=5)
        #entry for number of produced plots
        self.plotEntry=tk.Entry(self.headFrame, textvariable=self.entryVal, validate="key",
                             validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.plotEntry.pack(side=tk.LEFT, padx=10, pady=5, fill=tk.BOTH, expand=False)

        #BOTTOM FRAME with list boxes
        self.bottomFrame=tk.Frame(self)
        self.bottomFrame.pack(fill=tk.BOTH, expand=True)
        #list box all
        self.boxAll=tk.Listbox(self.bottomFrame, selectmode=tk.EXTENDED)
        self.boxSel=tk.Listbox(self.bottomFrame, selectmode=tk.EXTENDED)
        self.btnAdd=tk.Button(self.bottomFrame, text='>>', command=self.addSelection)
        self.btnRmv=tk.Button(self.bottomFrame, text='<<', command=self.rmvSelection)

        self.boxAll.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.boxSel.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.btnAdd.pack(side=tk.TOP, padx=5, pady=10)
        self.btnRmv.pack(side=tk.BOTTOM, padx=5, pady=10)

        #RUN FRAME with output directory and execute button
        self.runFrame=tk.Frame(self)
        self.runFrame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.btnOutDir=tk.Button(self.runFrame, text="Select output folder ...",
                                   command=self.chooseOutDir)
        self.outDirLabel=tk.Label(self.runFrame, text=self.outDir)
        self.btnRun=tk.Button(self.runFrame, text="RUN", command=self.runInitScript)

        self.btnOutDir.pack(side=tk.LEFT, padx=10, pady=10)
        self.outDirLabel.pack(side=tk.LEFT)
        self.btnRun.pack(side=tk.RIGHT, padx=10, pady=10)

    def onClose(self):
        """destructor"""
        pub.sendMessage('subWindowClosed')
        self.destroy()

    def onValidate(self, d, S):
        if(int(d) != 1): return True
        try:
            int(S)
        except ValueError:
            self.bell()
            return False
        else:
            return True

    def chooseDir(self):
        self.selDir=fd.askdirectory(mustexist=True)

        # list files in selected directories
        files=[]
        for dirpath, dirnames, filenames in os.walk(self.selDir):
            files.extend(filenames)
            break
        samplesF=sorted(list(filter(lambda x: 'pair1.truncated' in x, files)))
        samplesR=sorted(list(filter(lambda x: 'pair2.truncated' in x, files)))
        sampleNames=[i.split('_')[0] for i in samplesF]
        sampleNames=sorted(list(set(sampleNames)))

        self.boxAll.delete(0, tk.END)

        for name, fw, rv in zip(sampleNames, samplesF, samplesR):
            self.sampleList.append(sample(name, fw, rv, self.selDir))

        for i in self.sampleList:
            self.boxAll.insert(tk.END, i)

    def chooseOutDir(self):
        self.outDir=fd.askdirectory()
        self.outDirLabel.configure(text=self.outDir)

    def addSelection(self):
        if(os.path.isdir(self.selDir)):

            #Test if there is a selection
            selected=sorted(list(self.boxAll.curselection()), reverse=True)
            if not selected: return

            #get sample names from current selection in boxAll
            #add them into the new selection
            newSelection=list()
            for i in selected:
                newSelection.append(self.boxAll.get(i))
            for j in selected:
                self.boxAll.delete(j)

            #get sample names of previously selected items
            prevSelection=self.boxSel.get(first=0, last=tk.END)
            for k in prevSelection:
                newSelection.append(k)

            #sort new selection and add into boxSel
            self.boxSel.delete(first=0, last=tk.END)
            newSelection=sorted(newSelection)
            for l in newSelection:
                self.boxSel.insert(tk.END, l)

        return

    def rmvSelection(self):

        # Test if there is a selection
        selected=sorted(list(self.boxSel.curselection()), reverse=True)
        if not (selected): return

        # get sample names from current selection in boxAll
        # add them into the new selection
        newSelection=list()
        for i in selected:
            newSelection.append(self.boxSel.get(i))
        for j in selected:
            self.boxSel.delete(j)

        # get sample names of entries in boxAll
        prevSelection=self.boxAll.get(first=0, last=tk.END)
        for k in prevSelection:
            newSelection.append(k)

        # sort new selection and add into boxSel
        self.boxAll.delete(first=0, last=tk.END)
        newSelection=sorted(newSelection)
        for l in newSelection:
            self.boxAll.insert(tk.END, l)

        return

    def runInitScript(self):
        #check if all necessary inputs were made
        selected=list(self.boxSel.get(0, tk.END))
        if not selected:
            tk.messagebox.showinfo(title="Data missing",
                                   message="No samples selected!")
            return
        if (self.outDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Output directory missing!")
            return
        if (self.plotEntry.get() == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Number of plots missing!")
            return

        #write input file for R script 'input.R'
        if(not os.path.isdir(self.outDir)): os.mkdir(self.outDir)
        inputFilePath=os.path.join(self.outDir, "inputPaths.txt")
        inputFile=open(inputFilePath, 'w')

        for i in selected:
            for j in self.sampleList:
                if (j == i):
                    print(j.forwardPath, file=inputFile)
                    print(j.reversePath, file=inputFile)
                    break

        inputFile.close()

        #producing command line to run R Script

        process="Rscript"
        pathToScript="/project/genomics/Christoph/DADA2/input.R"
        inputArg=inputFilePath
        outputArg=self.outDir
        plotArg=self.plotEntry.get()

        commandLine=[process, pathToScript,
                     "-i", inputArg,
                     "-o", outputArg,
                     "-p", plotArg,
                     "-V", self.version[0]]

        #run the R script
        try:
            sp.check_call(commandLine, )
        except sp.CalledProcessError:
            tk.messagebox.showerror(title="Error in calling R Script",
                                    message="Execution of R Script input.R failed")
        else:
            tk.messagebox.showinfo(title="input.R",
                                   message="Execution of input.R finished")

class filterReads(tk.Toplevel):

    def __init__(self, version):
        """Constructor Select Files Frame"""
        tk.Toplevel.__init__(self)
        self.geometry('1200x300')
        self.version=version
        self.title('File selection and quality plots (V. ' + self.version[0] + ')')
        self.protocol('WM_DELETE_WINDOW', self.onClose)
        self.forwardReadsPaths=""
        self.reverseReadsPaths=""
        self.outDir=""

        self.initUI()

    def initUI(self):
        # HEAD FRAME with button and plot entry
        self.Frame=tk.Frame(self)
        self.Frame.grid()

        # button for forward and reverse reads file paths selection
        self.truncLLabel=tk.Label(self.Frame, text="cut before:", font="Helvetica 12")
        self.truncRLabel=tk.Label(self.Frame, text="cut after:", font="Helvetica 12 bold")
        self.minLenLabel=tk.Label(self.Frame, text="minimum length:", font="Helvetica 12")
        self.maxLenLabel=tk.Label(self.Frame, text="maximum length:", font="Helvetica 12")
        self.pathLabel=tk.Label(self.Frame, text="Path files in:", font="Helvetica 12 bold")
        self.labelFBtn=tk.Label(self.Frame, text="Forward reads: ", font="Helvetica 12 bold")
        self.labelRBtn=tk.Label(self.Frame, text="Reverse reads: ", font="Helvetica 12 bold")
        self.fileBtnF=tk.Button(self.Frame, text='Select file ...', command=self.chooseForward)
        self.fileBtnR=tk.Button(self.Frame, text='Select file ...', command=self.chooseReverse)
        self.labelFpath=tk.Label(self.Frame, text="Select foward path!")
        self.labelRpath=tk.Label(self.Frame, text="Select reverse path!")

        # Variables for truncation settings
        self.entryValLfwd=tk.StringVar()
        self.entryValLrev=tk.StringVar()
        self.entryValRfwd=tk.StringVar()
        self.entryValRrev=tk.StringVar()
        self.entryValLfwd.set("10")
        self.entryValLrev.set("10")
        self.entryValRfwd.set("300")
        self.entryValRrev.set("300")
        self.entryValminLenF=tk.StringVar()
        self.entryValminLenR=tk.StringVar()
        self.entryValmaxLenF=tk.StringVar()
        self.entryValmaxLenR=tk.StringVar()
        self.entryValminLenF.set("")
        self.entryValminLenR.set("")
        self.entryValmaxLenF.set("")
        self.entryValmaxLenR.set("")


        # entry for truncation of reads
        self.truncEntryLfwd=tk.Entry(self.Frame, textvariable=self.entryValLfwd, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.truncEntryLrev=tk.Entry(self.Frame, textvariable=self.entryValLrev, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.truncEntryRfwd=tk.Entry(self.Frame, textvariable=self.entryValRfwd, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.truncEntryRrev=tk.Entry(self.Frame, textvariable=self.entryValRrev, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.minLenFEntry=tk.Entry(self.Frame, textvariable=self.entryValminLenF, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.minLenREntry=tk.Entry(self.Frame, textvariable=self.entryValminLenR, validate="key",
                                     validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.maxLenFEntry=tk.Entry(self.Frame, textvariable=self.entryValmaxLenF, validate="key",
                                     validatecommand=(self.register(self.onValidate), '%d', '%S'))
        self.maxLenREntry=tk.Entry(self.Frame, textvariable=self.entryValmaxLenR, validate="key",
                                     validatecommand=(self.register(self.onValidate), '%d', '%S'))

        self.truncLLabel.grid(row=0, column=2)
        self.truncRLabel.grid(row=0, column=3)
        self.minLenLabel.grid(row=0, column=4)
        self.maxLenLabel.grid(row=0, column=5)
        self.pathLabel.grid(row=0, column=6, columnspan=5)
        self.labelFBtn.grid(row=1, column=0,padx=10, pady=5, sticky=tk.E)
        self.labelRBtn.grid(row=2, column=0, padx=10, pady=5, sticky=tk.E)
        self.fileBtnF.grid(row=1, column=1, padx=10, pady=5)
        self.fileBtnR.grid(row=2, column=1, padx=10, pady=5)
        self.labelFpath.grid(row=1, column=6, columnspan=5, padx=10)
        self.labelRpath.grid(row=2, column=6, columnspan=5, padx=10)
        self.truncEntryLfwd.grid(row=1, column=2)
        self.truncEntryLrev.grid(row=2, column=2)
        self.truncEntryRfwd.grid(row=1, column=3)
        self.truncEntryRrev.grid(row=2, column=3)
        self.minLenFEntry.grid(row=1, column=4)
        self.minLenREntry.grid(row=2, column=4)
        self.maxLenFEntry.grid(row=1, column=5)
        self.maxLenREntry.grid(row=2, column=5)

        #-----------------
        #separating line 1
        self.line1=tk.Canvas(self.Frame, width=1000, height=10)
        self.line1.grid(row=3, columnspan=11, sticky=tk.W + tk.E)

        self.line1.create_line(0,5,2000,5, fill="black", width=2)

        # maxError and quality settings
        self.maxErrorLabel=tk.Label(self.Frame, text="maximum expected errors:", font="Helvetica 12")
        self.entryValMaxEE=tk.StringVar()
        self.entryValMaxEE.set("")
        self.maxErrorEntry=tk.Entry(self.Frame, textvariable=self.entryValMaxEE, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))

        self.minQualLabel=tk.Label(self.Frame, text="minimum read quality:", font="Helvetica 12 bold")
        self.entryValMinQual=tk.StringVar()
        self.entryValMinQual.set("2")
        self.minQualEntry=tk.Entry(self.Frame, textvariable=self.entryValMinQual, validate="key",
                                  validatecommand=(self.register(self.onValidate), '%d', '%S'))

        self.maxErrorLabel.grid(row=4, column=1, padx=5, pady=10)
        self.maxErrorEntry.grid(row=4, column=2, padx=5, pady=10)
        self.minQualLabel.grid(row=4, column=3, padx=5, pady=10)
        self.minQualEntry.grid(row=4, column=4, padx=5, pady=10)

        #compress, verbose and derep settings

        self.compressVar=tk.IntVar()
        self.compressVar.set(1)
        self.verboseVar=tk.IntVar()
        self.verboseVar.set(1)
        self.derepVar=tk.IntVar()
        self.derepVar.set(0)
        self.compressCB=tk.Checkbutton(self.Frame, text="compress filtered FASTQs", var=self.compressVar)
        self.verboseCB=tk.Checkbutton(self.Frame, text="verbose output", var=self.verboseVar)
        self.derepCB=tk.Checkbutton(self.Frame, text="perform dereplication", var=self.derepVar, state=tk.DISABLED)

        self.compressCB.grid(row=6, column=2)
        self.verboseCB.grid(row=6, column=3)
        self.derepCB.grid(row=6, column=4)


        #-----------------
        #separating line 2
        self.line2=tk.Canvas(self.Frame, width=1000, height=10)
        self.line2.grid(row=7, columnspan=11, sticky=tk.E + tk.W)

        self.line2.create_line(0,5,2000,5, fill="black", width=2)

        #output directory and run button

        #output directory and execute button
        self.btnOutDir=tk.Button(self.Frame, text="Select output folder ...",
                                   command=self.chooseOutDir)
        self.outDirLabel=tk.Label(self.Frame, text=self.outDir)
        self.btnRun=tk.Button(self.Frame, text="RUN", command=self.runFilterScript)

        self.btnOutDir.grid(row=8, column=0, padx=10, pady=10)
        self.outDirLabel.grid(row=8, column= 1, columnspan=7)
        self.btnRun.grid(row=8, column=7, padx=10, pady=10, columnspan=3, sticky=tk.E + tk.W)

        #info about required and optional inputs
        self.infoLabel=tk.Label(self.Frame, text="INFO: required input in BOLD, optional input in STANDARD font",
                                  font="Helvetica 12 bold", fg="red", relief=tk.GROOVE)
        self.infoLabel.grid(row=9, column=2, columnspan=3)

    def chooseOutDir(self):
        self.outDir=fd.askdirectory()
        self.outDirLabel.configure(text=self.outDir)

    def chooseForward(self):
        self.forwardReadsPaths=fd.askopenfilename()
        fwdPath=".../" + "/".join(Path(self.forwardReadsPaths).parts[-3:-1])
        self.labelFpath.configure(text=fwdPath)

    def chooseReverse(self):
        self.reverseReadsPaths=fd.askopenfilename()
        revPath=".../" + "/".join(Path(self.reverseReadsPaths).parts[-3:-1])
        self.labelRpath.configure(text=revPath)

    def onValidate(self, d, S):
        if(int(d) != 1): return True
        try:
            int(S)
        except ValueError:
            self.bell()
            return False
        else:
            return True

    def runFilterScript(self):
        #check if all necessary inputs were made
        if ((self.forwardReadsPaths == "") | (self.reverseReadsPaths == "")):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Path files missing!")
            return

        if (self.outDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Output directory missing!")
            return

        if ((self.truncEntryRfwd.get() == "") | (self.truncEntryRrev.get() == "")):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Truncation settings incomplete!")
            return

        if (self.minQualEntry.get() == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Setting for minimum read quality score missing!")
            return

        #producing command line to run R Script

        process="Rscript"
        pathToScript="/project/genomics/Christoph/DADA2/filtering.R"
        fArg=self.forwardReadsPaths
        rArg=self.reverseReadsPaths
        truncLfwdArg=self.truncEntryLfwd.get() if int(self.truncEntryLfwd.get()) >= 10 else "10"
        truncLrevArg=self.truncEntryLrev.get() if int(self.truncEntryLrev.get()) >= 10 else "10"
        xArg=self.truncEntryRfwd.get()
        yArg=self.truncEntryRrev.get()
        oArg=self.outDir
        minFArg=self.minLenFEntry.get()
        minRArg=self.minLenREntry.get()
        maxFArg=self.maxLenFEntry.get()
        maxRArg=self.maxLenREntry.get()
        eArg=self.maxErrorEntry.get()
        qArg=self.minQualEntry.get()

        commandLine=[process, pathToScript,
                       "-f", fArg,
                       "-r", rArg,
                       "--truncLfwd", truncLfwdArg,
                       "--truncLrev", truncLrevArg,
                       "-x", xArg,
                       "-y", yArg,
                       "-o", oArg,
                       "-q", qArg,
                       "-V", self.version[0]]

        if not (eArg == ""): commandLine.append("-e"), commandLine.append(eArg)
        if not (minFArg == ""): commandLine.append("--minLenF"),  commandLine.append(minFArg)
        if not (minRArg == ""): commandLine.append("--minLenR"),  commandLine.append(minRArg)
        if not (maxFArg == ""): commandLine.append("--maxLenF"),  commandLine.append(maxFArg)
        if not (maxRArg == ""): commandLine.append("--maxLenR"),  commandLine.append(maxRArg)
        if not (self.compressVar.get() == 1): commandLine.append("-c")
        if not (self.verboseVar.get() == 1): commandLine.append("-v")
        if not (self.derepVar.get() == 1): commandLine.append("-d")

        #run the R script
        try:
            sp.check_call(commandLine)
        except sp.CalledProcessError:
            tk.messagebox.showerror(title="Error in calling R Script",
                                    message="Execution of R Script filtering.R failed")
        else:
            tk.messagebox.showinfo(title="filtering.R",
                                   message="Execution of filtering.R finished")

    def onClose(self):
        """destructor"""
        pub.sendMessage('subWindowClosed')
        self.destroy()

class denoiseReads(tk.Toplevel):

    def __init__(self, version):
        """Constructor Select Files Frame"""
        tk.Toplevel.__init__(self)
        #self.geometry('600x400')
        self.version=version
        self.title('File selection and quality plots (V. ' + self.version[0] + ')')
        self.protocol('WM_DELETE_WINDOW', self.onClose)
        self.filtered=""
        self.input=""
        self.outDir=""

        self.initUI()

    def initUI(self):
        # Frame with button and plot entry
        self.Frame=tk.Frame(self)
        self.Frame.grid()

            #file selection buttons
        self.filteredBtn=tk.Button(self.Frame, text="Select filtered fastq files ...",
                                  command=self.selFiltered, font="Helvetica 12")
        self.inputBtn=tk.Button(self.Frame, text="Select input file ...",
                                  command=self.selInput, font="Helvetica 12", state = tk.DISABLED)
        self.outpathBtn=tk.Button(self.Frame, text="Select output folder ...",
                                     command=self.selOutDir, font="Helvetica 12")
        self.labelFiltPath=tk.Label(self.Frame, text="Select directory of filtered fastqs", font="Helvetica 10")
        self.labelInput=tk.Label(self.Frame, text="Select input file", font="Helvetica 10")
        self.labelOutpath=tk.Label(self.Frame, text="Select output directory", font="Helvetica 10")

            #label for plot number entry
        self.plotVar=tk.StringVar()
        self.plotVar.set("5")
        self.plotLabel=tk.Label(self.Frame, text ="How many error model plots to produce?", font="Helvetica 12")
        #entry for number of produced plots
        self.plotEntry=tk.Entry(self.Frame, textvariable=self.plotVar, validate="key",
                             validatecommand=(self.register(self.onValidate), '%d', '%S'))

            #label for pooling and prevalence number entry
        self.poolingVar=tk.StringVar()
        self.poolingVar.set("0")
        self.poolLabel=tk.Label(self.Frame, text ="Pseudo-pooling, prevalence of sequences: (0 = no pooling)",
                                font="Helvetica 12")
        #entry for number of produced plots
        self.poolEntry=tk.Entry(self.Frame, textvariable=self.poolingVar, validate="key",
                             validatecommand=(self.register(self.onValidate), '%d', '%S'),
                                state= tk.DISABLED if self.version[0] < "1.8.0" else tk.NORMAL)

            #Check buttons for binary options
        self.mergeVar=tk.IntVar()
        self.mergeVar.set(0)
        self.mergeCB=tk.Checkbutton(self.Frame, text="override merging", var=self.mergeVar,
                                    font="Helvetica 10", state = tk.DISABLED)

        self.concatVar=tk.IntVar()
        self.concatVar.set(0)
        self.concatCB=tk.Checkbutton(self.Frame, text="concatenate instead of merge", var=self.concatVar,
                                    font="Helvetica 10")
        self.seqtabVar=tk.IntVar()
        self.seqtabVar.set(0)
        self.seqtabCB=tk.Checkbutton(self.Frame, text="override table creation", var=self.seqtabVar,
                                     font="Helvetica 10")
        self.chimeraVar=tk.IntVar()
        self.chimeraVar.set(0)
        self.chimeraCB=tk.Checkbutton(self.Frame, text="override chimera removal", var=self.chimeraVar,
                                      font="Helvetica 10")

            #run button
        self.runBtn=tk.Button(self.Frame, text="RUN", command=self.runDenoiseScript, font = "Helvetica 12")

        #place elements on Frame
            #buttons
        self.filteredBtn.grid(row=1, column=1, columnspan=2, pady=10, padx=5)
        self.inputBtn.grid(row=2, column=1, columnspan=2, pady=10, padx=5)
        self.outpathBtn.grid(row=3, column=1, columnspan=2, pady=10, padx=5)
        self.labelFiltPath.grid(row=1, column=3, padx=5)
        self.labelInput.grid(row=2, column=3, padx=5)
        self.labelOutpath.grid(row=3, column=3, padx=5)

            #plot entry
        self.plotLabel.grid(row=4, column=1, columnspan=2, pady=10, padx=5)
        self.plotEntry.grid(row=4, column=3, pady=10, padx=5)

            #pool entry
        self.poolLabel.grid(row=5, column=1, columnspan=2, pady=10, padx=5)
        self.poolEntry.grid(row=5, column=3, pady=10, padx=5)

            #checkbuttons
        self.mergeCB.grid(row=6, column=1, pady=10, padx=5)
        self.concatCB.grid(row=6, column=2, pady=10, padx=5)
        self.seqtabCB.grid(row=6, column=3, pady=10, padx=5)
        self.chimeraCB.grid(row=6, column=4, pady=10, padx=5)

            #run button
        self.runBtn.grid(row=7, column=3, pady=10, padx=5)

    def selFiltered(self):
        self.filtered=fd.askdirectory()
        self.labelFiltPath.configure(text=self.filtered)

    def selInput(self):
        self.input=fd.askopenfilename()
        self.labelInput.configure(text=self.input)

    def selOutDir(self):
        self.outDir=fd.askdirectory()
        self.labelOutpath.configure(text=self.outDir)

    def runDenoiseScript(self):
        #check if all necessary inputs were made
        if (self.filtered == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Path to filtered fastq files missing!")
            return

        # if (self.input == ""):
            # tk.messagebox.showinfo(title="Data missing",
            #                        message="Input file missing!")
            # return

        if (self.outDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Output directory missing!")
            return

        if (self.plotEntry.get() == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Number of error plots missing!")
            return

        #producing command line to run R Script

        process="qsub"
        pathToScript="/project/genomics/Christoph/DADA2/specialApplications/qsub_inferenceBig_DADA2.sh"
        fArg=self.filtered
        # iArg=self.input
        oArg=self.outDir
        pArg=self.plotVar.get()

        commandLine=[process, pathToScript,
                       "-f", fArg,
                       "-o", oArg,
                       "-p", pArg,
                       "-V", self.version[0],
                       "-r", self.poolingVar.get()
                     ]

        # if (self.mergeVar.get() == 1): commandLine.append("-m")
        if (self.seqtabVar.get() == 1): commandLine.append("-s")
        if (self.chimeraVar.get() == 1): commandLine.append("-c")
        if(self.concatVar.get() == 1): commandLine.append("-t")

        #run the R script
        try:
            sp.check_call(commandLine)
        except sp.CalledProcessError:
            tk.messagebox.showerror(title="Error in calling qsub script",
                                    message="Execution of script qsub_inference_DADA2.sh failed")
        else:
            tk.messagebox.showinfo(title="qsub_inference_DADA2.sh",
                                   message="Upload to server successful")

    def onValidate(self, d, S):
        if(int(d) != 1): return True
        try:
            int(S)
        except ValueError:
            self.bell()
            return False
        else:
            return True

    def onClose(self):
        """destructor"""
        pub.sendMessage('subWindowClosed')
        self.destroy()

class taxonomyReads(tk.Toplevel):

    def __init__(self, version):
        """Constructor Select Files Frame"""
        tk.Toplevel.__init__(self)
        #self.geometry('600x400')
        self.version=version
        self.title('File selection and quality plots (V. ' + self.version[0] + ')')
        self.protocol('WM_DELETE_WINDOW', self.onClose)
        self.input=""
        self.outDir=""

        self.initUI()

    def initUI(self):
        # HEAD FRAME with button and plot entry
        self.Frame=tk.Frame(self)
        self.Frame.grid()

        #file selection buttons
        self.inputBtn=tk.Button(self.Frame, text="Select input file ...",
                                  command=self.selInput, font="Helvetica 12")
        self.outpathBtn=tk.Button(self.Frame, text="Select output folder ...",
                                     command=self.selOutDir, font="Helvetica 12")
        self.labelInput=tk.Label(self.Frame, text="Select input file before continuing", font="Helvetica 10")
        self.labelOutpath=tk.Label(self.Frame, text="Select output directory before continuing", font="Helvetica 10")

        #data base selection
        self.labelDB=tk.Label(self.Frame, text="Database to use:", font="Helvetica 12 underline")

        DATABASES = [
            ("SILVA", "silva"),
            ("RDP", "rdp"),
            ("Green Genes", "gg"),
            ("UNITE (ITS)", "unite")
        ]

        self.dbVar=tk.StringVar()
        self.dbVar.set("silva")

        for idx, (text, db) in enumerate(DATABASES):
            rb=tk.Radiobutton(self.Frame, text=text, variable=self.dbVar,
                              value=db, indicatoron=1)
            rb.grid(row=idx+4, column=1, pady=10, padx=5, sticky=tk.W+tk.E)

        #Check buttons for binary options
        self.biomVar=tk.IntVar()
        self.biomVar.set(0)
        self.biomCB=tk.Checkbutton(self.Frame, text="create BIOM file", var=self.biomVar,
                                    font="Helvetica 10")
        self.treeVar=tk.IntVar()
        self.treeVar.set(0)
        self.treeCB=tk.Checkbutton(self.Frame, text="create phylogeny tree", var=self.treeVar,
                                     font="Helvetica 10")
        self.psVar=tk.IntVar()
        self.psVar.set(1)
        self.psCB=tk.Checkbutton(self.Frame, text="save data as phyloseq", var=self.psVar,
                                      font="Helvetica 10")

        #run button
        self.runBtn=tk.Button(self.Frame, text="RUN", command=self.runTaxonomyScript, font = "Helvetica 12")

        #place elements on Frame
            #buttons
        self.inputBtn.grid(row=1, column=1, pady=10, padx=5)
        self.outpathBtn.grid(row=2, column=1, pady=10, padx=5)
        self.labelInput.grid(row=1, column=2, padx=5)
        self.labelOutpath.grid(row=2, column=2, padx=5)

            #database label
        self.labelDB.grid(row=3, column=1, pady=10, padx=5)

            #checkbuttons
        self.biomCB.grid(row=4, column=2, pady=10, padx=5)
        self.treeCB.grid(row=5, column=2, pady=10, padx=5)
        self.psCB.grid(row=6, column=2, pady=10, padx=5)

            #run button
        self.runBtn.grid(row=7, column=2, pady=10, padx=5)

    def selInput(self):
        self.input=fd.askopenfilename()
        self.labelInput.configure(text=self.input)

    def selOutDir(self):
        self.outDir=fd.askdirectory()
        self.labelOutpath.configure(text=self.outDir)

    def runTaxonomyScript(self):
        #check if all necessary inputs were made
        if (self.input == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Input file missing!")
            return

        if (self.outDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Output directory missing!")
            return

        #producing command line to run R Script

        process="qsub"
        pathToScript="/project/genomics/Christoph/DADA2/qsub_taxonomy_DADA2.sh"
        iArg=self.input
        oArg=self.outDir
        dArg=self.dbVar.get()

        commandLine=[process, pathToScript,
                       "-i", iArg,
                       "-o", oArg,
                       "-d", dArg,
                       "-V", self.version[0]]

        if (self.biomVar.get() == 1): commandLine.append("-b")
        if (self.treeVar.get() == 1): commandLine.append("-t")
        if not (self.psVar.get() == 1): commandLine.append("-p")

        #run the R script
        try:
            sp.check_call(commandLine)
        except sp.CalledProcessError:
            tk.messagebox.showerror(title="Error in calling qsub script",
                                    message="Execution of script qsub_taxonomy_DADA2.sh failed")
        else:
            tk.messagebox.showinfo(title="qsub_taxonomy_DADA2.sh",
                                   message="Upload to server successful")

    def onClose(self):
        """destructor"""
        pub.sendMessage('subWindowClosed')
        self.destroy()

class phyloTree(tk.Toplevel):
    def __init__(self):
        """Constructor Select Files Frame"""
        tk.Toplevel.__init__(self)
        #self.geometry('600x400')
        self.title('Calculation of phylogenetic tree')
        self.protocol('WM_DELETE_WINDOW', self.onClose)
        self.input=""
        self.outDir=""

        self.initUI()

    def initUI(self):
        # HEAD FRAME with button and plot entry
        self.Frame=tk.Frame(self)
        self.Frame.grid()

        #file selection buttons
        self.inputBtn=tk.Button(self.Frame, text="Select input FASTA file ...",
                                  command=self.selInput, font="Helvetica 12")
        self.outpathBtn=tk.Button(self.Frame, text="Select output folder ...",
                                     command=self.selOutDir, font="Helvetica 12")
        self.labelInput=tk.Label(self.Frame, text="Select input FASTA before continuing", font="Helvetica 10")
        self.labelOutpath=tk.Label(self.Frame, text="Select output directory before continuing", font="Helvetica 10")

        #run button
        self.runBtn=tk.Button(self.Frame, text="RUN", command=self.runPhylotreeScript, font = "Helvetica 12")

        #place elements on Frame
            #buttons
        self.inputBtn.grid(row=1, column=1, pady=10, padx=5)
        self.outpathBtn.grid(row=2, column=1, pady=10, padx=5)
        self.labelInput.grid(row=1, column=2, padx=5)
        self.labelOutpath.grid(row=2, column=2, padx=5)

            #run button
        self.runBtn.grid(row=3, column=2, pady=10, padx=5)

    def selInput(self):
        self.input=fd.askopenfilename()
        self.labelInput.configure(text=self.input)

    def selOutDir(self):
        self.outDir=fd.askdirectory()
        self.labelOutpath.configure(text=self.outDir)

    def runPhylotreeScript(self):
        #check if all necessary inputs were made
        if (self.input == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Input FASTA missing!")
            return

        if (self.outDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Output directory missing!")
            return

        #producing command line to run R Script

        process="qsub"
        pathToScript="/project/genomics/Christoph/DADA2/qsub_phylotree_DADA2.sh"
        iArg=self.input
        oArg=self.outDir

        commandLine=[process, pathToScript,
                       "-i", iArg,
                       "-o", oArg]

        #run the R script
        try:
            sp.check_call(commandLine)
        except sp.CalledProcessError:
            tk.messagebox.showerror(title="Error in calling qsub script",
                                    message="Execution of script qsub_phylotree.sh failed")
        else:
            tk.messagebox.showinfo(title="qsub_phylotree.sh",
                                   message="Upload to server successful")

    def onClose(self):
        """destructor"""
        pub.sendMessage('subWindowClosed')
        self.destroy()

class sequenceTracker(tk.Toplevel):

    def __init__(self, version):
        """Constructor Select Files Frame"""
        tk.Toplevel.__init__(self)
        #self.geometry('600x400')
        self.version=version
        self.title('File selection and quality plots (V. ' + self.version[0] + ')')
        self.protocol('WM_DELETE_WINDOW', self.onClose)
        self.filteringDir=""
        self.dadaDir=""
        self.mergeDir=""
        self.seqtabDir=""
        self.chimeraDir=""
        self.outDir=""

        self.initUI()

    def initUI(self):
        """Initialise Elements for Subwindow"""

        # HEAD FRAME with button and plot entry
        self.Frame=tk.Frame(self)
        self.Frame.grid()

        #file selection buttons
        self.filteringBtn=tk.Button(self.Frame, text="Select filtering.RData ...",
                                  command=self.selFiltering, font="Helvetica 12")
        self.dadaBtn=tk.Button(self.Frame, text="Select dada.RData ...",
                                     command=self.selDada, font="Helvetica 12")
        self.mergeBtn=tk.Button(self.Frame, text="Select mergedReads.RData ...",
                                     command=self.selMerge, font="Helvetica 12")
        self.seqtabBtn=tk.Button(self.Frame, text="Select seqTabRaw.RData ...",
                                     command=self.selSeqtab, font="Helvetica 12")
        self.chimeraBtn=tk.Button(self.Frame, text="Select seqTabClean.RData ...",
                                     command=self.selChimera, font="Helvetica 12")
        self.outdirBtn=tk.Button(self.Frame, text="Select output directory ...",
                                     command=self.selOutdir, font="Helvetica 12")
        self.runBtn=tk.Button(self.Frame, text="RUN",
                              command=self.runTrackerScript, font="Helvetica 12")

        self.labelFilter = tk.Label(self.Frame, text="Select filtering.RData before continuing", font="Helvetica 10")
        self.labelDada = tk.Label(self.Frame, text="Select dada.RData before continuing", font="Helvetica 10")
        self.labelMerge = tk.Label(self.Frame, text="Select mergedReads.RData before continuing", font="Helvetica 10")
        self.labelSeqtab = tk.Label(self.Frame, text="Select seqTabRaw.RData before continuing", font="Helvetica 10")
        self.labelChimera = tk.Label(self.Frame, text="Select seqTabClean.RData before continuing", font="Helvetica 10")
        self.labelOutdir = tk.Label(self.Frame, text="Select output directory before continuing", font="Helvetica 10")

        self.filteringBtn.grid(row=1, column=1, pady=10, padx=5)
        self.dadaBtn.grid(row=2, column=1, pady=10, padx=5)
        self.mergeBtn.grid(row=3, column=1, pady=10, padx=5)
        self.seqtabBtn.grid(row=4, column=1, pady=10, padx=5)
        self.chimeraBtn.grid(row=5, column=1, pady=10, padx=5)
        self.outdirBtn.grid(row=6, column=1, pady=10, padx=5)
        self.runBtn.grid(row=7, column=2, pady=10, padx=5)

        self.labelFilter.grid(row=1, column=2, pady=10, padx=5)
        self.labelDada.grid(row=2, column=2, pady=10, padx=5)
        self.labelMerge.grid(row=3, column=2, pady=10, padx=5)
        self.labelSeqtab.grid(row=4, column=2, pady=10, padx=5)
        self.labelChimera.grid(row=5, column=2, pady=10, padx=5)
        self.labelOutdir.grid(row=6, column=2, pady=10, padx=5)

    def selFiltering(self):
        self.filteringDir=fd.askopenfilename()
        self.labelFilter.configure(text=self.filteringDir)

    def selDada(self):
        self.dadaDir=fd.askopenfilename()
        self.labelDada.configure(text=self.dadaDir)

    def selMerge(self):
        self.mergeDir=fd.askopenfilename()
        self.labelMerge.configure(text=self.mergeDir)

    def selSeqtab(self):
        self.seqtabDir=fd.askopenfilename()
        self.labelSeqtab.configure(text=self.seqtabDir)

    def selChimera(self):
        self.chimeraDir=fd.askopenfilename()
        self.labelChimera.configure(text=self.chimeraDir)

    def selOutdir(self):
        self.outDir=fd.askdirectory()
        self.labelOutdir.configure(text=self.outDir)

    def runTrackerScript(self):
        #check if all necessary inputs were made
        if (self.filteringDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Filtering.RData missing!")
            return

        if (self.dadaDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="dada.RData missing!")
            return

        if (self.mergeDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="mergedReads.RData missing!")
            return

        if (self.seqtabDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="seqTabRaw.RData missing!")
            return

        if (self.chimeraDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="seqTabClean.RData missing!")
            return

        if (self.outDir == ""):
            tk.messagebox.showinfo(title="Data missing",
                                   message="Output directory missing!")
            return

        #producing command line to run R Script

        process="Rscript"
        pathToScript="/project/genomics/Christoph/DADA2/tracker.R"
        fArg=self.filteringDir
        dArg=self.dadaDir
        mArg=self.mergeDir
        sArg = self.seqtabDir
        cArg = self.chimeraDir
        oArg = self.outDir

        commandLine=[process, pathToScript,
                     "-f", fArg,
                     "-d", dArg,
                     "-m", mArg,
                     "-s", sArg,
                     "-c", cArg,
                     "-o", oArg]

        #run the R script
        try:
            sp.check_call(commandLine)
        except sp.CalledProcessError:
            tk.messagebox.showerror(title="Error in calling R script",
                                    message="Execution of script tracker.R failed")
        else:
            tk.messagebox.showinfo(title="tracker.R",
                                   message="Script tracker.R finished successfully")

    def onClose(self):
        """destructor"""
        pub.sendMessage('subWindowClosed')
        self.destroy()

class mainFrame(tk.Frame):

    def __init__(self, parent):
        """Constructor Main Frame"""
        super().__init__()

        self.root=parent

        self.initUI()

        pub.subscribe(self.listener, 'subWindowClosed')


    def initUI(self):
        #tk.messagebox.showinfo(title="Updates ongoing",
        #                       message="Currently, the interface and pipeline are being updated. " +
        #                               "In case you experience any errors or shortcomings, please, let me know.\n" +
        #                               "Christoph")
        self.root.title('DADA2 startup')
        self.frame=tk.Frame()
        self.frame.pack()
        self.versionSelection=tk.StringVar()

        #choices for version selection
        # check which DADA2 installations are available
        self.versionsStable={}
        self.versionsDev={}
        with open('/project/genomics/Christoph/DADA2/package/versionsDADA2.txt', 'r') as f:
            next(f)  # skip headings
            reader = csv.reader(f, delimiter='\t')
            for version, path, status in reader:
                if status == 'stable':
                    self.versionsStable[version + ' (stable)'] = (version, path, status)
                elif status == 'experimental':
                    self.versionsDev[version + ' (experimental)'] = (version, path, status)

        #create list with all possible choices
        self.choices = [str(x) for x in self.versionsStable.keys()] + [str(y) for y in self.versionsDev.keys()]

        selBtn=tk.Button(self.frame, text='File selection and quality plots',
                            command=self.selectFrame)
        filtBtn=tk.Button(self.frame, text='Quality filtering of reads',
                            command=self.filterFrame)
        denoiseBtn=tk.Button(self.frame, text='Start denoising of reads',
                               command=self.denoiseFrame)
        taxnonmyBtn=tk.Button(self.frame, text='Taxonomic annotation',
                                command=self.taxonomyFrame)
        treeBtn=tk.Button(self.frame, text='Phylogenetic tree calculation',
                                command=self.phyloFrame)
        trackerBtn=tk.Button(self.frame, text='Tracking of sequence numbers',
                                command=self.sequenceTracker)
        versionLabel=tk.Label(self.frame, text="DADA2 version used:", font="Helvetica 10", )
        versionDD=tk.OptionMenu(self.frame, self.versionSelection, *self.choices)

        selBtn.pack(fill=tk.X, pady=10, expand=True)
        filtBtn.pack(fill=tk.X, pady=10, expand=True)
        denoiseBtn.pack(fill=tk.X, pady=10, expand=True)
        taxnonmyBtn.pack(fill=tk.X, pady=10, expand=True)
        treeBtn.pack(fill=tk.X, pady=10, expand=True)
        trackerBtn.pack(fill=tk.X, pady=10, expand=True)
        versionLabel.pack(fill=tk.X, pady=10, expand=True)
        versionDD.pack(fill=tk.X, pady=10, expand=True)

        #set default version to latest stable
        self.versionSelection.set(max(self.versionsStable.keys()))

    def show(self):
        """shows main frame"""
        self.root.deiconify()

    def hide(self):
        """hides main frame"""
        self.root.withdraw()

    def listener(self):
        """Pubsub listener function opens main frame after sub frames close"""
        self.show()

    def selectFrame(self):
        """opens selectFile window"""
        self.hide()
        if self.versionSelection.get() in self.versionsStable:
            subFrame = selectFile(version=self.versionsStable[self.versionSelection.get()])
        elif self.versionSelection.get() in self.versionsDev:
            subFrame = selectFile(version=self.versionsDev[self.versionSelection.get()])
        else:
            pub.sendMessage('subWindowClosed')
            tk.messagebox.showerror(title="DADA2 version unknown",
                                    message="Selected DADA2 version not available.")

    def filterFrame(self):
        """opens filterReads window"""
        self.hide()
        if self.versionSelection.get() in self.versionsStable:
            subFrame = filterReads(version=self.versionsStable[self.versionSelection.get()])
        elif self.versionSelection.get() in self.versionsDev:
            subFrame = filterReads(version=self.versionsDev[self.versionSelection.get()])
        else:
            pub.sendMessage('subWindowClosed')
            tk.messagebox.showerror(title="DADA2 version unknown",
                                    message="Selected DADA2 version not available.")

    def denoiseFrame(self):
        """opens denoisingReads window"""
        self.hide()
        if self.versionSelection.get() in self.versionsStable:
            subFrame = denoiseReads(version=self.versionsStable[self.versionSelection.get()])
        elif self.versionSelection.get() in self.versionsDev:
            subFrame = denoiseReads(version=self.versionsDev[self.versionSelection.get()])
        else:
            pub.sendMessage('subWindowClosed')
            tk.messagebox.showerror(title="DADA2 version unknown",
                                    message="Selected DADA2 version not available.")

    def taxonomyFrame(self):
        """opens taxonomyReads window"""
        self.hide()
        if self.versionSelection.get() in self.versionsStable:
            subFrame = taxonomyReads(version=self.versionsStable[self.versionSelection.get()])
        elif self.versionSelection.get() in self.versionsDev:
            subFrame = taxonomyReads(version=self.versionsDev[self.versionSelection.get()])
        else:
            pub.sendMessage('subWindowClosed')
            tk.messagebox.showerror(title="DADA2 version unknown",
                                    message="Selected DADA2 version not available.")

    def phyloFrame(self):
        """opens taxonomyReads window"""
        self.hide()
        subFrame = phyloTree()

    def sequenceTracker(self):
        """opens taxonomyReads window"""
        self.hide()
        if self.versionSelection.get() in self.versionsStable:
            subFrame = sequenceTracker(version=self.versionsStable[self.versionSelection.get()])
        elif self.versionSelection.get() in self.versionsDev:
            subFrame = sequenceTracker(version=self.versionsDev[self.versionSelection.get()])
        else:
            pub.sendMessage('subWindowClosed')
            tk.messagebox.showerror(title="DADA2 version unknown",
                                    message="Selected DADA2 version not available.")

if __name__ == '__main__':
    root=tk.Tk()
    root.geometry('200x400')
    app=mainFrame(root)
    root.mainloop()
