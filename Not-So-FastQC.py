#%%
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from Bio import SeqIO
import matplotlib.pyplot as plot
import pandas as pd
import numpy as np
from scipy.stats import norm
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

#No known issues
def perBaseSeqQuality(qualityScoreDataFrame):

    notebook.add(tab1, text='Per Base Sequence Quality')

    global frame1
    frame1 = tk.Frame(tab1)
    frame1.pack(side="top", expand=False)
   
    #Creates a data frame from the list and plots the data:
    f=plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.axhspan(28, 40, color='lightgreen')
    plot.axhspan(20, 28, color='lightyellow')
    plot.axhspan(0, 20, color='lightcoral')
    plot.xlabel('Position in read (bp)')
    plot.title('Quality score across all bases')
    plot.ylim(0,40)
    qualityScoreDataFrame.boxplot(showfliers=False, manage_ticks=True, 
    patch_artist=True, autorange=False, whis=(10,90))
    plot.xticks()
    plot.locator_params(axis='x', nbins=12)
    plot.grid()
    plot.close()

    canvas = FigureCanvasTkAgg(f, master = frame1)
    canvas.draw()
    canvas.get_tk_widget().pack() 

def sequenceLengthDistribution(sequenceList):
    notebook.add(tab2, text='Sequence Length Distribution')

    global frame2
    frame2 = tk.Frame(tab2)
    frame2.pack(side="top", expand=False)
    
    sequenceLengthList = [len(sequenceList[i]) for i in range(0, len(sequenceList))]
    align = 'mid'
    if sequenceLengthList.count(sequenceLengthList[0]) == len(sequenceLengthList):
        align = 'left'
    f=plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.xlabel('Sequence Length (bp)')
    plot.title('Distribution of sequence lengths over all sequences')
    n, bins, patches = plot.hist(sequenceLengthList, color='red', rwidth= 0.5, align=align)
    plot.xticks(np.round(bins))
    plot.grid()
    plot.close()
    
    canvas = FigureCanvasTkAgg(f, master = frame2)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def perSeqQualityScore(qualityScoreDataFrame):
    notebook.add(tab3, text='Per Sequence Quality Score')

    global frame3
    frame3 = tk.Frame(tab3)
    frame3.pack(side="top", expand=False)

    seqMeanQuality=qualityScoreDataFrame.mean(1)
    f=plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.xlabel('Mean Sequence Quality (Phred Score)')
    plot.title('Quality score distribution over all sequences')
    left, bins, right = plot.hist(seqMeanQuality, color='red', bins=20, align='left' )
    plot.xticks(np.round(bins))
    plot.locator_params(axis='x', nbins=40)
    #plot.xlim(0,40)
    plot.grid()
    plot.close()

    canvas = FigureCanvasTkAgg(f, master = frame3)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def perSeqGCcontent(seqGCcontent):
    notebook.add(tab4, text='Per Sequence GC-Content')

    global frame4
    frame4 = tk.Frame(tab4)
    frame4.pack(side="top", expand=False)

    f = range(0,100), plot.figure()
    f=plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.xlabel('GC-content')
    plot.title('GC-content distribution over all sequences')
    plot.hist(seqGCcontent, bins=25, density=True, color='red', align='left')
    
    #calculates mean and std for a probability density function (PDF):
    mean, std = norm.fit(seqGCcontent)

    #creates a set of evenly spaced numbers and plots it as a PDF into the histogram:
    xmin, xmax = plot.xlim()
    line = np.linspace(xmin, xmax)
    pdf = norm.pdf(line, mean, std)
    plot.plot(line, pdf, linewidth=5, color='blue', label='Theoretical Distribution')
    plot.legend(loc='upper right', fontsize=15)
    plot.grid()
    plot.close('all')

    canvas = FigureCanvasTkAgg(f, master = frame4)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def perBaseSeqContent(sequenceList, maxSeqLength, meanSeqLength):
    notebook.add(tab5, text='Per Base Sequence Content')

    global frame5
    frame5 = tk.Frame(tab5)
    frame5.pack(side="top", expand=False)

    iterationList, Tcontent, Acontent, Ccontent, Gcontent, = [], [], [], [], []
    #Counts the occurences of each base in each position:
    for i in range(0, maxSeqLength):
        for seq in sequenceList:
            try: iterationList.append(seq[i])
            except: continue 
        Tcontent.append(iterationList.count('T')/len(iterationList)*100)
        Acontent.append(iterationList.count('A')/len(iterationList)*100)
        Ccontent.append(iterationList.count('C')/len(iterationList)*100)
        Gcontent.append(iterationList.count('G')/len(iterationList)*100)
        iterationList=[]
        
    #Reduces the data so that the plot doesn't look too clamped 
    if meanSeqLength > 300:
        samplingFactor=int(maxSeqLength/100)
        x=list(range(0, maxSeqLength, samplingFactor))
        Tcontent = [Tcontent[i] for i in x]
        Acontent = [Acontent[i] for i in x]
        Ccontent = [Ccontent[i] for i in x]
        Gcontent = [Gcontent[i] for i in x]              

    #Plots the data:
    Xaxis, f = range(0,len(Tcontent)), plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.ylim(0,100)
    if maxSeqLength > 300: plot.xlabel('Position in read (bp, sampled)') 
    else: plot.xlabel('Position in read (bp)') 
    plot.title('Sequence content across all bases')
    plot.plot(Xaxis,Tcontent, color='red' ,label='%T',linewidth=3)
    plot.plot(Xaxis,Acontent, color='green',label='%A',linewidth=3)
    plot.plot(Xaxis,Gcontent, color= 'black',label='%G',linewidth=3)
    plot.plot(Xaxis,Ccontent,color='blue',label='%C',linewidth=3)
    plot.legend(loc='upper right', fontsize=15)
    plot.grid()
    plot.close()

    canvas = FigureCanvasTkAgg(f, master = frame5)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def perBaseNcontent(sequenceList, maxSeqLength):
    notebook.add(tab6, text="Per Base 'N' Content")

    global frame6
    frame6 = tk.Frame(tab6)
    frame6.pack(side="top", expand=False)

    iterationList, NcontentList, c =[], [], 0
    #Extracts the porcentage of Ns in the sequences:
    for i in range(0,maxSeqLength):
        for i in sequenceList:
            try: iterationList.append(i[c])
            except: continue
        NcontentList.append(iterationList.count('N')/len(iterationList)*100)
        iterationList=[]
        c+=1
    Xaxis=range(1, maxSeqLength+1)
    f=plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.ylim(0,100)
    plot.xlabel('Position in read (bp)')
    plot.title('N content across all bases')
    plot.plot(Xaxis, NcontentList, color='red' ,label='%N',linewidth=6)
    plot.grid()
    plot.close()

    canvas = FigureCanvasTkAgg(f, master = frame6)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def adapterContent(sequenceList, maxSeqLength):
    notebook.add(tab7, text="Adapter Content")

    global frame7
    frame7 = tk.Frame(tab7)
    frame7.pack(side="top", expand=False)

    iua, isa3, isa5, nts, ssa = [], [], [], [], []
    iuaPlot, isa3Plot, isa5Plot, ntsPlot, ssaPlot = [], [], [], [], []

    #find position of adapters:
    for seq in sequenceList:  
        if seq.find('CTGTCTCTTATA') != -1:
            iua.append(seq.find('AGATCGGAAGAG')) #Illumina Universal Adapter
            isa3.append(seq.find('TGGAATTCTCGG')) #Illumina Small RNA 3' Adapter
            isa5.append(seq.find('GATCGTCGGACT')) #Illumina Small RNA 5' Adapter
            nts.append(seq.find('CTGTCTCTTATA')) #Nextera Transposase Sequence
            ssa.append(seq.find('CGCCTTGGCCGT')) #SOLID Small RNA Adapter

    #calculates cumulative porcentages to plot:
    cumulative1, cumulative2, cumulative3, cumulative4, cumulative5, = 0, 0, 0, 0, 0, 
    for i in range(0, maxSeqLength-11):
        try:
            cumulative1 += iua.count(i)
            iuaPlot.append((cumulative1/len(sequenceList))*100) #Illumina Universal Adapter

            cumulative2 += isa3.count(i)
            isa3Plot.append((cumulative2/len(sequenceList))*100) #Illumina Small RNA 3' Adapter

            cumulative3 += isa5.count(i)
            isa5Plot.append((cumulative3/len(sequenceList))*100) #Illumina Small RNA 5' Adapter

            cumulative4 += nts.count(i)
            ntsPlot.append((cumulative4/len(sequenceList))*100) #Nextera Transposase Sequence

            cumulative5 += ssa.count(i)
            ssaPlot.append((cumulative5/len(sequenceList))*100) #SOLID Small RNA Adapter
        except: continue
    
    #Plots cumulative porcentages as lines:
    Xaxis, f = range(0, maxSeqLength-11), plot.figure()
    f.set_figwidth(20)
    f.set_figheight(12)
    plot.ylim(-0.3,100)   
    plot.title('Adapter content')
    plot.plot(Xaxis,iuaPlot, color='red' ,label='Illumina Universal Adapter',linewidth=2.5)
    plot.plot(Xaxis,isa3Plot, color='blue' ,label="Illumina Small RNA 3' Adapter",linewidth=2.5)
    plot.plot(Xaxis,isa5Plot, color='green' ,label="Illumina Small RNA 5' Adapter",linewidth=2.5)
    plot.plot(Xaxis,ntsPlot, color='black' ,label='Nextera Transposase Sequence',linewidth=2.5)
    plot.plot(Xaxis,ssaPlot, color='purple' ,label='SOLID Small RNA Adapter',linewidth=2.5)
    plot.legend(loc='upper right', fontsize=15)
    plot.grid()
    plot.close()

    canvas = FigureCanvasTkAgg(f, master = frame7)  
    canvas.draw()
    canvas.get_tk_widget().pack()

def start():
    
    #Sets a loading message:
    loadingMessage = tk.Label(window, text='LOADING... ', font="Times 30",)
    loadingMessage.place(x=300,y=5)

    #Opens file dialog and destroys loading message once process is done:
    message1.destroy()
    fileName = filedialog.askopenfilename(initialdir = "Users\luan_\Documents\School work\Programming project", title = "Select a File",  filetypes = (("FastQ files", "*.fastq  .fq"),  ("all files", "*.*")))
    loadingMessage.destroy()

    #Parse fastQ file and extracts needed information:
    maxSeqLength=0
    minSeqLength=1000000000
    qualityScoreList=[]
    seqList=[]
    seqCount=0
    seqLength=0
    data=SeqIO.parse(fileName, 'fastq')
    for i in data:
        seqCount+=1
        if len(i.seq) > maxSeqLength: maxSeqLength=len(i.seq)
        if len(i.seq) < minSeqLength: minSeqLength=len(i.seq)
        seqLength+=len(i.seq)
        qualityScoreList.append(i.letter_annotations['phred_quality'])
        seqList.append(i.seq)
    meanSeqLength = round(seqLength/seqCount)

    #Reduces the data to improves computing time if the file is too big:
    if checkbox_var.get() == '1':
            sequenceList = seqList
            qualityScoreDataFrame = pd.DataFrame(qualityScoreList)    
    else:
        if seqCount > 100000:
            samplingFactor=seqCount/10000
            x=list(range(0, seqCount, int(samplingFactor)))
            sequenceList = [seqList[i] for i in x]
            scoreData=pd.DataFrame([qualityScoreList[i] for i in x])
        else: 
            sequenceList = seqList
            scoreData = pd.DataFrame(qualityScoreList)

        #Reduces the quality data if the sequence is too long, but preserves the first 10 nt:
        if meanSeqLength > 300:   
            samplingFactor=int(maxSeqLength/100)    
            headbp=list(range(0,11))
            reducingIndex=list(range(11, maxSeqLength, samplingFactor))
            reducedData=headbp+reducingIndex
            try: qualityScoreDataFrame = scoreData[reducedData]
            except: qualityScoreDataFrame = scoreData
        else:
            qualityScoreDataFrame = scoreData
   
    #Extracts the GC content to be show on general information: 
    seqGCcontent=[]
    for seq in sequenceList:
        gc=(seq.count('G') + seq.count('C'))/len(seq)*100
        seqGCcontent.append(round(gc))
    GCcontent = max(set(seqGCcontent), key= seqGCcontent.count)

    #Extracts the sequence length to be show on general information: 
    if minSeqLength == maxSeqLength:
        sequenceLength = minSeqLength
    else:
        sequenceLength = str(minSeqLength)+' - '+str(maxSeqLength)

    #Generates general information to be shown on the first tab:
    global message2
    try: message2.destroy() 
    except:pass

    notebook.add(tab0, text='General Information')
    message2 = tk.Label(tab0, text='Open file: ' +fileName+'\n\n\n\n\n\n'
                    'Total sequences: '+str(seqCount)+'\n\n'+
                    'GC content: '+ str(GCcontent)+'%'+'\n\n'+         
                    'Sequence length: '+str(sequenceLength),             font='Arial 20 bold',)
    message2.pack(fill='both', expand=False)
    message2.place(x=20,y=100)

    #Overwrites current plots:
    try: 
        frame1.destroy()
        frame2.destroy()
        frame3.destroy()
        frame4.destroy()
        frame5.destroy()
        frame6.destroy()
        frame7.destroy()
    except:pass

    
    perBaseSeqQuality(qualityScoreDataFrame)
    sequenceLengthDistribution(sequenceList)
    perSeqQualityScore(qualityScoreDataFrame)
    perSeqGCcontent(seqGCcontent)
    perBaseSeqContent(sequenceList, maxSeqLength, meanSeqLength)
    perBaseNcontent(sequenceList, maxSeqLength)
    adapterContent(sequenceList, maxSeqLength)

#Creats UI window:
window = tk.Tk()
window.title('Not-So-FastQC')
window.geometry('2000x1300')

#Creates tabs for the plots:
notebook = ttk.Notebook(window)
notebook.pack(pady=10, expand=True)
notebook.place(x=0,y=100)
tab0 = ttk.Frame(notebook, width=2000, height=1000)
tab1 = ttk.Frame(notebook, width=2000, height=1000)
tab2 = ttk.Frame(notebook, width=2000, height=1000)
tab3 = ttk.Frame(notebook, width=2000, height=1000)
tab4 = ttk.Frame(notebook, width=2000, height=1000)
tab5 = ttk.Frame(notebook, width=2000, height=1000)
tab6 = ttk.Frame(notebook, width=2000, height=1000)
tab7 = ttk.Frame(notebook, width=2000, height=1000)
tab0.pack(fill='both', expand=1)
tab1.pack(fill='both', expand=1)
tab2.pack(fill='both', expand=1)
tab3.pack(fill='both', expand=1)
tab4.pack(fill='both', expand=1)
tab5.pack(fill='both', expand=1)
tab6.pack(fill='both', expand=1)
tab7.pack(fill='both', expand=1)

#Writes messages on the UI window;
message = tk.Label(window, text="Welcome to Not-So-FastQC!", font="Times 50 bold italic",)
message.pack()
message1 = tk.Label(window, text='\n\n\n\n\n\n\n\nNot-So-FastQC High Throughput Sequence QC Report \n Version 1.0 \n Developed by Luan Gardenalli, 2022-12.', font="Times 30",)
message1.pack(fill='both')

#Adds button and checkbox to UI window:
startButton = tk.Button(window, text='Open file', command=start, width=10, height=2)
startButton.place(x=5,y=5)
checkbox_var = tk.StringVar()
checkBox = ttk.Checkbutton(window, text= 'Bypass Data Reduction', variable=checkbox_var)
checkBox.place(x=100,y=5)


window.mainloop() 
#%%