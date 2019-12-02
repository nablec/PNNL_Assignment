#!/usr/bin/env python
# coding: utf-8



from collections import defaultdict
import numpy as np
import pandas as pd
import mygene,re
import matplotlib.pyplot as plt


#part2
def add_gene(phos,query,refseq):
    phos.insert(0,'Gene',phos["RefSeq"])
    for i,x in enumerate(refseq):
        phos.loc[phos["Gene"]==x,"Gene"] = query["symbol"][i]
    return phos

# part3
# read H_sapiens_RefSeq.fasta to a dictionary, ignore refseq which are not in the phosphopeptides.txt
def read_sequence(s,refseq):
    seq = ""
    ref = s[0].strip()[1:]
    dic = defaultdict()
    for line in s:
        line = line.strip()
        if line[0]==">":
            temp = line[1:]
            if seq!="":
                dic[ref]=seq
                seq = ""
            ref = temp
        elif ref not in refseq:
            continue
        else:
            seq +=line
    return dic
# df = pd.DataFrame(data=data_frame)


# part3 mapping
def add_site(phos,dic):
    phos["Site"] = phos["Gene"]
    for ref,seq in dic.items():
        pair = phos.loc[phos["RefSeq"]==ref,"Peptide":"Site"]
        site = "".join(map(str,pair["Site"].unique()))
        for peptide in pair["Peptide"]:
            protein = "".join(peptide.split(".")[1:-1])
            subseq = protein.replace("*","")
            pep_start = seq.find(subseq)
            protein_seq = "-"
            count=0
            for num in re.finditer(r"\*",protein):
                num = num.start()
                protein_seq += protein[num-1]+str(pep_start+num+count)
            pair.loc[pair["Peptide"]==peptide,"Site"]=site+protein_seq
    #     print(pair["Site"])
        phos.loc[phos["RefSeq"]==ref,"Site"]=pair["Site"]
    return phos

# part4
def find_series(sites):
    site_dic = defaultdict()
    cha,num = [],[]
    for site in sites:
        site=site.split("-")[-1]
        cha += re.findall("[A-Z]+",site)
        num += re.findall("[0-9]+",site)
    ser = list(map(lambda x:x[0]+x[1],zip(cha,num)))
    for i,x in enumerate(num):
        site_dic[int(x)]=ser[i]
    return site_dic

# plotting each sequence
def single_seq_plot(site_dic,seq,gene):
    
    y = np.array(list(site_dic.keys()))
    plt.figure(figsize=(20,0.7))
    plt.ylim(0,1) 
    plt.xlim(0,len(seq))
    plt.title(gene,size=20)
    plt.xlabel("amino acid position",size=10)
    figure = plt.bar(y,height=1,width=3*(len(seq)/2000))
    for i,bar in enumerate(figure):
        plt.annotate(list(site_dic.values())[i],xy=(bar.get_x(), 1),size=10, 
                     xycoords='data',xytext=(-20+(i%19)*(-20)*(-1*(i%2)), 30+(i%31)*15), textcoords='offset points',
                     bbox=dict(boxstyle="round", fc="0.7",alpha=0.5),
                     arrowprops=dict(arrowstyle="->",alpha=0.2))
#     plt.show()
#     plt.tight_layout()
    plt.savefig('picture/%s.png'%gene,dpi=100,bbox_inches = "tight")
    plt.show()
    plt.close()
#     print('%s.png saved'%gene)
    
def global_plot(phos,dic):        
    for ref,seq in dic.items():
        ref_list = phos.loc[phos["RefSeq"]==ref,"Gene":"Site"]
        gene = ref_list["Gene"].unique()[0]
        if type(gene)==float:
            continue
        sites = list(ref_list["Site"])
        site_dic = find_series(sites)
#  plot sequence
        single_seq_plot(site_dic,seq,gene)



if __name__ == "__main__":
# init mygene    
    mg = mygene.MyGeneInfo()
# read files
    phos = pd.read_csv('phosphopeptides.txt',sep = "\t",skiprows=0)
    H_sapiens = open("H_sapiens_RefSeq.fasta").readlines()
# unique RefSeq in phos dataframe
    refseq = np.unique(phos["RefSeq"])
    
# Question2     
# connect query to get gene symbol through refseq
    query = mg.querymany(refseq, scopes='refseq',fields = ["symbol"],as_dataframe=True)   
    phos = add_gene(phos,query,refseq)
    
# Question3        
# transfer sequence file to dictionary
    dic = read_sequence(H_sapiens,refseq)
    phos = add_site(phos,dic)
    phos.to_csv("New_phos.csv",header=True,sep=" ")
# Question4
#    global_plot(phos,dic)


