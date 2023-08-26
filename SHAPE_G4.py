import pandas as pds

df=pds.read_table('/cluster/home/zyw_jh/projects/structure/shapeG4/unsplice.info',sep="\t")
dd=dict(zip(df['rG4_intervals'],df['rG4_structural_class']))
dt=dict(zip(df['rG4_intervals'],df['Transcript_ids']))

fqfile=SeqIO.parse(open('/cluster/home/zyw_jh/projects/structure/shapeG4/hek293.g4.fa'),'fasta')
identifier=[]
sequence=[]
for seq1 in fqfile:
    identifier.append(seq1.id.split("::",1)[1].replace("(",":").replace(")",""))
    sequence.append(seq1.seq)
xdic=dict(zip(identifier,sequence))

out=open(r"/cluster/home/zyw_jh/projects/structure/shapeG4/G4.fa2type.tsv",'w+')
for key,value in dd.items():
    print(key+"\t"+value+"\t"+xdic[key]+"\t"+dt[key],file=out)
    
out.close()

tpm = pds.read_csv('/cluster/home/zyw_jh/projects/structure/DHX36/new_NAIdata/000_rtsc2react_sfalpha/abundance/wt_dmso_combined_TPM.cov1.csv')
tpm['tid']=[item[0] for item in tpm['transcript'].str.split(".")]
#tpm['transcript_filter']=tpm['transcript'].map(lambda x: str.split(x,'.')[0])

t2p=dict(zip(tpm['tid'],tpm['wt_dmso_combined_TPM']))

def get_max(x):
    largest_number = 0
    for number in x:
        if number is not None and largest_number < number:
            largest_number = number
    # print(largest_number)

data = []
a=0
out2=open("/cluster/home/zyw_jh/projects/structure/shapeG4/maxAbund.t2seq.tsv","w+")
for i in df['rG4_intervals']:
    tlist = dt[i].split(";")
    newlist = [item.split(".",1)[0] for item in tlist]
    sl = [t2p.get(t) for t in newlist]
    n2s = dict(zip(newlist,sl))
    newn2s = {k:v for k,v in n2s.items() if v is not None}
    #dict(zip([i+"\t"+maxt+"\t"+dd[i]],[str(xdic[i])]))   zip的两个参数都需要是list[]
    
    # remove none value from dict
    if bool(newn2s):
    #filter out empty dicts
        maxt = max(newn2s, key=newn2s.get)
        print(i+"\t"+maxt+"\t"+dd[i]+"\t"+str(xdic[i]),file=out2)
    else:
        a=a+1
        #count the number of empty dicts
     
out2.close()   
