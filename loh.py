import vcf
import vcf.filters
import sys
import numpy as np
from glob import glob
import os
import re


frac = 0.89

os.chdir("/home/mitchell/IW/hmm/freebayes")
debug = False
def printd(s):
    if(debug):
        print(s)


def key_func(s):
    return [int(x) if x.isdigit() else x for x in re.findall(r'\D+|\d+', s)]

def writeFile(fName, dic):
    f = open(fName, "w+")
    for key in sorted(dic, key=key_func):
        f.write(key + "\t" + str( dic[key] ) + "\n")
    f.close()

f = open("var.vcf", "rt")
s288c_vcf = vcf.Reader(f)
pos = {}
num = 0
for r in s288c_vcf:
    con1 = r.INFO['TYPE'][0] == 'snp'
    dp = float(r.INFO['DP'])
    ao = float( r.INFO['AO'][0] )

    if( con1 and dp > 50 and ao/dp >= 0.9):
        key = r.CHROM +"\t" 
        value = r.POS
        pos[key+str(value)] = "1\t0\t30\t0"
        num += 1 
printd(num)

files = glob("cis*vcf")
files.extend(glob('uv*vcf'))

for f in files:
    treated = vcf.Reader(open(f, "rt")) 
    pos2 = dict(pos)
    name = str(f).split(".")[0] + ".hmm" 
    
    for r in treated:
        key = r.CHROM +"\t"+ str(r.POS)
        dp = float( r.INFO['DP'] )
        ao = float (r.INFO['AO'][0] ) 
        con1 =  (ao/dp > frac or ao/dp < 1-frac)
        con2 = dp > 15
        state = str( int(con1 and con2) ) 
        if key in pos:
            pos2[key] = state + "\t" + str(ao) + "\t" + str(dp) +"\t"+ str(ao/dp)
    
    writeFile(name, pos2)

    printd(pos2.values())
    printd(len(pos2))





# hmm stuff
'''
p = []
x = []
for i in pos.values():
    p.append(i - 1)
    x.append(i)

obs = np.transpose(np.array( [x, p] ))
x = pos.values()
printd(type(obs))
print(obs)
m = hmm.GaussianHMM(n_components=2) 
m.fit( [obs] ) 
Z2 = m.predict(pos.values())
'''



