import MultiscaleMutationClusteringScript as M2C
from numpy import *	
import pylab as plt
import pandas as pd

data = pd.read_csv('E:/hotspot_test/m2c_data/thym_kmt2c.csv')


df = data[['S_ID','Cancer_type','gene_symbol','Amino_Acid_length','Protein_position','EFF']]
#將gene_symbol抓出來成list，去掉重複
gene = df['gene_symbol'].drop_duplicates().values.tolist()
#製作每個基因和其amino acid 長度
length = df[['gene_symbol','Amino_Acid_length']]
#將上一步做出來的data frame去掉重複做成tuple
tuples = tuple(length.drop_duplicates().itertuples(index=False))
'''
(Pandas(gene_symbol='MLLT10', Amino_Acid_length=1068),
 Pandas(gene_symbol='ERBB3', Amino_Acid_length=1342),
 Pandas(gene_symbol='CHD4', Amino_Acid_length=1937),
 Pandas(gene_symbol='FOXA1', Amino_Acid_length=472),
 Pandas(gene_symbol='ERBB2', Amino_Acid_length=1255),
 Pandas(gene_symbol='MTOR', Amino_Acid_length=2549),
 Pandas(gene_symbol='SF3B1', Amino_Acid_length=1304),
 Pandas(gene_symbol='PIK3CA', Amino_Acid_length=1068),
 Pandas(gene_symbol='MET', Amino_Acid_length=1408),
 Pandas(gene_symbol='EGFR', Amino_Acid_length=1210))

tuples[0]
Out[26]: Pandas(gene_symbol='MLLT10', Amino_Acid_length=1068)

tuples[0][0]
Out[27]: 'MLLT10'

tuples[0][1]
Out[28]: 1068
'''
#從tuple中一個一個基因抓出來找hotspot
for i in range(len(tuples)):
    
	#基因和其長度
    gene_symbol = tuples[i][0]
    interval = [0,tuples[i][1]] 
	
    #myall = df[(df['gene_symbol'] == gene_symbol)]
    #mylabel = df.set_index('gene_symbol')
    #把相對基因的all mutation抓出來
    myall = df[(df['gene_symbol'] == gene_symbol)]
    #myall = mylabel.loc[:tuples[i][0]]  
    #All_Mutations = []
    All_Mutations = myall['Protein_position'].values.tolist()
    All_Mutations = map(int, All_Mutations)
	
    #Synonymous_Mutations = []
    mysys = myall[(myall['EFF'] == 'synonymous_variant')]
    Synonymous_Mutations = mysys['Protein_position'].values.tolist()
    Synonymous_Mutations = map(int,Synonymous_Mutations)
	
    noise_estimate = 1.0*len(Synonymous_Mutations)/len(All_Mutations)
    bandwidths = [10,15,20]
    final_clusters, final_params, density_list = M2C.M2C(All_Mutations, noise_estimate, bandwidths = bandwidths, interval = interval, print_outs = False)

    numofmutation = int(len(All_Mutations)-len(Synonymous_Mutations))
    output_table = []
    for i in range(len(final_clusters)):
        muta = (str(gene_symbol)+" : "+str(min(final_clusters[i]))+"-"+str(max(final_clusters[i]))+":"+str(len(final_clusters[i])))
        output_table.append(muta)
	
    with open("E:/hotspot_test/m2c_data/%s.txt"%(gene_symbol), "w") as text_file:
        text_file.write("final_clusters_: {}".format(output_table))
    

	"""
	No Analysis Occurs Below - Results are just plotted using Matplotlib
	"""
    plt.figure(figsize=(20,15))
    color_list = ['blue', 'green', 'cyan', 'yellow', 'orange', 'red', 'purple']
    plt.subplot(211)
    plt.title("%s Mutation Histogram\nAll Mutations in gray with colored clusters overlayed on top\nNumber of nonsy_mutation =%i" %(gene_symbol,numofmutation))	#Create a histogram for all the data
    plt.hist(All_Mutations, bins = range(interval[1]), color = "gray")
	#Create a Histogram for each cluster on top of all the data.
    for i in range(len(final_clusters)):
        cluster = final_clusters[i]
        plt.hist(cluster, bins = range(interval[1]), color = color_list[i%len(color_list)], label = str(min(cluster))+"-"+str(max(cluster))+" : "+str(len(cluster)))
    plt.xlim(interval[0], interval[1])
    plt.ylabel("Mutation Count")
    plt.legend()

    plt.subplot(212)
    plt.title("KDE's and Final Mixture Model Gaussians")
    x_steps = arange(interval[0], interval[1], .1)
    #Plot each KDE
    for i in range(len(bandwidths)):
        density = density_list[i]
        h = bandwidths[i]
        plt.plot(x_steps, density, "--", label = "h="+str(h), color = str(.8*i/len(bandwidths)))
	#Plot Each Gaussian representing a cluster. Uniform noise is not plotted
    for i in range(len(final_params)):
        w, m, s = final_params[i]
        plt.plot(x_steps, [w/(sqrt(2*pi*s**2))*exp(-(x-m)**2/(2*s**2)) for x in x_steps], color_list[i%len(color_list)])
    plt.legend()
    plt.xlim(interval[0], interval[1])
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Mutation Density")
    #plt.show()
    plt.savefig('E:/hotspot_test/m2c_data/%s.png' % gene_symbol)
	