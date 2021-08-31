#!/usr/bin/env python
# coding: utf-8

# In[76]:


#Neat notebook for finding core reactions in Human1 model.
import cobra
import matplotlib.pyplot as plt
import gurobipy as gp
import pandas as pd
from pandas import DataFrame
import numpy as np


# In[230]:


#Upload Human1 model.
Human1 = "Human-GEM.xml"
model = cobra.io.read_sbml_model(Human1)
model


# In[243]:


x = model.reactions.get_by_id
y = model.metabolites.get_by_id


# In[232]:


#Check we can access the model's properties properly. 
print("Reactions:", len(model.reactions))
print("Metabolites:", len(model.metabolites))
print("Genes:", len(model.genes))


# In[233]:


#Finding biomass reactions, so I can pick objective function. 
for r in model.reactions:
    if "biomass" in r.name:
        print(r.id, ': ', r.name)


# In[234]:


#Running FBA, with biomass_human as objective function. 
biomass_rxn = model.reactions.get_by_id("biomass_human")
biomass_rxn.objective_coefficient = 1
biomass_rxn = model.reactions.get_by_id("HMR_10062")
biomass_rxn.objective_coefficient = 0
biomass_rxn = model.reactions.get_by_id("HMR_10063")
biomass_rxn.objective_coefficient = 0
biomass_rxn = model.reactions.get_by_id("HMR_10064")
biomass_rxn.objective_coefficient = 0
biomass_rxn = model.reactions.get_by_id("HMR_10065")
biomass_rxn.objective_coefficient = 0
solution = model.optimize()
print("FBA status:", solution.status)
print("FBA solution:", solution.objective_value)


# In[235]:


#Finding number of core reactions (reactions where flux is not zero, i.e. below or above zero). != is not equal to. 
core_reactions = []
all_other_reactions = []

for r in model.reactions:
    if solution.fluxes[r.id]>0:
        core_reactions.append(r.id)
        continue
    if solution.fluxes[r.id]<0:
        core_reactions.append(r.id)
        continue
    else:
        all_other_reactions.append(r.id)
        continue 
            
print(('Amount of core reactions: '),(len(core_reactions)))
print(('Head of core reactions: '),(core_reactions[0:5]))
print('\n')
print(('Other reactions (zero flux): '),(len(all_other_reactions)))
print(('Head of all other reactions: '),(all_other_reactions[0:5]))


# In[7]:


#Converting core reactions list into DF so it can be added as an extra column to Human1 CSV, ready for model filtering. 
core_reactionDF = pd.DataFrame(core_reactions)
core_reactionDF.to_csv(r'/Users/katemeeson/Desktop/core_reaction_ids.csv', index = False, header=True)


# In[236]:


#Upload CSV with core reactions added as new column. Now filter Human1 model to just show core reactions. 
filtering_model = pd.read_csv('/Users/katemeeson/Dropbox (The University of Manchester)/filtering_human1_for_core_reactions.csv')
filtering_model


# In[237]:


#Assigning series names to CoreReactions and ReactionID column. 
Core_reactions = filtering_model.iloc[:,0]
ReactionIDs = filtering_model.iloc[:,1]


# In[238]:


filtered_model = []

for i in range(13096): #number of Human1 reactions to search. 
    for x in range(840): #series of core reactions to filter model against. 
        if ReactionIDs[i]==Core_reactions[x]:
            print(filtering_model.iloc[i,1:4])
            print('\n')
            filtered_model.append(filtering_model.iloc[i,1:4])
        else:
            continue


# In[239]:


#Converting filtered model list into DF. Now have a DF of Human1 model, filtered according to core/non-zero flux reactions.
filtered_model_df = pd.DataFrame(filtered_model,columns=['ReactionID','UniprotID','ReactionRule'])
filtered_model_df


# In[240]:


filtered_model_df.iloc[[0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],:]


# In[244]:


#Checking what some of the reactions are. 
x('HMR_8472')


# In[245]:


y('m01424n')


# In[17]:


#Convert DF to CSV. 
filtered_model_df.to_csv(r'/Users/katemeeson/Desktop/filtered_model_df.csv', index = False, header=True)


# In[246]:


len(model.reactions)


# In[23]:


for rule in filtered_model_df.iloc[:,2]:
    print(rule)


# In[247]:


#Organise the reaction rules. 
ANDs = []
ORs = []
ANDORs = []
one_gene = []
no_gene = []

for r in model.reactions:
    if 'and' in r.gene_reaction_rule and 'or' not in r.gene_reaction_rule:
        ANDs.append(r)
    if 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule:
        ANDORs.append(r)
    if 'or' in r.gene_reaction_rule and 'and' not in r.gene_reaction_rule:
        ORs.append(r)
    if len(r.gene_reaction_rule) == 0:
        no_gene.append(r.id)
    if len(r.gene_reaction_rule) != 0:
        if 'or' in r.gene_reaction_rule:
            continue
        elif 'and' in r.gene_reaction_rule:
            continue
        else:
            one_gene.append(r.id)

print('AND rules: ', len(and_rules))
print('ANDOR rules: ', len(andor_rules))
print('OR rules: ', len(or_rules))
print('ONE GENE rules: ', len(one_gene))
print('NO GENE rules: ', len(no_gene))
print('Proportion of model not annotated: ', len(no_gene)/len(model.reactions)*100, '%')
print('Proportion of model which IS annotated: ', 100-(len(no_gene)/len(model.reactions)*100), '%')
print('Total Reactions = ', len(model.reactions))
print(len(or_rules)+len(and_rules)+len(andor_rules)+len(one_gene)+len(no_gene))
#Weakness is there's a whole load of reactions in the model which aren't covered by rules. 


# In[49]:


#Find if there are any in-built constraints. 

in_built_constraints = {} #Wavy brackets build dictionary. 

for r in model.reactions:
    if r.bounds != (0.0, 0.0) and r.bounds != (-1000.0, 0.0) and r.bounds != (-1000.0, 1000.0) and r.bounds != (0.0, 1000.0):
        in_built_constraints[r.id] = r.bounds
        print(r.id)

for key, val in in_built_constraints.items():
    print(key, val)
    
if len(in_built_constraints) == 0:
    print('No in-built constraints')


# In[248]:


#Running quick FBA to check model is feasible before you do an FVA and have to wait.
biomass_rxn = model.reactions.get_by_id("biomass_human")
biomass_rxn.objective_coefficient = 1
biomass_rxn = model.reactions.get_by_id("HMR_10062")
biomass_rxn.objective_coefficient = 0
biomass_rxn = model.reactions.get_by_id("HMR_10063")
biomass_rxn.objective_coefficient = 0
biomass_rxn = model.reactions.get_by_id("HMR_10064")
biomass_rxn.objective_coefficient = 0
biomass_rxn = model.reactions.get_by_id("HMR_10065")
biomass_rxn.objective_coefficient = 0
solution = model.optimize()
print("FBA status:", solution.status)
print("FBA solution:", solution.objective_value)

#FBA core reactions
core_reactions = []
all_other_reactions = []

for r in model.reactions:
    if solution.fluxes[r.id]>0:
        core_reactions.append(r.id)
        continue
    if solution.fluxes[r.id]<0:
        core_reactions.append(r.id)
        continue
    else:
        all_other_reactions.append(r.id)
        continue 
            
print(('Amount of FBA core reactions: '),(len(core_reactions)))


# In[249]:


#Trying an FVA just on biomass_human
print(cobra.flux_analysis.flux_variability_analysis(model, x('biomass_human')))


# In[250]:


#Running FVA. 
model.objective = x('biomass_human')
FVA = cobra.flux_analysis.flux_variability_analysis(model) #won't work


# In[251]:


#Code for finding core reactions using FVA (and FBA), when FVA is actually working. 
core_rxns = []

for r in model.reactions: 
    if FVA.at[r.id,'minimum']!=0 or FVA.at[r.id, 'maximum']!=0:
        #if FBA[r.id] !=0
        append.core_rxns[r.id]

print('Core reactions (FBA): ', len(core_reactions))
print(len('Core reactions (FVA): ', len(core_rxns)))

for r in core_reactions:
    if r not in core_rxns:
        print(r)


# In[252]:


#Check which reactions are annotated with genes. 
annotated = []
not_annotated = []

for r in core_reactions:
    if len(x(r).gene_reaction_rule) ==0:
        not_annotated.append(r)
    else:
        annotated.append(r)
        
print('Annotated reactions: ', len(annotated))
print('Reactions not annotated: ', len(not_annotated))
if len(core_reactions) == (len(annotated) + len(not_annotated)):
    print('Lengths match. All reactions searched')
print('Proportion of CORE reactions annotated: ', 100-(len(not_annotated)/len(core_reactions)*100), '%' )


# In[253]:


#To integrate, make dict to match genes to expression data. Then match expression to reactions as a library (nested dict).

#Uploading proteomics dataset and sorting into dataframe subsets. 
Myc_proteins = pd.read_csv('/Users/katemeeson/Dropbox (The University of Manchester)/UNIMAN Ovarian Cancer Project/July August modelling/ProteinOutput_FC_Group_Only.csv',index_col=0)
FC0 = Myc_proteins.iloc[:,0:3]
FC100 = Myc_proteins.iloc[:,[3,4,5]]
FC500 = Myc_proteins.iloc[:,[6,7,8]]


# In[254]:


#Converting df subsets into dictionaries
FC0_transposed = FC0.transpose() #orientation for dataframe to dictionary conversion. 
FC100_transposed = FC100.transpose()
FC500_transposed = FC500.transpose()

FC0_dict = FC0_transposed.to_dict('list')
FC100_dict = FC100_transposed.to_dict('list')
FC500_dict = FC500_transposed.to_dict('list')

print('FC0_dict: ', len(FC0_dict))
print('FC100_dict: ', len(FC100_dict))
print('FC500_dict: ', len(FC500_dict))


# In[255]:


#From Josef's coding. 


# Defining function to match expression data to reactions in the model
#this can be made better with a return function

def gene_expression_match(lis, data, lib): #lis = list (ANDs, ORs, ANDORs or ones), data = expression data,  
    for r in lis:                          #lib = library you are making
        ins= {}
        notins = {}
        for g in x(r).genes:
            for key, val in data.items():
                if key in str(g):
                    ins[str(g)] = val
                if key not in str(g):
                    notins[str(g)] = 'NO DATA'
        temp = {**notins, **ins}
        lib[r] = temp


# In[259]:


#The error seems to be when iterating through x(r).genes. Maybe it's because HMR_4137 has no genes annotated. 
#What format is Josef's model and proteomics data in? 
for r in ANDs:
    for g in x(r).genes:
        print(g)


# In[261]:


for r in data.items():
    print(r) #This line of the function does work so I don't think the error is here. 


# In[260]:


data = Myc_proteins #Trying to solve the function above as it is getting stuck at the first iteration below. 
data


# In[257]:


#From Josef's coding. 

# Matching genes and expression values to reactions
# Output is a nested dictionary, with key 1 = reaction ID, key 2 = gene ID and value = expression value or 'NO DATA' if
# no expression data available

FC0_ANDs = {}
FC100_ANDs = {}
FC500_ANDs = {}

gene_expression_match(ANDs, FC0_dict, FC0_ANDs)
gene_expression_match(ANDs, FC100_dict, FC100_ANDs)
gene_expression_match(ANDs, FC500_dict, FC500_ANDs)

FC0_ORs = {}
FC100_ORs = {}
FC500_ORs = {}

gene_expression_match(ORs, FC0_dict, FC0_ORs)
gene_expression_match(ORs, FC100_dict, FC100_ORs)
gene_expression_match(ORs, FC500_dict, FC500_ORs)

FC0_ANDORs = {}
FC100_ANDORs = {}
FC500_ANDORs = {}

gene_expression_match(ANDORs, FC0_dict, FC0_ANDORs)
gene_expression_match(ANDORs, FC100_dict, FC100_ANDORs)
gene_expression_match(ANDORs, FC500_dict, FC500_ANDORs)

FC0_ones = {}
FC100_ones = {}
FC500_ones = {}

gene_expression_match(one_gene, FC0_dict, FC0_ones)
gene_expression_match(one_gene, FC100_dict, FC100_ones)
gene_expression_match(one_gene, FC500_dict, FC500_ones)

#The problem in this might go back to it not being able to carry out FVA. 


# In[ ]:


#Work out the coverage of our dataset onto the core reactions. 


# In[ ]:


#Turn off glucose uptake (for example), and see how this affects the FBA solution.


# In[ ]:




