#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import matplotlib.pyplot as plt
import gurobipy as gp
import pandas as pd
from pandas import DataFrame
import numpy as np


# In[12]:


model = cobra.io.read_sbml_model('Human-GEM.xml')
model


# In[3]:


x = model.reactions.get_by_id
y = model.metabolites.get_by_id


# In[4]:


print("Reactions:", len(model.reactions))
print("Metabolites:", len(model.metabolites))
print("Genes:", len(model.genes))


# In[6]:


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

print('AND rules: ', len(ANDs))
print('ANDOR rules: ', len(ANDORs))
print('OR rules: ', len(ORs))
print('ONE GENE rules: ', len(one_gene))
print('NO GENE rules: ', len(no_gene))
print('Proportion of model not annotated: ', len(no_gene)/len(model.reactions)*100, '%')
print('Proportion of model which IS annotated: ', 100-(len(no_gene)/len(model.reactions)*100), '%')
print('Total Reactions = ', len(model.reactions))
print(len(ORs)+len(ANDs)+len(ANDORs)+len(one_gene)+len(no_gene))


# In[7]:


in_built_constraints = {} 

for r in model.reactions:
    if r.bounds != (0.0, 0.0) and r.bounds != (-1000.0, 0.0) and r.bounds != (-1000.0, 1000.0) and r.bounds != (0.0, 1000.0):
        in_built_constraints[r.id] = r.bounds
        print(r.id)

for key, val in in_built_constraints.items():
    print(key, val)
    
if len(in_built_constraints) == 0:
    print('No in-built constraints')


# In[8]:


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


# In[13]:


model.objective = x('biomass_human')
FVA = cobra.flux_analysis.flux_variability_analysis(model)


# In[14]:


Myc_proteins = pd.read_csv('/Users/katemeeson/Dropbox (The University of Manchester)/UNIMAN Ovarian Cancer Project/July August modelling/ProteinOutput_FC_Group_Only.csv',index_col=0)
FC0 = Myc_proteins.iloc[:,0:3]
FC100 = Myc_proteins.iloc[:,[3,4,5]]
FC500 = Myc_proteins.iloc[:,[6,7,8]]

FC0_transposed = FC0.transpose() #orientation for dataframe to dictionary conversion. 
FC100_transposed = FC100.transpose()
FC500_transposed = FC500.transpose()

FC0_dict = FC0_transposed.to_dict('list')
FC100_dict = FC100_transposed.to_dict('list')
FC500_dict = FC500_transposed.to_dict('list')

print('FC0_dict: ', len(FC0_dict))
print('FC100_dict: ', len(FC100_dict))
print('FC500_dict: ', len(FC500_dict))


# In[17]:


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


# In[33]:


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


# In[30]:


#troubleshooting the function. 
if 'HMR_4137' in ANDs:
    print('there')
#I think the problem is that the rule lists, e.g. 'ANDs' arent built so you can access the reaction IDs properly, so the function can't search properly. 


# In[40]:


print(no_gene[0:10])
print(ANDs[0:10]) #Can search through one_gene and no_gene but not ANDs, ORs and ANDORs. 


# In[32]:


x('HMR_4137').genes
#Could it be to do with this frozenset() situation? 


# In[26]:


#Have tried the function in a new notebook, and tried it for just one sample. Still doesn't work. 
for r in ANDs:
    ins= {}
    notins = {}
    for g in x(r).genes:
        for key, val in FC0_dict:
            if key in str(g):
                ins[str(g)] = val
            if key not in str(g):
                notins[str(g)] = 'NO DATA'
        temp = {**notins, **ins}
        lib[r] = temp


# In[ ]:




