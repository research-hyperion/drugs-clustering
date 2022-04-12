#!/usr/bin/env python
# coding: utf-8

# In[69]:


import csv
import xml.etree.ElementTree as ET
import time
import io
from tqdm import tqdm
import pandas as pd
import requests
import json
import logging


# In[70]:


URL_API_INTERACTIONS:str = "https://rxnav.nlm.nih.gov/REST/interaction/list.json"


# In[71]:


URL_API_RXCUI:str =  "https://rxnav.nlm.nih.gov/REST/rxcui.json"


# In[72]:


logging.basicConfig(filename='approved_drugs_log.log', level=logging.DEBUG)


# In[73]:


ERRORS_LIST = [
    (0, "SUCCESS"),
    (1, "Error in parsing RXCUI response JSON for drug {0}."),
    (2, "JSON RXCUI for drug {0} is empty."),
    (3, "Index Error while checking the interactions between {0} and {1}."),
    (4, "N/A"),
    (5, "Error in parsing interactions between {0}-{1} response JSON.")
]


# In[74]:


def get_rxcui(drugid:str):
    response = requests.get(URL_API_RXCUI,params={"idtype":"Drugbank","id":drugid})
    #print(response.json())
    try:
        json_response = response.json()['idGroup']['rxnormId']
        if not json_response:
            error_code = ERRORS_LIST[1]
            error_message = error_code[1].format(drugid)
            return error_code[0],error_message
        return 0,json_response.pop()
    except:
        error_code = ERRORS_LIST[2]
        error_message = error_code[1].format(drugid)
        return error_code[0],error_message
        


# In[75]:


def check_interactions(drug_f:str, drug_s:str):
    response = requests.get(URL_API_INTERACTIONS, params={"rxcuis":drug_f+"+"+drug_s,"source":"DrugBank"})
    try:
        interactionsTypeGroup = response.json()['fullInteractionTypeGroup']
        if len(interactionsTypeGroup) == 2:
            try:
                fullInteractionType = interactionsTypeGroup[1]['fullInteractionType'][0]
                interactionPair = fullInteractionType['interactionPair'][0]
                interactionSeverity = interactionPair['severity']
                return 0,interactionSeverity
            except IndexError:
                error_code = ERRORS_LIST[3]
                error_message = error_code[1].format(drug_f,drug_s)
                return error_code[0],error_message
        else:
            error_code = ERRORS_LIST[4]
            error_message = error_code[1].format(drug_f,drug_s)
            return error_code[0],error_message
    except:
            error_code = ERRORS_LIST[5]
            error_message = error_code[1].format(drug_f,drug_s)
            return error_code[0],error_message
        


# In[76]:


DRUGBANK_DATABASES:list() = [
    (0,"DB 3.0","database/3.0/drugbank.xml",2013),
    (1,"DB 4.1","database/4.1/drugbank.xml",2014), #not sure about year
    (2,"DB 4.2","database/4.2/drugbank.xml",2015), #not sure about year
    (3,"DB 4.3","database/4.3/drugbank.xml",2015),
    (4,"DB 4.5","database/4.5.0/drugbank.xml",2016),
    (5," DB 5.0","database/5.0.0/drugbank.xml",2016),
    (6,"DB 5.0.1","database/5.0.1/drugbank.xml",2016),
    (7,"DB 5.0.2","database/5.0.2/drugbank.xml",2016),
    (8,"DB 5.0.3","database/5.0.3/drugbank.xml",2016),
    (9,"DB 5.0.4","database/5.0.4/full_database.xml",2017),
    (10,"DB 5.0.5","database/5.0.5/full_database.xml",2017),
    (11,"DB 5.0.6","database/5.0.6/full_database.xml",2017),
    (12,"DB 5.0.7","database/5.0.7/full_database.xml",2017),
    (13,"DB 5.0.8","database/5.0.8/full_database.xml",2017),
    (14,"DB 5.0.9","database/5.0.9/full_database.xml",2017), 
    (15,"DB 5.0.10","database/5.0.10/full_database.xml",2017),
    (16,"DB 5.0.11","database/5.0.11/full_database.xml",2017),
    (17,"DB 5.1.0","database/5.1.0/full_database.xml",2018),
    (18,"DB 5.1.1","database/5.1.1/full_database.xml",2018),
    (19,"DB 5.1.2","database/5.1.2/full_database.xml",2018),
    (20,"DB 5.1.3","database/5.1.3/full_database.xml",2019),
    (21,"DB 5.1.4","database/5.1.4/full_database.xml",2019),
    (22,"DB 5.1.5","database/5.1.5/full_database.xml",2020),
    (23,"DB 5.1.6","database/5.1.6/full_database.xml",2020),
    (24,"DB 5.1.7","database/5.1.7/full_database.xml",2020),
    (25,"DB 5.1.8","database/5.1.8/full_database.xml",2021)
]


# In[77]:


#for db_entry in DRUGBANK_DATABASES:
for i in tqdm(range(len(DRUGBANK_DATABASES))):
    db_entry = DRUGBANK_DATABASES[i]
    database_name = db_entry[1].replace(".","_")
    database_name_trim = database_name.replace(" ","")
    input_csv_file = "database/approved_only/"+database_name_trim+"_approved_drugs.csv"
    output_csv_file = "database/interactions/"+database_name_trim+"_interaction_strength.csv"
    df = pd.read_csv(input_csv_file)
    #create header for dataframe
    out_df =  pd.DataFrame(columns=['DrugDB_Id', "InteractionDB_Id","Interaction_Strength"])
    #iterate on each row from csv file
    for index, row in df.iterrows():
        drug_db_id = row['DB_id']
        drug_interactions = row['DB_interactions']
        if not pd.isna(drug_interactions):
            drug_interactions_list = drug_interactions.split(';')
            #get rxcui for target drug
            drug_rxcui = get_rxcui(drug_db_id)
            #if the rxcui for the target drug is not found, just to the next row
            if drug_rxcui[0]!=0:
                logging.warning(drug_rxcui[1])
                continue
            for drug in drug_interactions_list:
                interaction_drug_rxcui = get_rxcui(drug)
                if interaction_drug_rxcui[0]!=0:
                    logging.warning(interaction_drug_rxcui[1])
                    continue
                #Check for interaction!
                interaction_strength = check_interactions(drug_rxcui[1], interaction_drug_rxcui[1])
                if interaction_strength[0]!=0:
                    logging.error(interaction_drug_rxcui[1])
                    continue
                out_df = df.append({'DrugDB_Id': drug_db_id,'InteractionDB_Id':drug,'Interaction_Strength':interaction_strength[1]}, ignore_index=True)
        else:
            logging.info("No interactions found")
    out_df.to_csv(output_csv_file, index=False)     
            


# In[ ]:




