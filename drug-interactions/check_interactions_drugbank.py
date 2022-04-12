import pandas as pd
import requests
import json
import sys
import getopt


url_find_rxcui:str = "https://rxnav.nlm.nih.gov/REST/rxcui.json?name={0}&search=0"
url_find_interaction:str = "https://rxnav.nlm.nih.gov/REST/interaction/list.json?rxcuis={0}+{1}&source=DrugBank"


def get_rxcui(drug:str)->str:
    if ' ' in drug:
        #replace with +
        drug = drug.replace(' ','+')
    api_url = url_find_rxcui.format(drug)
    response = requests.get(api_url)
    json_response = response.json()['idGroup']['rxnormId']
    if not json_response:
        return ""
    return json_response.pop()


def check_interactions(drug_f:str, drug_s:str)->str:
    api_url = url_find_interaction.format(drug_f, drug_s)    

    response = requests.get(api_url)
    interactionsTypeGroup = response.json()['fullInteractionTypeGroup']
    if len(interactionsTypeGroup) == 2:
        try:
            fullInteractionType = interactionsTypeGroup[1]['fullInteractionType'][0]
            interactionPair = fullInteractionType['interactionPair'][0]
            interactionSeverity = interactionPair['severity']
            return interactionSeverity
        except IndexError:
            return "Exception"
    else:
        return "N/A"


if __name__ == "__main__":
    input_file = ""
    output_file = ""
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, 'i:o:', ['input', 'output'])
    try: 
        if len(opts) == 0 or len(opts) > 2:
            print ('usage: check_interactions_drugbank.py -i <input_csv> -o <output_csv>')
            sys.exit(1)
        else:
            _,input_file = opts[0]
            _,output_file = opts[1]
    except getopt.GetoptError:
        print ('usage: check_interactions.py -i <input csv> -o <output_csv>')
        sys.exit(2)
    df = pd.read_csv(input_file)
    #df = df.head(1)
    severity = []
    for _, row in df.iterrows():
        first_drug = get_rxcui(row["Name"])
        second_drug = get_rxcui(row["Interaction"])
        #first_drug = "207106"
        #second_drug = "152923"
        if not first_drug  or not second_drug:
            severity.append(' ')
            continue

        interaction = check_interactions(first_drug, second_drug)
        print("{0}\t{1}\t{2}".format(row["Name"], row["Interaction"],interaction))
        severity.append(interaction)
    df["Severity"]=severity
    df.to_csv(output_file)
