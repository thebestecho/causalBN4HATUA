# Filename: process_network_array.py
import sys
import os
import pandas as pd
from causalnex.structure import StructureModel
from causalnex.network import BayesianNetwork
from causalnex.inference import InferenceEngine


def process_network(file_index, intervened_variable):
    arcs_file = f'./overall_arcs_csv/arcs_{file_index}.csv'
    imputed_file = f'./overall_imputed/imputed_{file_index}.csv'
    
    # Load arcs and create the structure model
    myarcs = pd.read_csv(arcs_file, delimiter=',', usecols=['from','to'])
    mygraph = StructureModel()
    for _, row in myarcs.iterrows():
        mygraph.add_edge(row['from'], row['to'])
    
    # Load the imputed data
    imputed = pd.read_csv(imputed_file)
    
    # Create and fit Bayesian network
    bn = BayesianNetwork(mygraph)
    bn = bn.fit_node_states(imputed)
    bn = bn.fit_cpds(imputed, method="BayesianEstimator", bayes_prior="K2")
        
    # Create an InferenceEngine
    ie = InferenceEngine(bn)
    
    # Stratified analysis for the specified variable
    if intervened_variable in imputed.columns:
        stratified = pd.crosstab(index=imputed['MDR'], columns=imputed[intervened_variable], normalize='columns')
    else:
        raise ValueError(f"{intervened_variable} is not a column in the provided data.")
    
    # Record necessary information
    results = {
        "MDR_whole": ie.query()["MDR"][1]
    }
    
    # Add stratified MDR rates to the results dictionary
    categories = imputed[intervened_variable].unique()
    for category in categories:
        if category in stratified.columns:  # Check if category exists in the stratified table
            results[f"stratified_MDR_{category}"] = stratified.loc[1, category] if 1 in stratified.index else 0
    
    # Hypothetical intervention analysis for the specified variable
    for category in categories:
        ie.reset_do(intervened_variable)
        intervention = {cat: 1.0 if cat == category else 0.0 for cat in categories}
        ie.do_intervention(intervened_variable, intervention)
        results[f"MDR_if_all_{category}"] = ie.query()["MDR"][1]
    
    ie.reset_do(intervened_variable)  # Reset intervention
    
    return results


if __name__ == "__main__":
    file_index = sys.argv[1]  # SLURM_ARRAY_TASK_ID passed as command-line argument
    intervened_variable = sys.argv[2]
    
    results = process_network(file_index, intervened_variable)

    # Create a directory named after the intervened variable if it doesn't exist
    results_dir = f'./overall_analysis_results/{intervened_variable}'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Output results to a CSV file within the directory named after the intervened variable
    output_csv_path = f'{results_dir}/results_{file_index}.csv'
    pd.DataFrame([results]).to_csv(output_csv_path)
    print(f"Analysis results for file index {file_index} have been saved to {output_csv_path}")