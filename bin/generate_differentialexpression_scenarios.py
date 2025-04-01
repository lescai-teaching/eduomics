#!/bin/env python

import os
import sys
import argparse
from google import genai

def generate_clinical_scenario(genes):
    """
    Constructs a prompt that directs Google Gemini to create a realistic clinical scenario
    for a patient whose condition is explained by abnormal gene expression.
    The prompt includes the provided gene list (as internal data) but instructs Gemini not to reveal them.
    """
    # Join the gene list into a comma-separated string.
    genes_data = ", ".join(genes)

    prompt = (
        "The following genetic data is provided as background for this case: Genes: [{}].\n\n"
        "Based on these genes being abnormally expressed, please construct a detailed, biologically and clinically "
        "realistic scenario describing a hypothetical patient. The clinical narrative should include features, family history, "
        "and diagnostic clues that align with a molecular abnormality involving these genes. However, do not reveal or mention "
        "the actual gene names in your final output. The final narrative should imply an underlying gene expression abnormality without "
        "disclosing the specific gene data."
    ).format(genes_data)

    api_key = os.getenv("GOOGLE_GENAI_API_KEY")
    if not api_key:
        print("Error: Please set your GOOGLE_GENAI_API_KEY environment variable.")
        sys.exit(1)

    client = genai.Client(api_key=api_key)
    response = client.models.generate_content(
        model="gemini-2.0-flash",
        contents=prompt
    )
    return response.text.strip()

def main():
    parser = argparse.ArgumentParser(
        description="Generate a clinical scenario based on differentially expressed genes using Google Gemini API."
    )
    parser.add_argument("--genes", required=True, help="Comma-separated list of genes (e.g., 'BRCA1,TP53,EGFR').")
    parser.add_argument("--output", required=True, help="Output file name for the generated scenario.")
    args = parser.parse_args()

    # Process the comma-separated gene list.
    genes = [gene.strip() for gene in args.genes.split(",")]
    scenario = generate_clinical_scenario(genes)
    with open(args.output, "w") as f:
        f.write(scenario)
    print(f"Scenario saved to {args.output}")

if __name__ == "__main__":
    main()
