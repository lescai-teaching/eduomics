#!/bin/env python

import os
import sys
import argparse
from google import genai

def generate_clinical_scenario(variant):
    """
    Constructs a prompt that directs Google Gemini to create a realistic clinical scenario
    for a patient whose phenotype is explained by a specific genetic variant.
    The prompt includes the variant (as internal data) but instructs Gemini not to reveal it.
    """
    prompt = (
        "The following genetic data is provided as background for this case: Variant: [{}].\n\n"
        "Based on this genetic variant, please construct a detailed, biologically and clinically realistic scenario "
        "describing a hypothetical patient. The clinical features, family history, and diagnostic clues should be consistent "
        "with the patient harboring a mutation that explains their condition. However, in your final output, do not reveal or "
        "mention the actual variant details. Provide a complete narrative that implies an underlying genetic abnormality without "
        "explicitly disclosing the variant data."
    ).format(variant)

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
        description="Generate a clinical scenario based on a genetic variant using Google Gemini API."
    )
    parser.add_argument("--variant", required=True, help="The genetic variant (e.g., 'chr22-1234-A-T').")
    parser.add_argument("--output", required=True, help="Output file name for the generated scenario.")
    args = parser.parse_args()

    scenario = generate_clinical_scenario(args.variant)
    with open(args.output, "w") as f:
        f.write(scenario)
    print(f"Scenario saved to {args.output}")

if __name__ == "__main__":
    main()
