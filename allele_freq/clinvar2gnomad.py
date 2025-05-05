import requests
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import time
from typing import List, Dict, Optional

# Configuration
Entrez.email = "viurir@example.com"  # Replace with your email
MIM_NUMBER = "600555"
OUTPUT_FILE = "clinvar_600555_allele_frequencies.csv"

def query_clinvar_variants(mim_number: str) -> List[str]:
    """Query ClinVar for variants associated with a MIM number and return Variation IDs."""
    try:
        handle = Entrez.esearch(db="clinvar", term=f"{mim_number}[MIM]", retmax=100)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error querying ClinVar: {e}")
        return []

def fetch_variant_details(variation_id: str) -> Dict[str, Optional[str]]:
    """Fetch details for a ClinVar variant, including rsID and genomic coordinates."""
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id={variation_id}&rettype=vcv&retmode=xml"
        response = requests.get(url)
        response.raise_for_status()
        
        # Parse XML response
        root = ET.fromstring(response.text)
        variant_info = {
            "VariationID": variation_id,
            "rsID": None,
            "Chromosome": None,
            "Position": None,
            "Ref": None,
            "Alt": None,
            "Gene": None,
            "VariantName": None
        }
        
        # Extract rsID
        for xref in root.findall(".//XRef[@DB='dbSNP']"):
            variant_info["rsID"] = xref.get("ID")
        
        # Extract genomic coordinates (from SequenceLocation)
        for loc in root.findall(".//SequenceLocation[@Assembly='GRCh38']"):
            variant_info["Chromosome"] = loc.get("Chr")
            variant_info["Position"] = loc.get("start")
            variant_info["Ref"] = loc.get("referenceAlleleVCF")
            variant_info["Alt"] = loc.get("alternateAlleleVCF")
        
        # Extract gene and variant name
        for gene in root.findall(".//GeneSymbol"):
            variant_info["Gene"] = gene.text
        for name in root.findall(".//Name"):
            variant_info["VariantName"] = name.text
        
        return variant_info
    except Exception as e:
        print(f"Error fetching details for VariationID {variation_id}: {e}")
        return {"VariationID": variation_id, "rsID": None, "Chromosome": None, "Position": None, "Ref": None, "Alt": None, "Gene": None, "VariantName": None}

def query_gnomad_frequency(rsid: Optional[str], chrom: Optional[str], pos: Optional[str], ref: Optional[str], alt: Optional[str]) -> Optional[float]:
    """Query gnomAD for allele frequency using rsID or genomic coordinates."""
    try:
        if rsid:
            # Try querying gnomAD by rsID
            url = f"https://gnomad.broadinstitute.org/api/variantid/{rsid}"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                return data.get("allele_frequency", None)
        
        if chrom and pos and ref and alt:
            # Query by genomic position (e.g., 16-2056789-G-A)
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
            url = f"https://gnomad.broadinstitute.org/api/variant/{variant_id}"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                return data.get("allele_frequency", None)
        
        return None
    except Exception as e:
        print(f"Error querying gnomAD for rsID {rsid} or {chrom}-{pos}-{ref}-{alt}: {e}")
        return None

def main():
    # Step 1: Query ClinVar for variant IDs
    print(f"Querying ClinVar for MIM {MIM_NUMBER}...")
    variation_ids = query_clinvar_variants(MIM_NUMBER)
    if not variation_ids:
        print("No variants found.")
        return
    
    # Step 2: Fetch details for each variant
    print(f"Fetching details for {len(variation_ids)} variants...")
    variants = []
    for vid in variation_ids:
        variant_info = fetch_variant_details(vid)
        variants.append(variant_info)
        time.sleep(0.1)  # Avoid overwhelming NCBI servers
    
    # Step 3: Query gnomAD for allele frequencies
    print("Querying gnomAD for allele frequencies...")
    for variant in variants:
        freq = query_gnomad_frequency(
            variant["rsID"],
            variant["Chromosome"],
            variant["Position"],
            variant["Ref"],
            variant["Alt"]
        )
        variant["AlleleFrequency"] = freq if freq is not None else "Not available"
        time.sleep(0.1)  # Avoid overwhelming gnomAD servers
    
    # Step 4: Save results to a DataFrame and CSV
    df = pd.DataFrame(variants)
    df = df[["VariationID", "VariantName", "Gene", "rsID", "Chromosome", "Position", "Ref", "Alt", "AlleleFrequency"]]
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Results saved to {OUTPUT_FILE}")
    
    # Print summary
    print(f"Processed {len(variants)} variants. {sum(1 for v in variants if v['AlleleFrequency'] != 'Not available')} had allele frequencies.")

if __name__ == "__main__":
    main()