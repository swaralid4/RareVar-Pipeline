"""Build disease-associated gene list for RareVar-Pipeline."""
import json

DISEASE_GENES = {
    "BRCA1": {"disease": "Hereditary Breast/Ovarian Cancer", "omim": "113705"},
    "BRCA2": {"disease": "Hereditary Breast/Ovarian Cancer", "omim": "600185"},
    "TP53":  {"disease": "Li-Fraumeni Syndrome",           "omim": "191170"},
    "CFTR":  {"disease": "Cystic Fibrosis",                  "omim": "602421"},
    "MECP2": {"disease": "Rett Syndrome",                    "omim": "300005"},
    "MLH1":  {"disease": "Lynch Syndrome",                   "omim": "120436"},
    "MSH2":  {"disease": "Lynch Syndrome",                   "omim": "609309"},
    "PTEN":  {"disease": "Cowden Syndrome",                  "omim": "601728"},
    "RB1":   {"disease": "Retinoblastoma",                   "omim": "614041"},
    "LDLR":  {"disease": "Familial Hypercholesterolemia",   "omim": "606945"},
}

with open("data/processed/disease_genes.json", "w") as f:
    json.dump(DISEASE_GENES, f, indent=2)

print(f"Saved {len(DISEASE_GENES)} disease genes.")
