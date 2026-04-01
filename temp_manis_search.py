import urllib.parse, urllib.request, xml.etree.ElementTree as ET

species_terms = [
    ("Manis javanica", "javanica"),
    ("Manis pentadactyla", "pentadactyla"),
]

base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

for sp, tag in species_terms:
    term = f'"{sp}"[Organism] AND ("rna-seq"[All Fields] OR transcriptome[All Fields] OR expression[All Fields])'
    url = f"{base}/esearch.fcgi?db=sra&term={urllib.parse.quote(term)}&retmax=30"
    with urllib.request.urlopen(url, timeout=60) as r:
        xml = r.read()
    root = ET.fromstring(xml)
    count = root.findtext("Count")
    ids = [e.text for e in root.findall("IdList/Id") if e.text]
    print(f"=== {sp} ===")
    print(f"Count: {count}")
    print("IDs:", ",".join(ids[:10]) if ids else "none")

