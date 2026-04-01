import csv
import io
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

species_terms = [
    "Manis javanica",
    "Manis pentadactyla",
]
base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
out_dir = Path("Expression/Pangolin")
out_dir.mkdir(parents=True, exist_ok=True)

for sp in species_terms:
    term = f'"{sp}"[Organism] AND ("rna-seq"[All Fields] OR transcriptome[All Fields] OR expression[All Fields])'
    search_url = f"{base}/esearch.fcgi?db=sra&term={urllib.parse.quote(term)}&retmax=1000"
    with urllib.request.urlopen(search_url, timeout=60) as r:
        root = ET.fromstring(r.read())
    ids = [e.text for e in root.findall("IdList/Id") if e.text]

    if not ids:
        print(f"No IDs for {sp}")
        continue

    id_blob = ",".join(ids)
    runinfo_url = f"{base}/efetch.fcgi?db=sra&id={id_blob}&rettype=runinfo&retmode=text"
    with urllib.request.urlopen(runinfo_url, timeout=120) as r:
        runinfo_text = r.read().decode("utf-8", errors="replace")

    safe = sp.replace(" ", "_")
    out_csv = out_dir / f"{safe}_sra_runinfo.csv"
    out_csv.write_text(runinfo_text, encoding="utf-8")

    # quick summary stats
    reader = csv.DictReader(io.StringIO(runinfo_text))
    rows = list(reader)
    runs = [x.get("Run", "") for x in rows if x.get("Run")]
    studies = sorted({x.get("SRAStudy", "") for x in rows if x.get("SRAStudy")})
    assays = sorted({x.get("LibraryStrategy", "") for x in rows if x.get("LibraryStrategy")})

    print(f"=== {sp} ===")
    print(f"SRA internal IDs: {len(ids)}")
    print(f"Run rows: {len(rows)}")
    print(f"Distinct studies: {len(studies)}")
    print(f"Studies (first 10): {', '.join(studies[:10])}")
    print(f"Assays: {', '.join(a for a in assays if a)}")
    print(f"Output: {out_csv}")

