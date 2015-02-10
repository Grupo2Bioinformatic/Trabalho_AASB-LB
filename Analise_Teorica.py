from Bio import Entrez
from Bio import Medline

Entrez.email = "brveloso10@gmail.com"
handle = Entrez.egquery(term = "Neisseria gonorrhoeae")
record = Entrez.read(handle)
 
for row in record["eGQueryResult"]:
    if row["DbName"]=="pubmed":
        total = row["Count"]
        print (len(total))
handle = Entrez.esearch(db = "pubmed", term = "Neisseria gonorrhoeae", retmax=total)
record = Entrez.read(handle)
idlist = record["IdList"]

handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
records = list(Medline.parse(handle))
record_results = open('artigos.txt', 'w')
for record in records:
    texts ="Title: " + str(record.get("TI", "?")) + "\n" + "Authors: " + str(record.get("AU", "?")) + "\n" + "Source: " + str(record.get("SO", "?")) + "\n"
    record_results.write(texts)
    record_results.write("\n")
record_results.close()