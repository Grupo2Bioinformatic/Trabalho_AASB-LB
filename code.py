from Bio import SeqIO
from Bio import Entrez


def get_sequence_file(seqstart,seqstop,file_name): #246001,468400
    Entrez.email="brveloso10@gmail.com"
    handle=Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="59800473", seq_start=seqstart, seq_stop=seqstop)
    for seq_record in SeqIO.parse(handle, "gb"):
        print(seq_record.id, seq_record.description[:50] + "...")
        print("Sequence length %i," % len(seq_record),)
        print("i% features," % len(seq_record.features),)
        print("from: %s" % seq_record.annotations["source"])

    record=SeqIO.write(seq_record, "%s.gb"%(file_name), "genbank")
    handle.close()
    return record
    
def read_sequence_file(file_name):
    record=SeqIO.read("%s.gb"%(file_name), "genbank")
    return record

def table(table_name):
    f = open('%s.txt'%(table_name), 'r')
    tabela=[]
    for line in f.readlines():
        tabela.append(line.split('\t'))
    f.close
    return tabela
    
def get_info_table(tabela,seqstart,seqstop):
    posicoes=[]
    localfim=[]
    i=1        
    while int(tabela[i][2])<int(seqstop):
        if int(tabela[i][2])>=int(seqstart):
            posicoes.append(int(tabela[i][2])-int(seqstart))
            localfim.append(int(tabela[i][3])-(int(seqstart)+1))
        i+=1
    return (posicoes,localfim)
 
def get_CDS(record):
    cds=[]
    for i in range(len(record.features)): 
        if record.features[i].type == "CDS": 
            cds.append(i)
        i=+1
    return cds

    
def get_gene(record):
    gene=[]
    for i in range(len(record.features)): 
        if (record.features[i].type == "gene"):
            gene.append(i)
        i=+1
    return gene


def  get_proteinID(record):
    CDS_proteinID =[]
    for i in record.features:
        if i.type == "CDS":
            CDS_proteinID.append(i.qualifiers["protein_id"][0])
    return CDS_proteinID

def  get_geneID(record):
    CDS_geneID =[]
    for i in record.features:
        if i.type == "CDS":
            CDS_geneID.append(i.qualifiers["db_xref"][1].strip("GeneID:"))
    return CDS_geneID
    
def  get_LocusTag(record):
    CDS_LocusTag =[]
    for i in record.features:
        if i.type == "gene":
            CDS_LocusTag.append(i.qualifiers["locus_tag"][0])
    return CDS_LocusTag
    
def  get_ECnumber(record):
    CDS_ECnumber =[]
    for i in record.features:
        if i.type == "CDS":
            if "EC_number" in i.qualifiers:
                CDS_ECnumber.append(i.qualifiers["EC_number"][0])
            else:
                CDS_ECnumber.append("None")
    return CDS_ECnumber
    
def  get_gene_name(record):
    Gene_name =[]
    for i in record.features:
        if i.type == "gene":
            if "gene" in i.qualifiers:
                Gene_name.append(i.qualifiers["gene"][0])
            else:
                Gene_name.append("None")
    return Gene_name

def get_traducao(record):
    dadostraducao=[]
    for j in record.features:
        if j.type=="CDS":
            dadostraducao.append((j.qualifiers["locus_tag"][0],j.qualifiers["translation"][0]))
            dadostraducao.append("Tamanho:")
            dadostraducao.append(len(j.qualifiers["translation"][0]))
    return(dadostraducao)
    
def get_funcao(record):
    dadosfuncoes=[]
    for i in record.features:
        if i.type == "CDS":
            if 'product' in i.qualifiers:
                dadosfuncoes.append((i.qualifiers["locus_tag"][0],i.qualifiers["product"][0]))
    return(dadosfuncoes)
    
def save(results):
    name=input("Indique o nome para o ficheiro:")
    file = open('%s.txt'%(name),'w')
    for i in range(len(results)):                                       
        file.writelines(str(results[i]))
        file.writelines("\n")
    file.close()

    
def main():
    res=input("Ja tem sequencia S/N? ")
    if res=="N" or res=="n":
        seqstart=input("Insira o inicio:")
        seqstop=input("Insira o fim:")
        file_name=input("Insira o nome para o ficheiro:")
        record=get_sequence_file(seqstart,seqstop,file_name)
    elif res=="S" or res=="s":
        file_name=input("Nome do ficheiro? ")
        seqstart=input("Insira o inicio:")
        seqstop=input("Insira o fim:")
        record =read_sequence_file(file_name)
    else:
        
        print ("Opção inválida")
    
    while True:

        op=input("""
        1 - Leitura da tabela
        2 - Leitura dos genes na tabela
        3 - Localização dos CDS na sequência
        4 - Localização dos Genes na sequência
        5 - Ver ProteinID
        6 - Ver GeneID 
        7 - Ver Locus_tag 
        8 - Ver EC_Number
        9 - Ver Nome do Gene
        10 - Ver tradução
        11 - Ver funções
        99 - Sair!
        O que pretende fazer?
        """)
        if op=='1':
            table_name=input('Indique o nome da tabela a ler:')
            tabela=table(table_name)
        elif op=='2':
            (posicoes,localfim)=get_info_table(tabela,seqstart,seqstop)
            print(posicoes)
        elif op=='3':
            record=SeqIO.read("%s.gb"%(file_name), "genbank")
            cds=get_CDS(record)
            print ("Numero de features CDS: " + str(len(cds)))
            print ("Localizacao: "+ str(cds))
        elif op=='4':
            gene=get_gene(record)
            print ("Numero de features Gene: " + str(len(gene)))
            print ("Localizacao: "+ str(gene))
        elif op=='5':
            proteinID=get_proteinID(record)
            print(proteinID)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(proteinID)
                print("Resultados gravados")
        elif op=='6': 
            geneID=get_geneID(record)
            print(geneID)            
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(geneID)
                print("Resultados gravados")
        elif op=='7':
            LocusTag=get_LocusTag(record)
            print(LocusTag)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(LocusTag)
                print("Resultados gravados")
        elif op=='8':           
            ECnumber=get_ECnumber(record)
            print(ECnumber)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(ECnumber)
                print("Resultados gravados")
        elif op=='9': 
            Gene_name=get_gene_name(record)
            print(Gene_name)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(Gene_name)
                print("Resultados gravados")
        elif op=='10': 
            trans=get_traducao(record)
            print(trans)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(trans)
                print("Resultados gravados")
        elif op=='11': 
            func=get_funcao(record)
            print(func)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(func)
                print("Resultados gravados")
        elif op=='99':
            print("\n Adeus!")
            break
        else:
            print ("Opção não valida! Tente novamente.")
    
#main
if __name__ == '__main__':
    main()
        
        
        
        
