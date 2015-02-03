from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW,NCBIXML


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
    return (record)
    
def read_sequence_file(file_name):
    record=SeqIO.read("%s.gb"%(file_name), "genbank")
    return (record)
 
def get_info(record):
    CDS_proteinID =[]
    GeneID=[]
    CDS_GI =[]
    LocusTag =[]
    CDS_ECnumber =[]
    Gene_name =[]
    dadostraducao=[]
    dadosfuncoes=[]
    CDS_note =[]
    size=[]
    for i in record.features: 
        if i.type == "CDS": 
            CDS_proteinID.append(i.qualifiers["protein_id"][0])
            CDS_GI.append(i.qualifiers["db_xref"][0].strip("GI:"))
            if "EC_number" in i.qualifiers:
                CDS_ECnumber.append(i.qualifiers["EC_number"][0])
            else:
                CDS_ECnumber.append("None")
            dadostraducao.append((i.qualifiers["locus_tag"][0],i.qualifiers["translation"][0]))
            dadostraducao.append("Tamanho:")
            dadostraducao.append(len(i.qualifiers["translation"][0]))
            dadosfuncoes.append((i.qualifiers["locus_tag"][0],i.qualifiers["product"][0]))
            size.append(len(i.qualifiers["translation"][0]))
            if "note" in i.qualifiers:
                CDS_note.append((i.qualifiers["locus_tag"][0],i.qualifiers["note"][0]))
            else:
                CDS_note.append("None")
        if i.type == "gene":
            GeneID.append(i.qualifiers["db_xref"][0].strip("GeneID:"))
            LocusTag.append(i.qualifiers["locus_tag"][0])
            if "gene" in i.qualifiers:
                Gene_name.append(i.qualifiers["gene"][0])
            else:
                Gene_name.append("None")
    return (CDS_proteinID,GeneID,CDS_GI,LocusTag,CDS_ECnumber,Gene_name,dadostraducao,dadosfuncoes,CDS_note,size)

def save(results):
    name=input("Indique o nome para o ficheiro:")
    file = open('%s.txt'%(name),'w')
    for i in range(len(results)):                                       
        file.writelines(str(results[i]))
        file.writelines("\n")
    file.close()    

def blast(GI,filename):
    result_handle = NCBIWWW.qblast("blastp", "swissprot", GI)
    save_file = open("%s.xml"%(filename), "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

def analisa(filename):
    result_handle = open("%s.xml"%(filename))
    E_VALUE_THRESH = 0.05
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print("****Alignment****")
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
            else:
                print("Os resultados encontrados não são considerados bons")
    result_handle.close()
    
def table(table_name):
    f = open('%s.txt'%(table_name), 'r')
    tabela=[]
    for line in f.readlines():
        tabela.append(line.split('\t'))
    f.close
    return (tabela)
    
def compare(table,seqstart,seqstop):
    p_ID_Table=[]
    for line in range(1,len(table)):
        if (table[line][2]>=seqstart and table[line][3]<=seqstop):
            p_ID_Table.append(table[line][8]) 
    return (p_ID_Table)    

def valida(Tab_ID,CDS_ID):
    valid=[]
    for p in range(len(CDS_ID)):
        for k in range(len(Tab_ID)):
            if CDS_ID[p]==Tab_ID[k]:
                valid.append((CDS_ID[p],"=Válido"))
    return valid
    
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
    record=SeqIO.read("%s.gb"%(file_name), "genbank")
    (CDS_proteinID,GeneID,CDS_GI,LocusTag,CDS_ECnumber,Gene_name,dadostraducao,dadosfuncoes,CDS_note,size)=get_info(record)
    while True:

        op=input("""
        1 - Ver ProteinID
        2 - Ver GeneID 
        3 - Ver GI
        4 - Ver Locus_tag 
        5 - Ver EC_Number
        6 - Ver Nome do Gene
        7 - Ver tradução
        8 - Ver funções
        9 - Ver descrições
        10 - Realizar Blast
        11 - Leitura da tabela
        12 - Validação
        99 - Sair!
        O que pretende fazer?
        """)
        if op=='1':
            print(CDS_proteinID)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(CDS_proteinID)
                print("Resultados gravados")

        elif op=='2':
            print(GeneID)            
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(GeneID)
                print("Resultados gravados")
        elif op=='3':
            print(CDS_GI)            
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(CDS_GI)
                print("Resultados gravados")
        elif op=='4':
            print(LocusTag)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(LocusTag)
                print("Resultados gravados")
        elif op=='5':
            print(CDS_ECnumber)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(CDS_ECnumber)
                print("Resultados gravados")
        elif op=='6': 
            print(Gene_name)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(Gene_name)
                print("Resultados gravados")
        elif op=='7':
            print(dadostraducao)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(dadostraducao)
                print("Resultados gravados")
        elif op=='8':           
            print(dadosfuncoes)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(dadosfuncoes)
                print("Resultados gravados")
        elif op=='9': 
            print(CDS_note)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(CDS_note)
                print("Resultados gravados")
        elif op=='10': 
            res=input("Indique o GI do gene:")
            file=input("Indique o nome para o ficheiro:")
            blast(res,file)
            esc=input("Deseja analisar o blast?(S/N):")
            if esc=="S" or esc=="s":
                analisa(file)
        elif op=='11':
            table_name=input('Indique o nome da tabela a ler:')
            tabela=table(table_name)
            print(tabela)
        elif op=='12':
            tab_ID=compare(tabela,seqstart,seqstop)
            valid=valida(tab_ID,CDS_proteinID)
            print(valid)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(valid)
                print("Resultados gravados")
        elif op=='13': 
            print(size)
            res=input("Deseja gravar os resultados num ficheiro?(S/N):")
            if res=="s" or res=="S":
                save(size)
                print("Resultados gravados")
        elif op=='99':
            print("\n Adeus!")
            break
        else:
            print ("Opção não valida! Tente novamente.")
    
#main
if __name__ == '__main__':
    main()
        
