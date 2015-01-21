from Bio import SeqIO
from Bio import Entrez

#Analise da Sequencia e das features presentes no NCBI

#Primeiro Ponto
    
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

#Segundo e Terceiro Ponto
#ProteinTable864_169534
def table(table_name):
    f = open('%s.txt'%(table_name), 'r')
    tab=f.readlines()
    tabela=[]
    for line in tab:
        tabela.append(tab[line].encode(utf-8).split('\t'))
    return tabela
    
def get_info_table(tabela,seqstart,seqstop):
    posicoes=[]
    localfim=[]
    i=0        
    while tabela[i][2]<seqstop:
        if tabela[i][2]>=seqstart:
            posicoes.append(tabela[i][2]-seqstart)
            localfim.append(tabela[i][3]-(seqstart+1))
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

    
##############informacao para a validacao#################
#obtem geneID e protein do ficheiro genbank para a sequencia a trabalhar

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
        

##############validacao#############      
#na alteracao feita foram utilizadas as posicoes determinadas anteriormente
#em vez de modificar o ficheiro com a tabela
  
def validation(tabela,posicoes,CDS_proteinID,CDS_geneID):
    valid=[]    
    valido=True
    k=0
    for k in range(len(posicoes)):
        j=posicoes[k] 
        if tabela[j][5] != CDS_geneID[k]or tabela[j][8] != CDS_proteinID[k]:
            valid.append(False)
        else:
            valid.append(valido)
    return valid

 
 
 
##########print,teste e validacao#############
def analise(cds,gene,posicoes, valid):
    teste=True    
    if len(posicoes)==len(cds):
        equal=True
    for i in valid:
        if valid[i]==False:
            teste=False           
    return(equal,teste)


    
    
    
    
    

    
    
    
    
    
#Terceiro Pont
    
def main():
    res=raw_input("Ja tem sequencia S/N? ")
    if res=="N" or res=="n":
        seqstart=input("Insira o inicio:")
        seqstop=input("Insira o fim:")
        file_name=raw_input("Insira o nome para o ficheiro:")
        record=get_sequence_file(seqstart,seqstop,file_name)
    elif res=="S" or res=="s":
        file_name=raw_input("Nome do ficheiro? ")
        seqstart=input("Insira o inicio:")
        seqstop=input("Insira o fim:")
        record =read_sequence_file(file_name)
    else:
        
        print "OpÃ§Ã£o invÃ¡lida"
    
    while True:

        op=input("""
        1 - Leitura da tabela
        2 - Leitura dos genes na tabela
        3 - LocalizaÃ§Ã£o dos CDS na sequÃªncia
        4 - LocalizaÃ§Ã£o dos Genes na sequÃªncia
        5 - Ver os ProteinID
        6 - Ver os GeneID 
        7 - ValidaÃ§Ã£o features
        8 - Ver EC_Number 
        8 - Ver Locus_tag
        9 - Ver TranslaÃ§Ã£o
        10 - Blast
        
        O que pretende fazer?
        """)
        if op==1:
            table_name=raw_input('Indique o nome da tabela a ler:')
            tabela=table(table_name)
            print table
        elif op==2:
            (posicoes,localfim)=get_info_table(tabela,seqstart,seqstop)
            print posicoes
            print localfim
        elif op==3:
            record=SeqIO.read("%s.gb"%(file_name), "genbank")
            cds=get_CDS(record)
            print "Numero de features CDS: " + str(len(cds))
            print "Localizacao: "+ str(cds)
        elif op==4:
            gene=get_gene(record)
            print ("Numero de features Gene: " + str(len(gene)))
            print ("Localizacao: "+ str(gene))
        elif op==5:
            proteinID=get_proteinID(record)
            print(proteinID)           
        elif op==6: 
            geneID=get_geneID(record)
            print(geneID)
        elif op==7:
            valid=validation(tabela,posicoes,proteinID,geneID)
            print valid
            (equal,teste)=analise(cds,gene,posicoes, valid)
            if equal==True:               
                print ("O numero de coding sequences corresponde as da informacao dada pela tabela")
            if teste==False:         
                print ("Existem features invalidas!")
            else:
                print ("As features obtidas sao validas!")
        elif op==9:
            print("\n Adeus!")
            break
        else:
            print "OpÃ§Ã£o nÃ£o valida! Tente novamente."
    
#main
if __name__ == '__main__':
    main()
